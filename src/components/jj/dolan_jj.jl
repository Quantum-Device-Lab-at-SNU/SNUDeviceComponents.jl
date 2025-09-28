import ...Utils: generate_undercut
using DeviceLayout
import DeviceLayout: Meta, Coordinate

"""
    layerpixels!{T}(c::Cell{T}, layers::AbstractMatrix{Int}, pixsize)
Given `layers`, a matrix of `Int`, make a bitmap of `Rectangle` where the GDS-II layer
corresponds to the number in the matrix. If the number is less than one, don't
write the rectangle. All the rectangles get rendered into cell `c`. The rectangles are all
in the first quadrant of the cell.
"""
function layerpixels!(c::Cell{T}, layers::AbstractMatrix{Int}, pixsize::S) where {T, S<:Coordinate}
    s = size(layers)
    for i in 1:s[1], j in 1:s[2]
        layers[i,j] < 0 && continue
        r = Rectangle(pixsize, pixsize)
        r += Point((j-1) * pixsize, (s[1]-i)*pixsize)
        render!(c, r, GDSMeta(layers[i,j]))
    end
    c
end


function bandage_hatching!(
    c::CoordinateSystem{T}, pix, pixpattern, ruc, uuc, undercut_meta::DeviceLayout.Meta;
    α = 0°
) where {T <: Coordinate}
    # create pixpattern
    layerpixels!(c, pixpattern, pix)

    # undercut regions to subtract from the pixpattern
    r1 = Rectangle(pix, uuc) + Point(zero(T), 2 * pix - uuc)
    r2 = Rectangle(ruc, pix) + Point(2 * pix - ruc, zero(T))
    r3 = Rectangle(ruc, uuc) + Point(2 * pix - ruc, 2 * pix - uuc)

    # replace polygons with the subtracted one
    plgs = to_polygons(difference2d(c.elements, Polygon{T}[r1,r2,r3]))
    c.elements = plgs
    c.element_metadata = fill(element_metadata(c)[1], length(c.elements))

    # rotate and create undercuts
    uc_plgs = generate_undercut(plgs, zero(T), ruc, uuc, zero(T); α=α)
    render!(c, uc_plgs, undercut_meta)
    return c
end


function heal_jj_uc!(c::CoordinateSystem{T}, jj_meta::Meta, uc_meta::Meta) where {T<:Coordinate}
    flatten!(c)

    # Flatten and take union of full dose elements
    jj_inds = findall(c.element_metadata .== jj_meta)
    jj_polys = to_polygons(union2d(view(c.elements, jj_inds)))

    deleteat!(c.elements, jj_inds)
    deleteat!(c.element_metadata, jj_inds)
    append!(c.elements, jj_polys)
    append!(c.element_metadata, fill(jj_meta, length(jj_polys)))

    # subtract full dose elements from uc elements 
    uc_inds = findall(c.element_metadata .== uc_meta)
    uc_polys = to_polygons(difference2d(union2d(view(c.elements, uc_inds)), jj_polys))

    deleteat!(c.elements, uc_inds)
    deleteat!(c.element_metadata, uc_inds)
    append!(c.elements, uc_polys)
    append!(c.element_metadata, fill(uc_meta, length(uc_polys)))

    return c
end


"""
    jjPadTopSide!(c, pix, uc, uuc, ruc, ext, padh, bandaid_buffer, jj_meta, uc_meta; α=0°)
- `c`: Cell for top-side JJ pad (should be empty)
- `pix`: size of the "pixel" of the bandaid hatching
- `uc`: safety margin (additional undercut beyond what is needed from geometry)
- `uuc`: Upper undercut
- `ruc`: Right undercut
- `ext`: Width of rectangular landing pad
- `padh`: Height of rectangular landing pad
- `bandaid_buffer`: How far the bandaid should extend beyond the patch (ea. side)
- `α`: Relative angle of the josephson junction with respect to qubit direction.
"""
function jjPadTopSide!(
    c::Cell{T}, pix, uc, uuc, ruc, ext, padh, bandaid_buffer,
    jj_meta::DeviceLayout.Meta, uc_meta::DeviceLayout.Meta,
    bandage_meta::DeviceLayout.Meta; α=0°
) where {T<:Coordinate}

    # generate pixcell for creating hatching pattern
    pixcell = Cell(uniquename("pixcell"), nm)
    bandage_hatching!(
        pixcell, pix, [jj_meta.layer jj_meta.layer; -1 jj_meta.layer],
        ruc, uuc, uc_meta; α = α
    )

    # Part 1. Big contact pad
    pad_rect = Rectangle(ext, padh)
    render!(c, pad_rect, GDSMeta(JJ_LAYER))
    #   generate undercut for the big contact pad (rotated by angle α)
    ucp = generate_undercut(pad_rect, zero(T), ruc, uuc, zero(T); α=α)
    render!(c, ucp, GDSMeta(UC_LAYER))

    # Part 2. Hatching Pattern
    # width and height of layer pixel
    lp_bound = bounds(pixcell)
    lpw, lph = width(lp_bound), height(lp_bound)

    ya, yb, yc = promote(padh, 2 * pix, padh + lph)
    yrange = ustrip(ya):ustrip(yb):ustrip(yc)

    x0 = lp_bound.ll[1]
    xa, xb, xc = promote(-x0, 2 * pix, ext - lpw - x0)
    xrange = ustrip(xa):ustrip(xb):ustrip(xc)

    nc, nr = length(xrange), length(yrange)
    push!(c.refs,
        CellArray(pixcell, Point(xrange[1] * unit(xa), yrange[1] * unit(ya));
        nc = nc, nr = nr,
        dc = Point(step(xrange) * unit(xa), zero(T)),
        dr = Point(zero(T), step(yrange) * unit(ya))))

    # Part 3. Bandaid
    patch = Rectangle(
        (nc - 1) * 2 * pix + lpw + 2 * bandaid_buffer,
        (nr - 1) * 2 * pix + lph + 2 * bandaid_buffer
    )

    render!(c,
        Rounded(patch, 0.3um) + Point(-bandaid_buffer, padh - bandaid_buffer),
        bandage_meta
    )

    heal_jj_uc!(c, jj_meta, uc_meta)

    # Center in x (up to undercut and bandaid); bottom edge of elements will sit at y = 0
    c -= Point(ext/2, zero(T))
    return c
end

function loc_idx(loc::Symbol)
    if loc == :left
        return -1
    elseif loc == :right
        return +1
    else
        throw(ArgumentError("`loc` must be either `:left` or `:right`"))
    end
end

"""
    jjPadBotSide!(pix, uc, uuc, ruc, ext, padh, bandaid_buffer,
        armlength, armwidth, lr, jj_meta, uc_meta; α=0)
- `c`: Cell for top-side JJ pad (should be empty)
- `pix`: size of the "pixel" of the bandaid hatching
- `uc`: safety margin (additional undercut beyond what is needed from geometry)
- `uuc`: Upper undercut
- `ruc`: Right undercut
- `ext`: Width of rectangular landing pad
- `padh`: Height of rectangular landing pad
- `bandaid_buffer`: How far the bandaid should extend beyond the patch (ea. side)
- `armlength`: Length of the arm connecting the landing pad and the hatching pattern
- `armwidth`: Width of the arm connecting the landing pad and the hatching pattern
- 'loc': `:left` for left jj, `:right` for right jj
- `jj_direction`: Orientation of the josephson junction with respect to qubit direction.
    (0: 0°, 1: 90°, 2: 180°, and 3: 270°)
"""
function jjPadBotSide!(
    c::Cell{T}, pix, uc, uuc, ruc, ext, padh, bandaid_buffer,
    armlength, armwidth, loc::Symbol, jj_meta::DeviceLayout.Meta,
    uc_meta::DeviceLayout.Meta, bandage_meta::DeviceLayout.Meta; α=0
) where {T <: Coordinate}

    lr = loc_idx(loc)

    # generate pixcell for creating hatching pattern
    pixcell = Cell(uniquename("pixcell"), nm)
    bandage_hatching!(
        pixcell, pix, [jj_meta.layer jj_meta.layer; -1 jj_meta.layer],
        ruc, uuc, uc_meta; α = α
    )

    # Part 1. Big contact pad
    pad_rect = Rectangle(ext, padh)
    render!(c, pad_rect, jj_meta)
    #   generate undercut for the big contact pad (rotated by angle α)
    ucp = generate_undercut(pad_rect, zero(T), ruc, uuc, zero(T); α=α)
    render!(c, ucp, uc_meta)

    # Part 2. Arm
    arm_rect = (
        Rectangle(armlength, armwidth)
        + ifelse(lr == -1, Point(-armlength, zero(T)), Point(ext, zero(T)))
    )
    render!(c, arm_rect, jj_meta)
    #   generate undercut for the arm (rotated by angle α)
    ucp = generate_undercut(arm_rect, zero(T), ruc, uuc, zero(T); α=α)
    render!(c, ucp, uc_meta)

    # Part 3. Hatching Pattern
    lp_bound = bounds(pixcell)
    lpw, lph = width(lp_bound), height(lp_bound)  # width and height of layer pixel

    ya, yb, yc = promote(- 2 * pix, -2 * pix, -lph - 2 * pix)
    yrange = ustrip(ya):ustrip(yb):ustrip(yc)

    x0 = lp_bound.ll[1]
    xa, xb, xc = promote(-x0, 2 * pix, ext - lpw - x0)
    xrange = ustrip(xa):ustrip(xb):ustrip(xc)
    right = xrange[end] * unit(xa) + lpw
    left = xrange[1] * unit(xa)

    # position to add patch to
    patch_adj_point = ifelse(
        lr == -1,
        Point(-armlength, zero(T)),
        Point(ext + armlength -lpw - 2 * pix * (length(xrange) - 1), zero(T))
    )

    # pixcell
    push!(c.refs,
        CellArray(pixcell, Point(xrange[1] * unit(xa), yrange[1] * unit(ya)) +
            patch_adj_point;
            dc = Point(step(xrange) * unit(xa), zero(T)),
            dr = Point(zero(T), step(yrange) * unit(ya)),
            nc = length(xrange), nr = length(yrange)))

    # Part 4. Bandaid
    patch = Rectangle(
        (length(xrange) - 1) * 2 * pix + lpw + 2 * bandaid_buffer,
        (length(yrange) - 1) * 2 * pix + lph + 2 * bandaid_buffer
    )

    render!(c,
        Rounded(patch, 0.4um) + Point(-bandaid_buffer, -height(patch) + bandaid_buffer) + patch_adj_point,
        bandage_meta
    )

    # top edge of elements will sit at y = uuc
    c -= Point(ext/2, padh)

    # Flatten and take union of full dose elements, subtract jj layer elements from uc
    heal_jj_uc!(c, jj_meta, uc_meta)

    return c
end

"""
    dolan_jj!{T}(c::Cell{T}, m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
        jj_meta::Meta, uc_meta::Meta)
Explanation of parameters:

```                <-w3->
                   |████|  ◬
                   |████|  |
                   |████|  |
                   |█r3█|  l3
                   |████|  |
                   |████|  |
                   |████|  |
 ◬   ––––––––––––––+████|  ▿
 w2  ████████r2█████████|
 ▿ ◬ –––––––––––––––––––+
   b          <–l2–>
   |   <–-l4+w1–>
   ▿    ________
       |█r4█████|  ◬    | w4
          |██|     |
          |██|     |
     <-m->|r1|     l1
          |██|     |
          |██|     |
          |██|     ▿
          <w1>
          X origin
```

Other parameters:

- `c`: Cell for junctions (should be empty)
- `uc`: Amount of undercut to extend uniformly around the features
in the diagram. Note that no undercut extends beyond the top & bottom of figure.
"""
function dolan_jj!(c::Cell{T}, m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
        jj_meta::DeviceLayout.Meta, uc_meta::DeviceLayout.Meta) where {T <: Coordinate}

    # start drawing JJ
    r1 = Rectangle(Point(zero(T), zero(T)), Point(w1, l1 - w4))          # bottom vertical rectangle
    r2 = Rectangle(Point(-m , l1 + b), Point(w1 + l2 + w3,l1 + b + w2))  # top horizontal rectangle
    r3 = Rectangle(Point(w1 + l2, l1 + b + w2), Point(w1 + l2 + w3, l1 + b + w2 + l3)) # top vertical rectangle
    r4 = Rectangle(Point(-l4 / 2, l1 - w4), Point(w1 + l4 / 2, l1))      # bottom horizontal rectangle

    p14 = union2d(r1, r4) # join r1, r4 into one polygon: bottom fully exposed region
    p23 = union2d(r2, r3) # join r2, r3 into one polygon: top fully exposed region

    # Generate undercut
    ucb2 = uc + t1 * tan(ϕ)
    u1 = offset(Rectangle(Point(zero(T), zero(T)), Point(w1,l1+b)), uc)[1]
    u23 = offset(p23, uc)[1]
    uca = Rectangle(Point(- m - uc, l1 + b + w2), Point(w1 + l2, uc + l1 + b + w2 + t1 * tan(Θ)))
    ucb1 = Rectangle(Point(w1, zero(T)), Point(w1 + uc + t1 * tan(ϕ), l1 + b))
    ucb2 = Rectangle(Point(w1 + l2 + w3, l1 + b - uc), Point(w1 + l2 + w3 + uc + t1 * tan(ϕ), l1 + b + w2 + l3))
    u = union2d(u1, u23)
    u = union2d(u, uca)
    u = union2d(u, ucb1)
    u = union2d(u, ucb2)

    # Remove undercut sticking out top and bottom
    u = intersect2d(u, Rectangle(Point(T(-1000000), zero(T)), Point(T(1000000),l1+b+w2+l3)))
    u = difference2d(u, p14)
    u = difference2d(u, p23)

    # Rendering
    render!(c, u, uc_meta)
    render!(c, p14, jj_meta)
    render!(c, p23, jj_meta)

    c -= 
    return c
end

"""
    convenience function for creating a JJ with pad
"""
function dolan_jj_with_pad!(
    c::Cell{T}, m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
    pix, uuc, ruc, ext, buffer, qubit_bottom_gap, loc::Symbol,
    jj_meta::DeviceLayout.Meta, uc_meta::DeviceLayout.Meta,
    bandage_meta::DeviceLayout.Meta; α=0
) where {T<:Coordinate}

    jj_height = l3 + w2 + b + l1
    gap_between_leads_for_jjs = jj_height + 2 * buffer

    # integer to represent the relative direction of JJ with respect to qubit
    jj_direction = fld(mod(α + π/4, 2π), π/2)

    # integer to represent the location (left or right) of the JJ with respect to qubit
    lr = loc_idx(loc)

    # height of the pad
    padh = ifelse(
        jj_direction % 2 == 0,
        (qubit_bottom_gap) / 2 + 2 * buffer,
        (qubit_bottom_gap + gap_between_leads_for_jjs) / 2 + 2 * buffer
    )

    armlength = ifelse(
        jj_direction % 2 ==0,
        9um, 9um - (lead_width + gap_between_leads_for_jjs)
    )
    armwidth = pix

    # create jj pads
    jjpadtop = Cell(uniquename("jjpadtop"), nm)
    jjpadtop = jjPadTopSide!(
        jjpadtop, pix, uc, uuc, ruc, ext, padh, buffer,
        jj_meta, uc_meta, bandage_meta; α=α
    )
    jjpadbot = Cell(uniquename("jjpadbot"), nm)
    jjpadbot = jjPadBotSide!(
        jjpadbot, pix, uc, uuc, ruc, ext, padh, buffer, armlength, armwidth, loc,
        jj_meta, uc_meta, bandage_meta; α=α
    )
    # Create Josephson Junctions
    jj = Cell(uniquename("jj"), nm)
    dolan_jj!(
        jj, m, b, w1, w2, w3, w4, l1 + w1 * abs(sin(α)),
        l2, l3 + w3 * abs(sin(α)), l4, uc, t1, Θ, ϕ,
        GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER)
    )

    jj_bound = bounds(jj)
    jw, jh = width(jj_bound), height(jj_bound)


    toppad_loc = ifelse(
        jj_direction % 2 == 0,
        Point(-lr * lead_width/2, jh - (2 * uc + ruc) * abs(sin(α))),
        Point(-lr * (jh + ext) / 2, (jw/2 - lead_width/2 - buffer))
    )

    botpad_loc = ifelse(
        jj_direction % 2 == 0,
        Point(lr * lead_width/2, zero(T)),
        Point(lr * (ext + jh) / 2, (jw/2 + lead_width/2 + buffer))
    )

    # jj_rot_offset = Point((jh * abs(sin(α))) / 2, -w1 * abs(sin(α)))
    jj_rot_offset = Point((w1 + l2 + w3 - m) / 2, jj_height / 2)
    # add jj and pads to the reference 
    # addref!(c, jj, jj_rot_offset; rot = α)
    addref!(c, jj - jj_rot_offset, jj_rot_offset; rot = α)
    addref!(c, jjpadtop, toppad_loc)
    addref!(c, jjpadbot, botpad_loc)

    heal_jj_uc!(c, jj_meta, uc_meta)

end