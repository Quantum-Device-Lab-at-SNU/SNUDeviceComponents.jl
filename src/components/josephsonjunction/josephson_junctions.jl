# work in progress

using DeviceLayout
import DeviceLayout: Coordinate, CoordinateSystem
import DeviceLayout: Meta

"""
    layerpixels!{T}(c::CoordinateSystem{T}, layers::AbstractMatrix{Int}, pixsize)
Given `layers`, a matrix of `Int`, make a bitmap of `Rectangle` where the GDS-II layer
corresponds to the number in the matrix. If the number is less than one, don't
write the rectangle. All the rectangles get rendered into cell `c`. The rectangles are all
in the first quadrant of the cell.
"""
function layerpixels!(c::CoordinateSystem{T}, layers::AbstractMatrix{Int}, pixsize::S) where {T, S<:Coordinate}
    s = size(layers)
    for i in 1:s[1], j in 1:s[2]
        layers[i,j] < 0 && continue
        r = Rectangle(pixsize, pixsize)
        r += Point((j-1)*pixsize, (s[1]-i)*pixsize)
        render!(c, r, Rectangles.Plain(), GDSMeta(layers[i,j]))
    end
    c
end

"""
    jj!{T}(c::CoordinateSystem{T}, m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
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
function jj!(c::CoordinateSystem{T}, m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
        jj_meta::Meta, uc_meta::Meta) where {T}

    # start drawing JJ
    r1 = Rectangle(Point(zero(T), zero(T)), Point(w1,l1-w4))
    r4 = Rectangle(Point(-l4/2, l1-w4), Point(w1+l4/2, l1))
    r2 = Rectangle(Point(-m,l1+b), Point(w1+l2+w3,l1+b+w2))
    r3 = Rectangle(Point(w1+l2, l1+b+w2), Point(w1+l2+w3, l1+b+w2+l3))

    # join r2, r3 into one polygon
    p14 = clip(ClipTypeUnion, r1, r4)[1]
    p23 = clip(ClipTypeUnion, r2, r3)[1]

    # Generate undercut
    ucb2 = uc + t1*tan(ϕ)
    u1 = offset(Rectangle(Point(zero(T), zero(T)), Point(w1,l1+b)), uc)[1]
    u23 = offset(p23, uc)[1]
    uca = Rectangle(Point(-m-uc, l1+b+w2), Point(w1+l2, uc+l1+b+w2+t1*tan(Θ)))
    ucb1 = Rectangle(Point(w1, zero(T)), Point(w1+uc+t1*tan(ϕ), l1+b))
    ucb2 = Rectangle(Point(w1+l2+w3, l1+b-uc), Point(w1+l2+w3+uc+t1*tan(ϕ), l1+b+w2+l3))
    u = clip(ClipTypeUnion, u1, u23)[1]
    u = clip(ClipTypeUnion, u, uca)[1]
    u = clip(ClipTypeUnion, u, ucb1)[1]
    u = clip(ClipTypeUnion, u, ucb2)[1]

    # Remove undercut sticking out top and bottom
    u = clip(ClipTypeIntersection, u,
        Rectangle(Point(T(-1000000), zero(T)), Point(T(1000000),l1+b+w2+l3)))[1]
    u = clip(ClipTypeDifference, u, p14)
    uclip = clip(ClipTypeDifference, u, [p23]) # could be an array
    layers = fill(uc_meta, length(uclip))
    push!(uclip, p14, p23)
    push!(layers, jj_meta, jj_meta)

    # Horizontally center the shapes
    cen = center(bounds(uclip))
    p = Point(getx(cen),zero(getx(cen)))
    uclip .-= StaticArrays.Scalar{typeof(p)}((p,))

    # Render them
    for (x,y) in zip(uclip, layers)
        render!(c, x, Polygons.Plain(), y)
    end

    c
end

"""
    bandage_hatching!{T}(c::Cell{T}, pix, pixpattern, ruc, uuc, undercut_meta;
        jj_direction=0)
This function is helpful in generating a hatching pattern so that two layers of
aluminum can be deposited at different steps and electrically connected by ion milling
and subsequent evaporation of a bandage layer without breaking vacuum.

- `pix`: Size of a pixel in generating the hatching pattern. Passed to `layerpixels!`.
- `pixpattern`: A matrix to pass to `layerpixels`, e.g. [1 1; -1 1]
- `ruc`: Right undercut dimension
- `uuc`: Upper undercut dimension
- `undercut_meta`: Undercut `Meta` info.
- `jj_direction`: Orientation of the josephson junction with respect to qubit direction.
    (0: 0°, 1: 90°, 2: 180°, and 3: 270°)
"""
function bandage_hatching!(c::CoordinateSystem{T}, pix::T, pixpattern::T, ruc::T, uuc::T,
    undercut_meta::Meta; jj_direction=0) where {T <: Coordinate}
    layerpixels!(c, pixpattern, pix)
    plgs = clip(Clipper.ClipTypeUnion, polygon.(elements(c)), Polygon{T}[])
    c.elements = CellPolygon.(plgs, [meta(elements(c)[1])])

    r1 = Rectangle(pix, uuc) + Point(zero(T), 2 * pix - uuc)
    r2 = Rectangle(ruc, pix) + Point(2 * pix - ruc, zero(T))
    r3 = Rectangle(ruc, uuc) + Point(2 * pix - ruc, 2 * pix-uuc)

    plgs2 = clip(Clipper.ClipTypeDifference, polygon.(elements(c)), Polygon{T}[r1,r2,r3])
    origin = Point(zero(T), zero(T))
    if jj_direction == 1
        origin = Point(-2 * pix, zero(T))
    elseif jj_direction == 2
        origin = Point(-2 * pix, -2 * pix)
    elseif jj_direction == 3
        origin = Point(zero(T), -2 * pix)
    end
    trans = Translation(-origin) ∘ Rotation(π/2 * jj_direction)
    c.elements = CellPolygon.([trans(p) for p in plgs2], [meta(elements(c)[1])])

    for r in (r1, r2, r3)
        render!(c, trans(r), undercut_meta)
    end
    c
end

"""
    jjPadTopSide(pixcell, pix, uc, uuc, ruc, ext, padh, bandaid_buffer;
        jj_direction=0)
- `pixcell`: Cell containing a "pixel" of the bandaid hatching
- `pix`: size of the "pixel"
- `uc`: safety margin (additional undercut beyond what is needed from geometry)
- `uuc`: Upper undercut
- `ruc`: Right undercut
- `ext`: Width of rectangular landing pad
- `padh`: Height of rectangular landing pad
- `bandaid_buffer`: How far the bandaid should extend beyond the patch (ea. side)
- `jj_direction`: Orientation of the josephson junction with respect to qubit direction.
    (0: 0°, 1: 90°, 2: 180°, and 3: 270°)
"""
function jjPadTopSide(
    pixcell::CoordinateSystem{T}, pix::S, uc::S, uuc::S, ruc::S, ext::S, padh::S,
    bandaid_buffer::S; jj_direction=0) where {T, S <: Coordinate}
    jjpadtop = Cell(uniquename("jjpadtop"), nm)
    # initialize pix cells
    pixtop = Cell(uniquename("pixtop"), nm)
    pixleft = Cell(uniquename("pixleft"), nm)
    pixright = Cell(uniquename("pixright"), nm)
    pixbot = Cell(uniquename("pixbot"), nm)
    # if jj_direction = 1, (uuc, ruc) becomes (ruc, duc)
    # if jj_direction = 2, (uuc, ruc) becomes (duc, luc)
    # if jj_direction = 3, (uuc, ruc) becomes (luc, uuc)
    if jj_direction == 0
        render!(pixbot, Rectangle(Point(pix, -uuc), Point(2*pix - ruc, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixtop, Rectangle(Point(pix - uc, 0.0um), Point(2*pix, uuc)), GDSMeta(UC_LAYER))
        render!(pixright, Rectangle(Point(0.0um, pix - uc), Point(ruc, 2*pix - uuc + uc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 1
        render!(pixbot, Rectangle(Point(uuc, -uuc), Point(pix, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixbot, Rectangle(Point(0.0um, -uuc), Point(uuc, 0.0um)), GDSMeta(UC_LAYER))
        render!(pixleft, Rectangle(Point(-uuc, pix-uc), Point(0.0um, 2 * pix)), GDSMeta(UC_LAYER))
        render!(pixtop, Rectangle(Point(uuc - uc, 0.0um), Point(pix + uc, ruc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 2
        render!(pixbot, Rectangle(Point(ruc, -uuc), Point(pix, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixbot, Rectangle(Point(0.0um, -uuc), Point(ruc, 0.0um)), GDSMeta(UC_LAYER))
        render!(pixleft,Rectangle(Point(-ruc, uuc - uc), Point(0.0um, pix + uc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 3
        render!(pixbot, Rectangle(Point(pix, -uuc), Point(2 * pix - uuc, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixright, Rectangle(Point(0.0um, 0.0um), Point(uuc, pix + uc)), GDSMeta(UC_LAYER))
        render!(pixbot, Rectangle(Point(2 * pix - uuc, -uuc), Point(2*pix, 0.0um)), GDSMeta(UC_LAYER))
    end

    uc_array = circshift(Any[0.0um, uuc, ruc, 0.0um], -jj_direction)
    # Big contact pad
    render!(jjpadtop, Rectangle(ext, padh),
        Rectangles.Undercut(push!(copy(uc_array), GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER))...))

    # left, top, right, bottom
    luc, tuc, ruc, buc = uc_array

    # Hatching
    bnd = bounds(pixcell)

    # width and height of layer pixel
    lpw, lph = width(bnd), height(bnd)

    ya, yb, yc = promote(padh + uuc, lph, padh + lph + uuc)
    yrange = ustrip(ya):ustrip(yb):ustrip(yc)
    top = yrange[end] * unit(ya) + lph
    bottom = yrange[1] * unit(ya)

    xa, xb, xc = promote(0.0um, lpw, ext - lpw)
    xrange = ustrip(xa):ustrip(xb):ustrip(xc)
    right = xrange[end] * unit(xa) + lpw
    left = xrange[1] * unit(xa)

    push!(jjpadtop.refs,
        CellArray(pixcell, Point(xrange[1] * unit(xa), yrange[1] * unit(ya));
        nc = length(xrange), nr = length(yrange),
        dc = Point(step(xrange) * unit(xa), 0.0um),
        dr = Point(0.0um, step(yrange) * unit(ya))))
    # pixtop
    push!(jjpadtop.refs, CellArray(pixtop, Point(xrange[1] * unit(xa), top);
        nc = length(xrange), nr = 1,
        dc = Point(step(xrange) * unit(xa), 0.0um),
        dr = Point(0.0um, 0.0um)))
    # pixleft
    push!(jjpadtop.refs, CellArray(pixleft, Point(left, yrange[1] * unit(ya));
        nc = 1, nr = length(yrange),
        dc = Point(0.0um, 0.0um),
        dr = Point(0.0um, step(yrange) * unit(ya))))
    # pixright
    push!(jjpadtop.refs, CellArray(pixright, Point(right, yrange[1] * unit(ya));
        nc = 1, nr = length(yrange),
        dc = Point(0.0um, 0.0um),
        dr = Point(0.0um, step(yrange) * unit(ya))))

    # pixbot
    push!(jjpadtop.refs, CellArray(pixbot, Point(xrange[1] * unit(xa), bottom);
        nc = length(xrange), nr = 1,
        dc = Point(step(xrange) * unit(xa), 0.0um),
        dr = Point(0.0um, 0.0um)))

    # Bandaid
    patch = Rectangle(ext + ruc + 2 * bandaid_buffer,
        length(yrange) * lph + 2 * uuc + 2 * bandaid_buffer)
    render!(jjpadtop,
        patch + Point(-bandaid_buffer, padh - bandaid_buffer),
        Rectangles.Rounded(0.3um), GDSMeta(BRIDGE_FEET_LAYER))

    # Center in x (up to undercut and bandaid); bottom edge of elements will sit at y = 0
    jjpadtop -= Point(ext/2, 0um)
    # jjpadtop -= Point(ext/2, gety(center(jjpadtop)))
    # jjpadtop += Point(0um, height(bounds(jjpadtop))/2)
    jjpadtop
end

"""
    jjPadBotSide(pixcell, pix, uc, uuc, ruc, ext, padh, bandaid_buffer,
        armlength, armwidth, lr; jj_direction=0)
- `pixcell`: Cell containing a "pixel" of the bandaid hatching
- `pix`: size of the "pixel"
- `uc`: safety margin (additional undercut beyond what is needed from geometry)
- `uuc`: Upper undercut
- `ruc`: Right undercut
- `ext`: Width of rectangular landing pad
- `padh`: Height of rectangular landing pad
- `bandaid_buffer`: How far the bandaid should extend beyond the patch (ea. side)
- `armlength`: Length of the arm connecting the landing pad and the hatching pattern
- `armwidth`: Width of the arm connecting the landing pad and the hatching pattern
- 'lr': `-1` for left jj, `+1` for right jj
- `jj_direction`: Orientation of the josephson junction with respect to qubit direction.
    (0: 0°, 1: 90°, 2: 180°, and 3: 270°)
"""
function jjPadBotSide(pixcell, pix, uc, uuc, ruc, ext, padh, bandaid_buffer,
    armlength, armwidth, lr; jj_direction=0)

    jjpadbot = Cell(uniquename("jjpadbot"), nm)

    # initialize pix cells
    pixtop = Cell(uniquename("pixtop"), nm)
    pixleft = Cell(uniquename("pixleft"), nm)
    pixright = Cell(uniquename("pixright"), nm)
    pixbot = Cell(uniquename("pixbot"), nm)
    # if jj_direction = 1, (uuc, ruc) becomes (ruc, duc)
    # if jj_direction = 2, (uuc, ruc) becomes (duc, luc)
    # if jj_direction = 3, (uuc, ruc) becomes (luc, uuc)
    if jj_direction == 0
        # render!(pixbot, Rectangle(Point(pix, -uuc), Point(2*pix - ruc, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixtop, Rectangle(Point(pix - uc, 0.0um), Point(2*pix, uuc)), GDSMeta(UC_LAYER))
        render!(pixright, Rectangle(Point(0.0um, pix - uc), Point(ruc, 2*pix - uuc + uc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 1
        # render!(pixbot, Rectangle(Point(ruc, -uuc), Point(pix, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixleft, Rectangle(Point(-uuc, pix-uc), Point(0.0um, 2 * pix)), GDSMeta(UC_LAYER))
        render!(pixtop, Rectangle(Point(uuc - uc, 0.0um), Point(pix + uc, ruc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 2
        # render!(pixbot, Rectangle(Point(ruc, -uuc), Point(pix, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixbot, Rectangle(Point(0.0um, -uuc), Point(pix + uc, 0.0um)), GDSMeta(UC_LAYER))
        render!(pixleft,Rectangle(Point(-ruc, uuc - uc), Point(0.0um, pix + uc)), GDSMeta(UC_LAYER))
    elseif jj_direction == 3
        # render!(pixbot, Rectangle(Point(pix, -uuc), Point(2*pix - ruc, 0.0um)), GDSMeta(JJ_LAYER))
        render!(pixright, Rectangle(Point(0.0um, 0.0um), Point(uuc, pix + uc)), GDSMeta(UC_LAYER))
        render!(pixbot, Rectangle(Point(pix - uc, -ruc), Point(2*pix - uuc + uc, 0.0um)), GDSMeta(UC_LAYER))
    end

    b = bounds(pixcell)
    lpw, lph = width(b), height(b)

    # Big contact pad and "noodly appendages"
    uc_array = circshift(Any[0.0um, uuc, ruc, 0.0um], -jj_direction)

    render!(jjpadbot, Rectangle(ext, padh),
        Rectangles.Undercut(push!(copy(uc_array), GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER))...))
    render!(jjpadbot, Rectangle(armlength, armwidth) +
        ifelse(lr == -1, Point(-armlength, 0.0um), Point(ext, 0.0um)),
        Rectangles.Undercut(push!(copy(uc_array), GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER))...))

    # left, top, right, bottom
    luc, tuc, ruc, buc = uc_array

    # Hatching
    b = bounds(pixcell)
    lpw, lph = width(b), height(b)

    ya, yb, yc = promote(-3 * lph + uuc, lph, -lph + uuc + 20nm)
    yrange = ustrip(ya):ustrip(yb):ustrip(yc)
    top = yrange[end] * unit(ya) + lph
    bottom = yrange[1] * unit(ya)

    xa, xb, xc = promote(0.0um, lpw, ext - lpw)
    xrange = ustrip(xa):ustrip(xb):ustrip(xc)
    right = xrange[end] * unit(xa) + lpw
    left = xrange[1] * unit(xa)

    patch_adj_point = ifelse(lr == -1,
        Point(-armlength, pix-uuc), Point(ext + armlength - lpw * length(xrange), pix-uuc))

    # pixcell
    push!(jjpadbot.refs,
        CellArray(pixcell, Point(xrange[1] * unit(xa), yrange[1] * unit(ya)) +
            patch_adj_point;
            dc = Point(step(xrange)*unit(xa), 0.0um), dr = Point(0.0um, step(yrange)*unit(ya)),
            nc = length(xrange), nr = length(yrange)))

    # pixtop
    push!(jjpadbot.refs, CellArray(pixtop, Point(xrange[1] * unit(xa), top) +
        patch_adj_point;
        nc = length(xrange), nr = 1,
        dc = Point(step(xrange) * unit(xa), 0.0um),
        dr = Point(0.0um, 0.0um)))
    # pixleft
    push!(jjpadbot.refs, CellArray(pixleft, Point(left, yrange[1] * unit(ya)) +
        patch_adj_point;
        nc = 1, nr = length(yrange),
        dc = Point(0.0um, 0.0um),
        dr = Point(0.0um, step(yrange) * unit(ya))))
    # pixright
    push!(jjpadbot.refs, CellArray(pixright, Point(right, yrange[1]*unit(ya)) +
        patch_adj_point;
        nc = 1, nr = length(yrange),
        dc = Point(0.0um, 0.0um),
        dr = Point(0.0um, step(yrange) * unit(ya))))
    # pixbot
    push!(jjpadbot.refs, CellArray(pixbot, Point(xrange[1]*unit(xa), bottom) +
        patch_adj_point;
        nc = length(xrange), nr = 1,
        dc = Point(step(xrange) * unit(xa), 0.0um),
        dr = Point(0.0um, 0.0um)))

    # Bandaid
    patch = Rectangle(length(xrange) * lpw + 2 * bandaid_buffer,
        length(yrange) * lph + 2 * bandaid_buffer - uuc)
    Rec_round = 0.4um
    render!(jjpadbot,
        patch + Point(-bandaid_buffer, -height(patch) - Rec_round * 2) + patch_adj_point,
        Rectangles.Rounded(Rec_round), GDSMeta(BRIDGE_FEET_LAYER))

    # top edge of elements will sit at y = uuc
    jjpadbot -= Point(ext/2, padh)

    # Flatten and take union of full dose elements
    flatten!(jjpadbot)
    inds = findall(x->layer(x) == JJ_LAYER, jjpadbot.elements)
    res = clip(Clipper.ClipTypeUnion, polygon.(view(jjpadbot.elements, inds)), Polygon{typeof(1.0nm)}[],
        pfs=Clipper.PolyFillTypePositive)
    deleteat!(jjpadbot.elements, inds)
    append!(jjpadbot.elements, CellPolygon.(res, [GDSMeta(JJ_LAYER)]))

    jjpadbot
end


"""
jjComplt(m, b, w1, w2, w3, w4, l1, l2, l3, l4, uc, t1, Θ, ϕ,
    qubit_bottom_gap, jj_height, gap_between_leads_for_jjs, pix,
    uuc, ruc, junc_pad_spacing, lead_width, z_qb_dist, Leg_totH, n; jj_direction = 0)
The function to put all peices for juctions together.
    n is to distinguish between different cells.
Return the jj_cell
"""
function jjComplt(m, b, w1L, w1R, w2L, w2R, w3, w4, l1, l2, l3L, l3R, l4, uc, t1, Θ, ϕ,
    qubit_bottom_gap, jj_height, gap_between_leads_for_jjs, pix,
    uuc, ruc, junc_pad_spacing, lead_width, z_qb_dist, Leg_totH, n; jj_direction = 0)

    c = Cell("hatching", nm)
    pixcell = bandage_hatching!(c, pix, pixpattern, ruc, uuc,
        GDSMeta(UC_LAYER); jj_direction=jj_direction)
    jjCellL = Cell("jjCellL", nm)
    jjCellR = Cell("jjCellR", nm)
    jj!(jjCellL, m, b, w1L, w2L, w3, w4, l1, l2, l3L, l4, uc, t1, Θ, ϕ,
        GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER))
    jj!(jjCellR, m, b, w1R, w2R, w3, w4, l1, l2, l3R, l4, uc, t1, Θ, ϕ,
        GDSMeta(JJ_LAYER), GDSMeta(UC_LAYER))
    jj_bound = bounds(jjCellL)
    jw, jh = width(jj_bound), height(jj_bound)
    # println("func jw: $jw, jh: $jh ")
    jj_rot_offset = Point(0um, 0um)
    # jj origin offset under rotation
    if jj_direction == 2
        jj_rot_offset += Point(0um, jh)
    elseif jj_direction == 1
        jj_rot_offset += Point(jh/2,jw/2)
    elseif jj_direction == 3
        jj_rot_offset += Point(-jh/2, jw/2)
    end
    left_jj_with_pad = Cell("left_junction_with_pad$n", nm)
    right_jj_with_pad = Cell("right_junction_with_pad$n", nm)
    # jj with pads
    for (j, lr) in zip((left_jj_with_pad, right_jj_with_pad), (-1, 1))
        # jj
        if lr == -1
            push!(j.refs, CellReference(jjCellL, jj_rot_offset;
                rot=jj_direction * π/2))
        elseif lr == +1
            push!(j.refs, CellReference(jjCellR, jj_rot_offset;
                rot=jj_direction * π/2))
        end

        ext = 2 * lead_width + 1um
        padh = ifelse(jj_direction % 2 ==0, (qubit_bottom_gap - jj_height)/2 + 0.5um,
            (qubit_bottom_gap + gap_between_leads_for_jjs) / 2 + 0.5um)
        bandaid_buffer = 0.5um #ifelse(jj_direction % 2 ==0, 0.5um, 0.7um)
        armlength = ifelse(jj_direction % 2 ==0,
            9um, 9um - (lead_width + gap_between_leads_for_jjs))
        armwidth = pix

        jjpadtop = jjPadTopSide(pixcell, pix, uc, uuc, ruc, ext, padh + 2um,
            bandaid_buffer; jj_direction=jj_direction)
        jjpadbot = jjPadBotSide(pixcell, pix, uc, uuc, ruc, ext, padh,
            bandaid_buffer, armlength, armwidth, lr; jj_direction=jj_direction)

        # Junction pads
        toppad_loc = ifelse(jj_direction % 2 == 0,
            Point(-lr * lead_width/2, jh),
            Point(-lr * (jh + ext) / 2, (jw/2 - lead_width/2 - 0.5um)) )
        botpad_loc = ifelse(jj_direction % 2 == 0,
            Point(lr * lead_width/2, 0.0um),
            Point(lr * (ext + jh) / 2, (jw/2 + lead_width/2 + 0.5um)))

        push!(j.refs, CellReference(jjpadtop, toppad_loc))
        push!(j.refs, CellReference(jjpadbot, botpad_loc))
        flatten!(j)

        # union of undercuts
        inds = findall(x->layer(x) == UC_LAYER, j.elements)
        res_ = clip(Clipper.ClipTypeUnion, polygon.(view(j.elements, inds)), Polygon{typeof(1.0um)}[],
            pfs=Clipper.PolyFillTypePositive)
        deleteat!(j.elements, inds)
        append!(j.elements, CellPolygon.(res_, [GDSMeta(UC_LAYER)]))

        # no overlaps of UC_LAYER and LAYER
        inds = findall(x->layer(x) == UC_LAYER, j.elements)
        inds2 = findall(x->layer(x) == JJ_LAYER, j.elements)
        res_ = clip(Clipper.ClipTypeDifference, polygon.(view(j.elements, inds)),
            polygon.(view(j.elements, inds2)))
        deleteat!(j.elements, inds)
        append!(j.elements, CellPolygon.(res_, [GDSMeta(UC_LAYER)]))
    end
    jj_hori_offset = ifelse(jj_direction % 2 == 0,
        (junc_pad_spacing + 1.5 * lead_width),
        (junc_pad_spacing + 2 * lead_width + 0.5um + jj_height / 2))
    jj_vert_offset = ifelse(jj_direction % 2 == 0,
        -(qubit_bottom_gap - gap_between_leads_for_jjs)/2 -
            (gap_between_leads_for_jjs - jj_height)/2 - jj_height,
            -(qubit_bottom_gap + jw)/2)
    # println("function")
    jj_cell = Cell("jj_with_pads$n", nm)
    push!(jj_cell.refs,
        CellReference(left_jj_with_pad,
            Point(-jj_hori_offset, jj_vert_offset)))
    push!(jj_cell.refs,
        CellReference(right_jj_with_pad,
            Point(jj_hori_offset, jj_vert_offset)))
    jj_cell
end