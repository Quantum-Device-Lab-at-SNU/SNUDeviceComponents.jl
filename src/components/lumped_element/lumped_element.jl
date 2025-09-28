module LumpedElement

using DeviceLayout
import DeviceLayout: Coordinate, CoordinateSystem, Meta

"""
    meander_wire_path!(p::Path{T}, t, s, n::Int, leg, r, lead; first_turn=:right) where {T<:Coordinate}

Constructs a meander (serpentine) wire pattern along an existing `Path` object `p`.

# Arguments
- `p::Path{T}`: The path object (from DeviceLayout) that will be modified in place.
- `t`: Trace width.
- `s`: Spacing (gap) between adjacent legs of the meander.
- `n::Int`: Number of meander turns (zig-zag repetitions).
- `leg`: Length of each leg segment between turns.
- `r`: Bend radius. 
    - If `r == 0`, sharp 90° turns are created using `corner!`.
    - If `r > 0`, smooth rounded turns are created using `turn!` with radius `r`.
- `lead`: Additional straight length before and after the meander.
- `first_turn`: Direction of the first turn, either `:right` (default) or `:left`.

# Behavior
- The function draws a meander starting with a straight lead segment.
- For each of the `n` turns:
  - The trace alternates direction (left/right) depending on `first_turn`.
  - Each leg is connected by either sharp or rounded turns, depending on `r`.
  - End legs (first and last) are shortened by half the trace width (and bend radius, if applicable) to ensure symmetric ports.
- Spacing `s` is inserted between successive legs.

# Returns
- Returns the modified path `p` with the meander geometry appended.

# Notes (Sanity Checks)
- `t > 0`, `s > 0`, `leg > 0`, and `n ≥ 2`.
- For sharp bends (`r == 0`): ensure `s ≥ t` to avoid trace overlap.
- For rounded bends (`r > 0`):
  - Require `r ≥ t/2` for manufacturable curvature.
  - Require `s ≥ t` and `s ≥ 2r` to fit a U-turn without self-intersection.
  - Require `leg ≥ 2r + t` to leave positive straight length between turns.
"""
function meander_wire_path!(
    p::Path{T},
    t, s, n::Int, leg, r, lead;
    first_turn = :right
) where {T<: Coordinate}

    # -------------
    # Sanity checks
    # -------------
    @assert t > zero(T) "Trace width (t) must be > 0"
    @assert s > zero(T) "Spacing (s) must be > 0"
    @assert leg > zero(T) "Leg length (leg) must be > 0"
    @assert n ≥ 2 "Number of turns (n) must be ≥ 2"

    if r == zero(T)
        @assert s ≥ t "For sharp bends (r=0), spacing (s) must be ≥ trace width (t) to avoid overlap"
        @assert leg ≥ t "For sharp bends, leg length must be ≥ trace width"
    elseif r > zero(T)
        @assert r ≥ t/2 "Bend radius must be ≥ t/2 for manufacturability"
        @assert s ≥ t "Spacing must be ≥ t to allow a U-turn without overlap"
        @assert s ≥ 2r "Spacing must be ≥ 2r to allow a U-turn"
        @assert leg ≥ 2r + t "Leg length must be ≥ 2r + t to leave straight length between bends"
    end

    # -------------------
    # Geometry generation
    # -------------------

    sty = Paths.Trace(t)
    # when r is zero, we use corner! to represent sharp turns
    if r == zero(T)
        straight!(p, lead + t / 2, sty)

        for i = 1:n
            sgn = (first_turn == :right) ? (-1) ^ i : (-1) ^ (i + 1)

            corner!(p, sgn * π/2, Paths.SimpleTraceCorner())

            if i in [1, n]
                straight!(p, leg / 2 - t / 2, sty)
            else
                straight!(p, leg - t)
            end
            corner!(p, - sgn * π/2, Paths.SimpleTraceCorner())

            if i < n
                straight!(p, s)
            end
        end
        straight!(p, lead + t / 2, sty)
    
    # when the specified `r` is not zero, use turn! to create rounded corners
    elseif r > zero(T)
        straight!(p, lead + t / 2 - r, sty)
        for i = 1:n
            sgn = (first_turn == :right) ? (-1) ^ i : (-1) ^ (i + 1)

            turn!(p, sgn * π/2, r)

            if i in [1, n]
                straight!(p, leg / 2 - 2 * r - t / 2, sty)
            else
                straight!(p, leg - 2 * r - t)
            end
            turn!(p, - sgn * π/2, r)

            if i < n
                straight!(p, s - 2 * r)
            end
        end
        straight!(p, lead + t / 2 - r, sty)
    end

    return p
end


"""
    meander_inductor!(c::CoordinateSystem{T}, clearance, t, s, n::Int, leg, r, lead, meta; first_turn=:right) where {T<:Coordinate}

Generate a **negatively-defined meander inductor** within the given `Cell` by subtracting
the drawn wire geometry from a rectangular ground opening.

# Arguments
- `c::CoordinateSystem{T}`: Target cell where the final meander inductor (negative definition) will be added.
- `clearance`: Vertical clearance (opening height) of the ground rectangle that accommodates the meander.
- `t`: Trace width of the meander wire.
- `s`: Spacing (gap) between adjacent meander legs.
- `n::Int`: Number of meander turns (zig-zag repetitions).
- `leg`: Length of each leg segment between turns.
- `r`: Bend radius.  
  - If `r == 0`, corners are sharp (90°).  
  - If `r > 0`, corners are rounded with radius `r`.
- `lead`: Extra straight section length at both ends of the meander.
- `meta`: Rendering metadata (layer, datatype, etc.).
- `first_turn`: Direction of the first turn, `:right` (default) or `:left`.

# Behavior
1. Builds a temporary child cell (`c0`) containing the **positively-defined** meander wire
   using [`meander_wire_path!`](@ref).
2. Unions the wire polygons into a single positive geometry (`wire_poly`).
3. Defines a ground rectangle of width `d` (the total horizontal span of the meander)
   and height `clearance`, centered vertically.
4. Subtracts (`difference2d`) the wire geometry from this rectangle to produce the
   **negative ground opening** that accommodates the meander.
5. Renders the result into the target cell `c`.

# Returns
- The modified `Cell` `c`, now containing the negatively-defined meander opening.

# Notes
- Use this function when your PDK expects the inductor to be represented as a
  **ground opening** (negative definition), rather than a positive metal trace.
- Ensure that `clearance` is chosen large enough so that the meander does not clip
  outside of the opening.
"""
function meander_inductor!(
    c::CoordinateSystem{T},
    clearance, t, s, n::Int, leg, r, lead, meta::Meta;
    first_turn=:right, clearance_corner_r=zero(T)) where {T<:Coordinate}

    c0 = Cell(uniquename("meander_wire"), unit(zero(T)))

    p = Path(unit(zero(T)))
    meander_wire_path!(p, t, s, n, leg, r, lead; first_turn=first_turn)
    render!(c0, p, meta)
    # take the union of the polygons forming the meander wire (positive)
    wire_poly = union2d(c0.elements)
    d = width(bounds(c0)) # horizontal length of the full meander section
    # rectangle to subtract from
    r = Rounded(Rectangle(d, clearance), clearance_corner_r) - Point(zero(T), clearance / 2)
    render!(c, difference2d(r, wire_poly), meta)
    return c
end

end # end of module