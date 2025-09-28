module CPWResonator

using DeviceLayout
import DeviceLayout: Meta, GDSMeta, Coordinate
import Unitful: unit

"""
    gap!(p::Path{T}, cpwr_trace, cpwr_gap) where {T<:Coordinate}

Insert a **capacitive gap** at the current end of the path `p` for a CPW line.

This draws a short straight segment whose *trace* is widened to `cpwr_trace + 2*cpwr_gap`
and whose *length* is `cpwr_gap`, producing an open-ended capacitive break that matches the
surrounding CPW geometry.

# Arguments
- `p`: A `DeviceLayout.Paths.Path` being constructed.
- `cpwr_trace`: Center trace width of the CPW (same units as the path).
- `cpwr_gap`: Gap to the ground planes; also used as the gap length here.

# Returns
The mutated path `p`.

# Notes
This is typically used at the **open** end(s) of a resonator or coupler to realize a capacitive break.
"""
function gap!(p::Path{T}, cpwr_trace, cpwr_gap) where {T<:Coordinate}
    straight!(p, cpwr_gap, Paths.Trace(cpwr_trace + 2*cpwr_gap))
    return p
end

cpwr_error_msg = "CPWResonator Error : "

"""
    short_error()

Throw a standardized error indicating that the **specified resonator length
is too short** to complete the next geometric feature (e.g., a bend, taper, or straight).

This is used internally by the segment builders when `remaining_length` is insufficient.
"""
function short_error()
    error_msg =
        string(cpwr_error_msg,
            "specified length of the resonator too short.")
    error(error_msg)
end


"""
    cpwr_segment_vertical!(
        p::Path{T}, remaining_length, res_coupler_len,
        res_coupler_len2, res_radius, res_straight_length,
        coupler_trace, coupler_gap, cpwr_trace, cpwr_gap,
        open_start::Bool, open_end::Bool
    ) where {T<:Coordinate}

Advance the path `p` by **one vertical meander segment** of a coplanar-waveguide (CPW)
half-wave resonator, including its input **coupler** and subsequent **meanders**.

This routine mutates `p` by appending the next geometric feature based on the current
segment count (`length(p)`), honoring open-ended capacitive gaps at the start/end, and
consuming as much of `remaining_length` as required for the next feature. It is designed to
be called repeatedly until the target total length is reached.

# Arguments
- `p`: The path under construction (mutated in place).
- `remaining_length`: Length left to draw (including any end-gap allowance if applicable).
- `res_coupler_len`: Horizontal coupler straight length near the feed line.
- `res_coupler_len2`: Secondary straight/taper run used to **taper** from coupler CPW to resonator CPW.
- `res_radius`: Bend radius used for 90°/180° turns.
- `res_straight_length`: Nominal long straight used in the meander sections.
- `coupler_trace`, `coupler_gap`: CPW dimensions for the **coupler** section.
- `cpwr_trace`, `cpwr_gap`: CPW dimensions for the **resonator** proper.
- `open_start`, `open_end`: If `true`, create a capacitive `gap!` at the start/end.

# Behavior
- First few calls build the **coupler** (open gap, horizontal straight, 90° bend, taper).
- Subsequent calls build a vertical **meander** using straights and ±π turns.
- If `remaining_length` is insufficient for the next required feature, the function either:
  - draws a shortened feature and, if `open_end` is `true`, terminates with a `gap!`, or
  - throws `short_error()` when an open-end termination is not allowed in that context.

# Returns
The mutated path `p`.
"""
function cpwr_segment_vertical!(p::Path{T}, remaining_length, res_coupler_len,
    res_coupler_len2, res_radius, res_straight_length, coupler_trace, coupler_gap,
    cpwr_trace, cpwr_gap, open_start::Bool, open_end::Bool) where {T<:Coordinate}
    # coplanar waveguide resonator style
    coupler_style = Paths.CPW(coupler_trace, coupler_gap)
    taper_style = Paths.TaperCPW(coupler_trace, coupler_gap, cpwr_trace, cpwr_gap)
    cpwr_style = Paths.CPW(cpwr_trace, cpwr_gap)
    # Coupler section
    (length(p) <= 4) ? begin
        if length(p) == 0
            # gap at the start if `open_start` is `true`
            open_start ? gap!(p, coupler_trace, coupler_gap) : nothing
            # coupler (horizontal)
            (remaining_length > (res_coupler_len + open_end * coupler_gap)) ?
            begin
                straight!(p, res_coupler_len, coupler_style)
                simplify!(p)
            end : short_error()
        elseif length(p) == 1
            # 90° bend
            (remaining_length > (π * res_radius / 2 + open_end * coupler_gap)) ?
            turn!(p, -π/2, res_radius, coupler_style) : short_error()
        elseif length(p) == 2
            # keep distance from coupler
            (remaining_length > (res_coupler_len2 + open_end * cpwr_gap)) ?
            begin
                straight!(p, res_coupler_len2 / 3, coupler_style)

                taper_len = res_coupler_len2 / 3
                straight!(p, taper_len, taper_style)
                straight!(p, res_coupler_len2 / 3, cpwr_style)
                simplify!(p, (length(p)-2):length(p))
            end : short_error()
        elseif length(p) == 3
            (remaining_length > (π * res_radius / 2 + open_end * cpwr_gap)) ?
            turn!(p, -π/2, res_radius, cpwr_style) :
            begin 
                turn!(p, -(remaining_length - open_end * cpwr_gap) / res_radius)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : short_error()
            end
        elseif length(p) == 4
            # cpw
            (remaining_length > ((res_straight_length - 2 * res_radius + res_coupler_len) / 2 + open_end * cpwr_gap)) ?
            straight!(p, (res_straight_length - 2 * res_radius + res_coupler_len) / 2, cpwr_style) :
            begin
                straight!(p, remaining_length - open_end * cpwr_gap)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : short_error()
            end
        end
    # end of coupler
    end : begin
    # meander for resonator starts
        if length(p) % 2 == 0
            (remaining_length > (res_straight_length + open_end * cpwr_gap)) ?
            straight!(p, res_straight_length) :
            begin
                straight!(p, remaining_length - open_end * cpwr_gap)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        elseif length(p) % 4 == 1
            (remaining_length > (π * res_radius + open_end * cpwr_gap)) ?
            turn!(p, π, res_radius) :
            begin
                turn!(p, (remaining_length - open_end * cpwr_gap) /
                    res_radius, res_radius)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        elseif length(p) % 4 == 3
            (remaining_length > (π * res_radius + cpwr_gap)) ?
            turn!(p, -π, res_radius) :
            begin
                turn!(p, -(remaining_length - open_end * cpwr_gap) /
                    res_radius, res_radius)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        end
    end
    return p
end

"""
    cpwr_segment_horizontal!(
        p::Path{T}, remaining_length, res_coupler_len, res_coupler_len2,
        res_radius, res_straight_length, coupler_trace, coupler_gap,
        cpwr_trace, cpwr_gap, open_start::Bool, open_end::Bool
    ) where {T<:Coordinate}

Advance the path `p` by **one horizontal meander segment** of a CPW half-wave
resonator, including coupler construction and subsequent meanders laid out primarily
along the horizontal axis.

# Arguments
Same as [`cpwr_segment_vertical!`](@ref), but the segment sequencing produces a
horizontal progression (useful when `is_vertical = false` in the top-level API).

# Behavior
- Builds: start gap (optional) → horizontal coupler → 90° bend → taper to resonator CPW →
  initial straight → meander using straights and ±π turns.
- Includes a small epsilon guard: if `remaining_length < cpwr_gap`, a short straight with
  a neutral `Trace` is drawn to consume the residue (useful for floating-point leftovers).
- On insufficient `remaining_length`, it either shortens and optionally terminates with an
  end gap, or throws `short_error()` depending on context.

# Returns
The mutated path `p`.
"""
function cpwr_segment_horizontal!(p::Path{T}, remaining_length,
    res_coupler_len, res_coupler_len2, res_radius, res_straight_length,
    coupler_trace, coupler_gap, cpwr_trace, cpwr_gap, open_start::Bool,
    open_end::Bool) where {T<:Coordinate}
    # coplanar waveguide resonator style
    coupler_style = Paths.CPW(coupler_trace, coupler_gap)
    taper_style = Paths.TaperCPW(coupler_trace, coupler_gap, cpwr_trace, cpwr_gap)
    cpwr_style = Paths.CPW(cpwr_trace, cpwr_gap)
    # Coupler section
    (length(p) <= 3) ? begin
        if length(p) == 0
            # capacitive gap
            open_start ? gap!(p, coupler_trace, coupler_gap) : nothing
            # coupler (horizontal)
            (remaining_length > (res_coupler_len + open_end * coupler_gap)) ?
            begin
                straight!(p, res_coupler_len, coupler_style)
                simplify!(p)
            end : short_error()
        elseif length(p) == 1
            # 90° bend
            (remaining_length > (π * res_radius / 2 + open_end * cpwr_gap)) ?
            turn!(p, -π/2, res_radius, coupler_style) : short_error()
        elseif length(p) == 2
            # keep distance from coupler
            (remaining_length > (res_coupler_len2 + open_end * cpwr_gap)) ?
            begin
                straight!(p, res_coupler_len2 / 3, coupler_style)
                taper_len = res_coupler_len2 / 3
                straight!(p, taper_len, taper_style)
                straight!(p, res_coupler_len2 / 3, cpwr_style)
                simplify!(p, (length(p)-2):length(p))
            end : short_error()
        elseif length(p) == 3
            # cpw
            (remaining_length > (res_radius + open_end * cpwr_gap)) ?
            straight!(p, res_radius, cpwr_style) : short_error()
        end
    end : begin
        # to compensate for any epsilon float amount left in the remaining_length
        # will fix later
        (remaining_length < cpwr_gap) ?
            straight!(p, remaining_length,
                Paths.Trace(open_end * (cpwr_trace + 2 * cpwr_gap))) :
        if length(p) % 2 == 0
            (remaining_length > (res_straight_length + open_end * cpwr_gap)) ?
            straight!(p, res_straight_length) :
            begin
                straight!(p, remaining_length - open_end * cpwr_gap)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        elseif length(p) % 4 == 1
            (remaining_length > (π * res_radius + open_end * cpwr_gap)) ?
            turn!(p, -π, res_radius) :
            begin
                turn!(p, -(remaining_length - open_end * cpwr_gap) /
                    res_radius, res_radius)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        elseif length(p) % 4 == 3
            (remaining_length > (π * res_radius + open_end * cpwr_gap)) ?
            turn!(p, π, res_radius) :
            begin
                turn!(p, (remaining_length - open_end * cpwr_gap) /res_radius, res_radius)
                open_end ? gap!(p, cpwr_trace, cpwr_gap) : nothing
            end
        end
    end
    return p
end

"""
    cpw_resonator!(
        c::Cell{T}, length, res_coupler_len, res_coupler_len2, res_radius,
        res_straight_length, cpwr_trace, cpwr_gap,
        meta::DeviceLayout.Meta = GDSMeta(0, 0);
        open_start::Bool = true, open_end::Bool = true,
        is_vertical::Bool = true,
        coupler_trace = cpwr_trace, coupler_gap = cpwr_gap,
        path_output::Bool = false
    ) where {T<:Coordinate}

Draw a **CPW resonator** with a simple capacitive input coupler and a
meandered body, and render it into the given `Cell`. Optionally return the `Path`.

This function repeatedly calls either [`cpwr_segment_vertical!`](@ref) or
[`cpwr_segment_horizontal!`](@ref), depending on `is_vertical`, until the total drawn
path length (minus optional end gaps) equals the requested `length` within the
segment discretization.

# Required Arguments
- `c`: The `DeviceLayout.Cell` to render into (mutated in place).
- `length`: Target resonator **electrical length** to draw (unitful or native to your
  coordinate type). The algorithm will attempt to match this length; small residuals may
  occur due to bends/tapers.
- `res_coupler_len`, `res_coupler_len2`: Straight lengths (see segment docstrings) defining
  the coupler and taper region.
- `res_radius`: Bend radius for all turns.
- `res_straight_length`: Straight length used in meanders.
- `cpwr_trace`, `cpwr_gap`: CPW dimensions of the **resonator** proper.

# Keywords
- `meta`: A `DeviceLayout.Meta` (e.g., `GDSMeta(layer, datatype)`) used by `render!`.
- `open_start`, `open_end`: If `true`, place a capacitive `gap!` at the start/end.
- `is_vertical`: If `true`, build vertical meanders; if `false`, horizontal meanders.
- `coupler_trace`, `coupler_gap`: CPW dimensions for the **coupler** (defaults to resonator CPW).
- `path_output`: If `true`, return the constructed `Path`; otherwise return the mutated `Cell`.

# Returns
- If `path_output == false` (default): the mutated `Cell c`.
- If `path_output == true`: the constructed `Path`.

# Warnings
If the final drawn length differs from the requested `length` (after subtracting any
`open_start`/`open_end` gaps), a warning is emitted:
`"CPWResonator Warning: Resonator length mismatch with the drawn path"`.

# Examples
```julia
using DeviceLayout, Unitful
using .CPWResonator

c = Cell{GDSMeta}("cpw_res")
cpw_resonator!(c,
    500u"µm",          # length
    40u"µm",           # res_coupler_len
    60u"µm",           # res_coupler_len2 (with taper)
    20u"µm",           # res_radius
    120u"µm",          # res_straight_length (meander straight)
    10u"µm", 6u"µm";   # cpwr_trace, cpwr_gap
    open_start=true, open_end=true, is_vertical=true,
    coupler_trace=8u"µm", coupler_gap=6u"µm",
    meta=GDSMeta(1, 0)
)

# Or, to inspect the path without rendering:
p = cpw_resonator!(c,
    500u"µm", 40u"µm", 60u"µm", 20u"µm", 120u"µm", 10u"µm", 6u"µm";
    path_output=true)
"""
function cpw_resonator!(c::Cell{T}, length, res_coupler_len,
    res_coupler_len2, res_radius, res_straight_length, cpwr_trace, cpwr_gap,
    meta::DeviceLayout.Meta=GDSMeta(0, 0); open_start::Bool=true, open_end::Bool=true,
    is_vertical::Bool=true, coupler_trace=cpwr_trace,
    coupler_gap=cpwr_gap, path_output::Bool=false) where {T<:Coordinate}
    p = Path(unit(zero(T)))
    remaining_length = length + open_start * coupler_gap + open_end * cpwr_gap
    if is_vertical == true
        while remaining_length > zero(T)
            cpwr_segment_vertical!(p, remaining_length, res_coupler_len,
                res_coupler_len2, res_radius, res_straight_length, coupler_trace,
                coupler_gap, cpwr_trace, cpwr_gap, open_start, open_end)
            remaining_length = (length + open_start * coupler_gap +
                open_end * cpwr_gap - pathlength(p))
        end
    elseif is_vertical == false
        while remaining_length > zero(T)
            cpwr_segment_horizontal!(p, remaining_length, res_coupler_len,
                res_coupler_len2, res_radius, res_straight_length, coupler_trace,
                coupler_gap, cpwr_trace, cpwr_gap, open_start, open_end)
            remaining_length = (length + open_start * coupler_gap +
                open_end * cpwr_gap - pathlength(p))
        end
    end

    render!(c, p, meta)
    if length != (pathlength(p) - open_start * coupler_gap - open_end * cpwr_gap)
        warn("CPWResonator Warning: Resonator length mismatch with the drawn path")
    end

    if path_output
        return p
    else
        return c
    end
end

# end of module
end