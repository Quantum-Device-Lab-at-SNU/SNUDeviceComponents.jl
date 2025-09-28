# Functions useful for ExamplePDK but not suitable/ready for base DeviceLayout/SchematicDrivenLayout

module Utils
using DeviceLayout
import DeviceLayout: Coordinate
"""
    generate_undercut(ps::AbstractArray{<:Polygon{T}}, luc, ruc, uuc, duc; α=0°) where {T<:Coordinate}

Compute a **directional undercut region** (an anisotropic offset “ring”) around one
or more polygons.

This constructs an oriented sweep/offset of the union of `ps` by expanding
(left/right/up/down) by the specified amounts in a local frame rotated by angle `α`,
then subtracts the original geometry. The result is the *undercut shell only*
(i.e., the grown area minus the original).

# Arguments
- `ps`: Collection of polygons to be undercut (all are unioned before processing).
- `luc`: Left undercut (≥ 0). Horizontal expansion toward negative local-x.
- `ruc`: Right undercut (≥ 0). Horizontal expansion toward positive local-x.
- `uuc`: Up undercut (≥ 0). Vertical expansion toward positive local-y.
- `duc`: Down undercut (≥ 0). Vertical expansion toward negative local-y.
- `α`: Rotation angle (in degrees or `°` units), defining the local sweep frame.
       `α = 0°` means left/right are along global −x/+x, up/down along +y/−y.

# Behavior
1. `ps` are unioned → `p0`.
2. Rotate by `-α` to enter the local frame.
3. Apply directional sweeps:
   - sweep left by `luc`, then right by `ruc`,
   - then up by `uuc`, then down by `duc`.
4. Rotate back by `+α` to global frame.
5. Subtract original (`p0`) from the expanded result → **undercut shell**.

Equivalent interpretation: an **anisotropic Minkowski expansion** by a
rectangle of half-widths `(ruc, luc, uuc, duc)` in the rotated local frame,
followed by removing the original interior.

# Returns
- A polygon set representing only the **undercut region** (expanded minus original).

# Notes
- Typical use: model etch bias / metal pull-back, solder-mask or ground-cut openings
  with different margins on each side; or generate clearance rings.
- All undercut distances should be non-negative. If a value is zero, no expansion
  occurs on that side.
- If you need the **fully grown shape** instead of the shell, use the intermediate
  expanded geometry (before subtraction) or compute `union2d(ps) ⊕ rect` directly.
"""
function generate_undercut(
    ps::AbstractArray{<:Polygon{T}}, luc, ruc, uuc, duc; α=0°
) where {T<:Coordinate}
    # rectangle that represents a set of displacement vectors to be applied to the original polygon

    p0 = union2d(ps)
    p1 = transform(p0, Rotation(-α))
    p2 = sweep_poly(p1, Point(-luc, zero(T)))
    p3 = sweep_poly(p2, Point(ruc, zero(T)))
    p4 = sweep_poly(p3, Point(zero(T), uuc))
    p5 = sweep_poly(p4, Point(zero(T), -duc))
    p6 = transform(p5, Rotation(α))

    undercut = difference2d(p6, p0)
    return undercut
end

"""
    generate_undercut(p::Polygon{T}, luc, ruc, uuc, duc; α=0°) where {T<:Coordinate}

Convenience wrapper for a single polygon.

Equivalent to:
`generate_undercut(Polygon{T}[p], luc, ruc, uuc, duc; α=α)`

See [`generate_undercut(::AbstractArray{<:Polygon})`](@ref) for details.
"""
generate_undercut(p::Polygon{T}, luc, ruc, uuc, duc; α=0°) where {T<:Coordinate} = generate_undercut(
    Polygon{T}[p], luc, ruc, uuc, duc; α=α
)

"""
    generate_undercut(p::Rectangle{T}, luc, ruc, uuc, duc; α=0°) where {T<:Coordinate}

Convenience wrapper for a rectangle.

The rectangle is first converted to polygons via `to_polygons(p)` and then
passed to the array method.

See [`generate_undercut(::AbstractArray{<:Polygon})`](@ref) for details.
"""
generate_undercut(p::Rectangle{T}, luc, ruc, uuc, duc; α=0°) where {T<:Coordinate} = generate_undercut(
    to_polygons(p), luc, ruc, uuc, duc; α=α
)


end # end of module