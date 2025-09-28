module AlignmentMarker

using DeviceLayout
import DeviceLayout: Meta, Coordinate, AbstractCoordinateSystem

# const KANC_MARKER_PITCH = 
# const KANC_MARKER_SIZE = 

"""
    cross_marker!(c::CoordinateSystem{T}, width::S, size::S, meta::Meta) where {T, S <: Coordinate}

Renders a symmetric 2D cross-shaped marker to the layout cell `c`, commonly used for lithographic alignment.

# Parameters
- `c`: Target coordinate system or layout cell.
- `width`: Width of each arm (limb) of the cross.
- `size`: Total size of the cross (i.e., extent of the bounding square in both X and Y directions).
- `meta`: Metadata specifying layer and datatype information.

# Description
The cross is centered at the origin `(0, 0)` and consists of a vertical and horizontal bar intersecting at the center. Each arm spans the full `size` along one axis and has a constant thickness of `width`.

# ASCII Diagram

    _ |<--------- size ----------->|
    |            |▓▓▓▓|
    |            |▓▓▓▓|
    size         |▓▓▓▓|
    |            |<  >| width
    |            |▓▓▓▓|
    | |▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓| width
    | |▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓|
    |            |▓▓▓▓|╲     
    |            |▓▓▓▓| ╲--- (x, y) = (0, 0)
    |            |▓▓▓▓|
    |            |▓▓▓▓|
    _            |▓▓▓▓|

"""
function cross_marker!(c::AbstractCoordinateSystem{T}, width::S, size::S, meta::Meta) where {T, S <: Coordinate}
    render!(c, simple_cross(width, size), meta)
end

"""
    ell_marker!(c::CoordinateSystem{T}, width::S, size::S, meta::Meta) where {T, S <: Coordinate}

Renders an L-shaped (ell) marker to the layout cell `c`, typically used as an alignment marker in lithography.

# Parameters
- `c`: Target coordinate system or layout cell.
- `width`: Width of the limbs of the L-shape.
- `size`: Outer extent (height and width) of the L-marker bounding box.
- `meta`: Metadata specifying layer and datatype information.

# Description
The ell marker consists of two perpendicular rectangular bars forming an "L" shape, with the inner corner located at the origin `(0, 0)`. The horizontal limb extends rightward and the vertical limb extends upward from the origin. Each bar has thickness `width` and length `size`.

# ASCII Diagram

    _   |<----size---->|
    |   |▓▓|
    |   |▓▓|
   size |▓▓|
    |   |<>| width
    |   |▓▓|
    v   |▓▓▓▓▓▓▓▓▓▓▓▓▓▓| width
        ^ (x,y) = (0,0) at the lower-left corner of the L.

"""
function ell_marker!(c::AbstractCoordinateSystem{T}, pitch::S, size::S, meta::Meta) where {T, S <: Coordinate}
    render!(c, simple_ell(width, size), meta)
end

end # end of module