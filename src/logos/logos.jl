module Logos

using DeviceLayout
import DeviceLayout: Coordinate, AbstractCoordinateSystem
import DeviceLayout: Meta, GDSMeta, μm
import FileIO: load

export snu_logo!

"""
    snu_logo!(c::AbstractCoordinateSystem{T}, w::S, meta:: Meta) where {T, S <: Coordinate}

Renders polygons to cell `c` to construct a Seoul National University logo pattern.

 - `w`: horizontal width of the desired logo
"""
function snu_logo!(c::AbstractCoordinateSystem{T}, w::S, meta:: Meta) where {T, S <: Coordinate}
    # fancy snu logo
    snu_pattern = load(joinpath(@__DIR__, "snu_pattern.gds"))["snu_pattern"]

    bounding_rectangle = bounds(snu_pattern)
    w0 = width(bounding_rectangle)          # width of the original pattern in the GDS file

    mag_factor = w / w0                     # resize the pattern to the specified width
    snu_pattern -= center(bounding_rectangle)
    snu_pattern *= mag_factor
    render!(c, snu_pattern.elements, meta)
end

end # end of module