module Components

import DeviceLayout: AbstractCoordinateSystem, Coordinate, Meta
using DeviceLayout

export grayscale_bridge!, cross_marker!, ell_marker!

include("alignment_markers.jl")

"""
    grayscale_bridge!(
    c::AbstractCoordinateSystem{T},
    steps::Integer,
    foot_width::S,
    foot_height::S,
    span::S,
    metas::AbstractVector{<:Meta}
    ) where {T,S<:Coordinate}

Renders polygons to cell `c` to construct a bridge based on grayscale lithography.
Can be called twice to construct a bridge with a hopover.
- `steps`: Number of discrete steps in the bridge height.
- `foot_width`: Width of the foot.
- `foot_height`: Height of the bridge foot (in the plane).
- `span`: Length of the bridge span.
- `metas`: A vector of `Meta` objects with length `steps + 1`. The first object is the
  `Meta` for the bridge foot and the others are for the bridge steps, in order.

```
    <-f.w.->                          <-f.w.->
    ******************************************  ---
    |+++++|                            |+++++|   |
    |+++++|                            |+++++|   |
    |+++++|                            |+++++|   |
    |+++++|<-----------span----------->|+++++|   | foot height
    |+++++|                            |+++++|   |
    |+++++|                            |+++++|   |
    |+++++|                            |+++++|   |
    ******************************************  ---
```

The profile of the bridge is given as `h(x) = H*(2 - cosh(2*cosh^(-1)(2) x/L))` (inverted
catenary arch) where `H` is the height of bridge and `L` is the total span of bridge.
The inverse of h(x) is taken to compute the size of each layer.
"""
function grayscale_bridge!(
    c::AbstractCoordinateSystem{T},
    steps::Integer,
    foot_width::S,
    foot_height::S,
    span::S,
    metas::AbstractVector{<:Meta}
    ) where {T,S<:Coordinate}
    @assert length(metas) == steps + 1

    # Bridge Profile
    yinv_gnd(y, span) = span/2 * acosh(-y + 2) / acosh(2)
    ypos = range(0, stop=1, length=steps + 1)

    # render bridge feet layer
    render!(
        c, Rectangle(Point(-(span/2 + foot_width), -foot_height/2),
        Point(-span/2, foot_height/2)), metas[1]
    )
    render!(
        c, Rectangle(Point(span/2, -foot_height/2),
        Point(span/2 + foot_width, foot_height/2)), metas[1]
    )

    x = zero(T)
    for (i, val) in enumerate(ypos)
        if i > 1
            xpos = yinv_gnd(val, span)
            render!(c, Rectangle(Point(x, -foot_height/2),
                Point(xpos, foot_height/2)), metas[i])
            render!(c, Rectangle(Point(-x, -foot_height/2),
                Point(-xpos, foot_height/2)), metas[i])
            x = xpos
        else
            x = yinv_gnd(val, span)
        end
    end
    c
end

end # end of module