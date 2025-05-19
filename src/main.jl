

include("SNUDeviceComponents.jl")

using DeviceLayout
import DeviceLayout: μm, nm
using .SNUDeviceComponents

c = Cell("pattern", nm)
# Logos.snu_logo!(c, 200μm, meta = GDSMeta(0, 0))


steps = 10
foot_width = 4.0μm
foot_height = 10.0μm
span = 40.0μm
metas = GDSMeta.(range(0, 10))

Components.grayscale_bridge!(c, steps, foot_width, foot_height, span, metas)