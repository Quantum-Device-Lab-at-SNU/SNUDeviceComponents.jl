
using DeviceLayout
import DeviceLayout: μm, nm
using SNUDeviceComponents

c = Cell("meander_inductor_test", nm)
Components.LumpedElement.meander_inductor!(
    c, 80.0μm, 2.0μm, 5.0μm, 21, 50.0μm, 2.0μm, 20.0μm, GDSMeta(0);
    first_turn=:left, clearance_corner_r=10.0μm)

