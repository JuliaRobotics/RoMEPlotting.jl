# test multi kde plot, issue 65

using RoME, RoMEPlotting
using Test
# Gadfly.set_default_plot_size(35cm, 25cm)


@testset "test plotKDE array (issue 65)" begin


fg = initfg()

addVariable!(fg, :x0, ContinuousScalar)
addVariable!(fg, :x1, ContinuousScalar)

addFactor!(fg, [:x0], Prior(Normal()))
addFactor!(fg, [:x0;:x1], LinearRelative(Normal()))

ensureAllInitialized!(fg)

plotKDE(fg, ls(fg))


end
