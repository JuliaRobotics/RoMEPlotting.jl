using Base: Test
using KernelDensityEstimatePlotting
using IncrementalInference
using RoME


println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
include(joinpath(dirname(@__FILE__),"..","..","IncrementalInference","test","fourdoortest.jl"))
println("Success")


println("[TEST] plot functions...")
using Gadfly
# draw all beliefs
DOYTICKS = false
xx,ll = ls(fg)
msgPlots = drawHorBeliefsList(fg, xx, gt=gt,nhor=2);
evalstr = ""
for i in 1:length(msgPlots)
    evalstr = string(evalstr, ",msgPlots[$(i)]")
end
pl = eval(parse(string("vstack(",evalstr[2:end],")")));
println("Success")


# warn("plotMCMC needs ImageMagick on osx, not running test yet.")
# plotMCMC(tree, :x1, show=false)
# println("Success")


# Tests from RoME

# drawPosesLandms(fg);
