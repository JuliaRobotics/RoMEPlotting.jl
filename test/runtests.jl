using Test
using KernelDensityEstimatePlotting
using IncrementalInference
using RoME
using RoMEPlotting

##


include("testPose2Point2Plotting.jl")
include("testPlotKDE.jl")
include("testPlotSLAMandPose.jl")

# @warn "plotMCMC needs ImageMagick on osx, not running test yet."
# plotMCMC(tree, :x1, show=false)
# println("Success")



# test ellipse function
# Z = MvNormal([0.0;0],[9.0 5; 5 4])
# eX = covEllipseParameterized(Z)


## legacy


# println("[TEST] with local Graphs.jl dictionary and arrays only (multicore)...")
# include(joinpath(dirname(@__FILE__),"..","..","IncrementalInference","test","fourdoortest.jl"))
# println("Success")
#
# println("[TEST] plot functions...")
# using Gadfly
# # draw all beliefs
# DOYTICKS = false
# xx,ll = ls(fg)
# msgPlots = drawHorBeliefsList(fg, xx, gt=gt,nhor=2);
# pl = vstack(msgPlots...);
# # Gadfly.draw(PDF("/tmp/test.pdf",15cm,30cm),pl)
# println("Success")

