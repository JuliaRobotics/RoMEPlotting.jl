# test plotSLAM2D and plotPose features

using RoME, RoMEPlotting
using Test

##
@testset "test plotSLAM2D features" begin
##

fg = generateGraph_Hexagonal()

plotSLAM2D(fg, drawhist=false, drawPoints=false, drawContour=false, drawEllipse=false)
plotSLAM2D(fg, drawhist=false, drawPoints=false, drawContour=false, drawEllipse=true)
plotSLAM2D(fg, drawhist=false, drawPoints=false, drawContour=true, drawEllipse=false)
plotSLAM2D(fg, drawhist=false, drawPoints=false, drawContour=true, drawEllipse=true)
plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=false, drawEllipse=false)
plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=false, drawEllipse=true)
plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=true, drawEllipse=false)
plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=true, drawEllipse=true)

plotSLAM2D(fg, drawhist=true, drawPoints=false, drawContour=false, drawEllipse=false)
plotSLAM2D(fg, drawhist=true, drawPoints=false, drawContour=false, drawEllipse=true)
plotSLAM2D(fg, drawhist=true, drawPoints=false, drawContour=true, drawEllipse=false)
plotSLAM2D(fg, drawhist=true, drawPoints=false, drawContour=true, drawEllipse=true)
plotSLAM2D(fg, drawhist=true, drawPoints=true, drawContour=false, drawEllipse=false)
plotSLAM2D(fg, drawhist=true, drawPoints=true, drawContour=false, drawEllipse=true)
plotSLAM2D(fg, drawhist=true, drawPoints=true, drawContour=true, drawEllipse=false)
plotSLAM2D(fg, drawhist=true, drawPoints=true, drawContour=true, drawEllipse=true)

plotSLAM2D(fg; variableList=[:x1;:x5;:l1])

# plot runs but doesnt draw the levels at time of wrting the test.  Good to test the API regardless.
plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=true, drawEllipse=false, levels=5)


## test with save and load of factor graph

saveDFG("/tmp/rp_test", fg)
fg_ = loadDFG("/tmp/rp_test.tar.gz")

plotSLAM2D(fg, drawhist=false, drawPoints=true, drawContour=true, drawEllipse=false, levels=5)
# plotSLAM2D(fg_, drawhist=true, drawPoints=true, drawContour=true, drawEllipse=true)

# gets deleted lower down
# Base.rm("/tmp/rp_test.tar.gz")

##
end


@testset "test CJL Docs example with plotSLAM2DPoses" begin
##

fg = generateGraph_Hexagonal()

# Draw the (x,y) marginal estimated belief contour for :x0, :x2, and Lx4
pl = plotKDE(fg, [:x0; :x2; :x4], c=["red";"green";"blue"], levels=2, dims=[1;2])

# add a few fun layers
pl3 = plotSLAM2DPoses(fg, regexPoses=r"x\d", from=3, to=3, drawContour=false, drawEllipse=true)
pl5 = plotSLAM2DPoses(fg, regexPoses=r"x\d", from=5, to=5, drawContour=false, drawEllipse=true, drawPoints=false)
pl_ = plotSLAM2DPoses(fg, drawContour=false, drawPoints=false, dyadScale=0.001, to=5)
union!(pl.layers, pl3.layers)
union!(pl.layers, pl5.layers)
union!(pl.layers, pl_.layers)

# change the plotting coordinates
pl.coord = Coord.Cartesian(xmin=-10,xmax=20, ymin=-1, ymax=25)

# save the plot to SVG and giving dedicated (although optional) sizing
pl |> SVG("/tmp/test.svg", 25cm, 15cm)

Base.rm("/tmp/test.svg")

##
end



##
@testset "test plotPose features" begin
##

fg = generateGraph_Hexagonal()

X1 = getBelief(fg, :x1)

plotPose(Pose2, X1)

plotPose(fg, :x1)
plotPose(fg, :x1, levels=1, scale=0.3)
plotPose(fg, [:x1;:x2], levels=2, c=["green"; "magenta"])

@warn "fix plotPose legend keyword"
# plotPose(fg, [:x1;:x2], levels=2, c=["green"; "magenta"], legend=["hello"; "world"])


##

fg_ = loadDFG("/tmp/rp_test.tar.gz")
Base.rm("/tmp/rp_test.tar.gz")

plotPose(fg, :x1)


##
end