# test Pose2DPoint2D constraint evaluation function

using LinearAlgebra, Statistics
using RoME
# using KernelDensityEstimate
using RoMEPlotting, Gadfly
using Test

# to export PDF
using Gadfly
using Cairo

##

@testset "Prepare a 2D factor graph with poses and points..." begin

##

N = 100
fg = initfg()


initCov = Matrix(Diagonal([0.03;0.03;0.001]))
odoCov = Matrix(Diagonal([3.0;3.0;0.01]))

# Some starting position
addVariable!(fg, :x0, Pose2, N=N)
initPosePrior = PriorPose2(MvNormal(zeros(3), initCov))
addFactor!(fg,[:x0], initPosePrior)

# and a second pose
addVariable!(fg, :x1, Pose2, N=N)
ppc = Pose2Pose2(MvNormal([50.0;0.0;pi/2], odoCov))
addFactor!(fg, [:x0; :x1], ppc)

# test evaluation of pose pose constraint
pts = approxConv(fg, :x0x1f1, :x1)


# perform inference
tree = solveTree!(fg)


# check that yaw is working
addVariable!(fg, :x2, Pose2, N=N)
ppc = Pose2Pose2(MvNormal([50.0;0.0;0.0], odoCov))
addFactor!(fg, [:x1;:x2], ppc)


# new landmark
l1 = addVariable!(fg, :l1, Point2, N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2Point2BearingRange{Uniform, Normal}(Uniform(-pi,pi),Normal(rhoZ1,1.0))
addFactor!(fg, [:x0;:l1], ppr)


# add a prior to landmark
pp2 = PriorPoint2(MvNormal([10.0;0.0], Matrix(Diagonal([1.0;1.0]))))

f5 = addFactor!(fg,[:l1], pp2)

# do inference
tree = solveTree!(fg)

println("test Pose2D plotting")

plotSLAM2DPoses(fg);
plotSLAM2DLandmarks(fg);
plotSLAM2D(fg);

p1 = getBelief(fg, :l1)
# pts = getVal(fg, :l1)
# p1= kde!(pts)
p1c = getBelief(fg, :x0)

plotBelief( p1 , dimLbls=["x";"y";"z"]) # |> PDF("/tmp/test.pdf")


plotBelief( [p1c;p1] , dimLbls=["x";"y";"z"],c=["red";"black"],levels=3, dims=[1;2])
# plotBelief( [marginal(p1c,[1;2]);marginal(p1,[1;2])] , dimLbls=["x";"y";"z"],c=["red";"black"],levels=3, dims=[1;2])
p1c = deepcopy(p1)

plotBelief( marginal(getBelief(getVariable(fg, :x2)),[1;2]) , dimLbls=["x";"y";"z"])

axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']

# @warn "Reinsert draw test.pdf"
plotBelief( p1, dimLbls=["x";"y";"z"], axis=axis) |> PDF("/tmp/test.pdf",30cm,20cm)

#
Base.rm("/tmp/test.pdf")

##

end




#
