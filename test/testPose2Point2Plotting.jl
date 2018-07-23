# test Pose2DPoint2D constraint evaluation function

using RoME, IncrementalInference, Distributions
using KernelDensityEstimate
using RoMEPlotting, Gadfly
using Base.Test

begin

println("Prepare a 2D factor graph with poses and points...")

N = 75
fg = initfg()


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, :x0, Pose2, N=N)
# v1 = addNode!(fg, :x0, zeros(3,1), diagm([1.0;1.0;0.1]), N=N)
initPosePrior = PriorPose2(zeros(3,1), initCov, [1.0])
f1  = addFactor!(fg,[v1], initPosePrior)

# and a second pose
v2 = addNode!(fg, :x1, Pose2, N=N)
# v2 = addNode!(fg, :x1, vectoarr2([50.0;0.0;pi/2]), diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;pi/2]), odoCov, [1.0])
f2 = addFactor!(fg, [v1;v2], ppc)

# test evaluation of pose pose constraint
pts = evalFactor2(fg, f2, v2.index)

# @show ls(fg)

tree = wipeBuildNewTree!(fg)
inferOverTreeR!(fg, tree,N=N)
# inferOverTree!(fg, tree, N=N)

# check that yaw is working
v3 = addNode!(fg, :x2, zeros(3,1), diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;0.0]), odoCov, [1.0])
f3 = addFactor!(fg, [v2;v3], ppc)


# new landmark
l1 = addNode!(fg, :l1, zeros(2,1), diagm([1.0;1.0]), N=N)
# and pose to landmark constraint
rhoZ1 = norm([10.0;0.0])
ppr = Pose2Point2BearingRange{Uniform, Normal}(Uniform(-pi,pi),Normal(rhoZ1,1.0))
f4 = addFactor!(fg, [v1;l1], ppr)


# add a prior to landmark
pp2 = PriorPoint2D([10.0;0.0], diagm([1.0;1.0]), [1.0])

f5 = addFactor!(fg,[l1], pp2)

ensureAllInitialized!(fg)
tree = wipeBuildNewTree!(fg)
[inferOverTree!(fg, tree, N=N) for i in 1:2]

println("test Pose2D plotting")

drawPoses(fg);
drawPosesLandms(fg);

pts = getVal(fg, :l1)

p1= kde!(pts)
p1c = getVertKDE(getVert(fg, :x0))
plotKDE( p1 , dimLbls=["x";"y";"z"])

plotKDE( [marginal(p1c,[1;2]);marginal(p1,[1;2])] , dimLbls=["x";"y";"z"],c=["red";"black"],levels=3)
p1c = deepcopy(p1)

plotKDE( marginal(getVertKDE(fg, :x2),[1;2]) , dimLbls=["x";"y";"z"])

axis = [[1.5;3.5]';[-1.25;1.25]';[-1.0;1.0]']
warn("Reinsert draw test.pdf")
# Gadfly.draw( PDF("test.pdf",30cm,20cm),
#       plotKDE( p1, dimLbls=["x";"y";"z"], axis=axis)  )
# #
# Base.rm("test.pdf")

end




#
