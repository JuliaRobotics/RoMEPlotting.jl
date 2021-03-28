module RoMEPlotting

using Cairo, Fontconfig

using Reexport
@reexport using Gadfly
@reexport using Colors

using Statistics, LinearAlgebra
using StatsBase
using Compose
using Dates
using DistributedFactorGraphs
using KernelDensityEstimate, KernelDensityEstimatePlotting
using IncrementalInference, RoME
using DocStringExtensions
using ApproxManifoldProducts
using Requires
# import ApproxManifoldProducts: mmd! # future dependency

import Gadfly: plot
import KernelDensityEstimatePlotting: plot, drawHorDens, plotKDE
import KernelDensityEstimatePlotting: getColorsByLength


# assuming this is a good size for everybody
@info "Assuming plot size, Gadfly.set_default_plot_size(30cm,20cm)"
Gadfly.set_default_plot_size(30cm,20cm)

export
  # Associated with IncrementalInference
  investigateMultidimKDE,
  drawHorDens, # from KDEPlotting
  plotKDE,
  plotVariable2D,
  plotKDEofnc,
  plotKDEresiduals,
  plotMCMC,
  plotUpMsgsAtCliq,
  plotPriorsAtCliq,
  investigateMultidimKDE,
  plot,
  whosWith,
  plotUpMsgAtCliq,
  dwnMsgsAtCliq,
  plotPose2DMC!,
  mcmcPose2D!,
  plotLbl,
  predCurrFactorBeliefs,
  plotFactorBeliefs,
  localProduct,
  plotLocalProduct,
  plotLocalProductCylinder,
  plotTreeProductUp,
  plotTreeProductDown,
  plotCliqUpMsgs,
  saveplot,
  animateVertexBelief,
  getColorsByLength,

  # Associated with RoME
  plotTrajectoryArrayPose2,
  togglePrtStbLines,
  plotLsrScanFeats,
  saveImgSeq,
  stbPrtLineLayers!,
  plotSLAM2D,
  plotSLAM2DPoses,
  plotSLAM2DLandmarks,
  plotPose,
  plotMarginalContour,
  accumulateMarginalContours,
  plotPose3Pairs,
  progressExamplePlot,
  plotTrckStep,
  plotPose2Vels,
  plotProductVsKDE,
  plotPairVariables,
  plotPairPose2,
  plotVariableGivenFactor,
  plotCliqDownMsgs,
  plotFactor,
  plotFactorValues,
  plotFactorMeasurements,
  reportFactors


# EXPERIMENTAL
const AbstractMatrix__{T} = Union{AbstractArray{T,2}, Adjoint{T,<:AbstractArray{T,2}}}


using Requires

# will be overwritten if flux is present (dont make const)
PlotTypesPose2 = Union{Type{Pose2Pose2}, Type{Pose2Point2BearingRange}, Type{Pose2Point2Range}, Type{Pose2Point2Bearing}}
ExtendedPose2Pose2Types = Pose2Pose2


include("SolverVisualization.jl")
include("RobotViz.jl")
include("CurrentlyUnmaintained.jl")
include("PlotHexUtils.jl")
include("PlotFactors.jl")
include("PlotVariables.jl")
include("PlotFactorsReload.jl")
include("Deprecated.jl")
include("NeedsFixing.jl")

function __init__()
  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" include("FluxSpecificFeatures.jl")
  @require ImageMagick="6218d12a-5da1-5696-b52f-db25d2ecc6d1" include("imagemagick.jl")
end

end
