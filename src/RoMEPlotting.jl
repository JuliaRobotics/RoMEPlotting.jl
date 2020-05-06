module RoMEPlotting

using Reexport
@reexport using Gadfly
@reexport using Colors

using Cairo, Fontconfig

using Statistics, LinearAlgebra
using StatsBase
using Compose
using Dates
using DistributedFactorGraphs
using KernelDensityEstimate, KernelDensityEstimatePlotting
using IncrementalInference, RoME
using DocStringExtensions
using ApproxManifoldProducts
# import ApproxManifoldProducts: mmd! # future dependency


import KernelDensityEstimatePlotting: plot, drawHorDens, plotKDE
import IncrementalInference: CliqGibbsMC, DebugCliqMCMC
import Gadfly: plot

# assuming this is a good size for everybody
@info "Assuming plot size, Gadfly.set_default_plot_size(30cm,20cm)"
Gadfly.set_default_plot_size(30cm,20cm)

export
  # Associated with IncrementalInference
  investigateMultidimKDE,
  drawHorDens, # from KDEPlotting
  spyCliqMat,
  plotKDE,
  plotVariable2D,
  plotKDEofnc,
  plotKDEresiduals,
  plotMCMC,
  plotUpMsgsAtCliq,
  plotPriorsAtCliq,
  investigateMultidimKDE,
  draw,
  plot,
  whosWith,
  plotUpMsgAtCliq,
  dwnMsgsAtCliq,
  plotPose2DMC!,
  mcmcPose2D!,
  drawLbl,
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
  plotFactorMeasurements,
  reportFactors


# EXPERIMENTAL
const AbstractMatrix__{T} = Union{AbstractArray{T,2}, Adjoint{T,<:AbstractArray{T,2}}}



include("SolverVisualization.jl")
include("RobotViz.jl")
include("PlotHexUtils.jl")
include("PlotFactors.jl")
include("Deprecated.jl")

end
