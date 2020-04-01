module RoMEPlotting

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
# import ApproxManifoldProducts: mmd! # future dependency


import KernelDensityEstimatePlotting: plot, drawHorDens, plotKDE
import IncrementalInference: CliqGibbsMC, DebugCliqMCMC
import Gadfly: plot


export
  # Associated with IncrementalInference
  investigateMultidimKDE,
  drawHorDens,
  drawHorBeliefsList,
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
  drawUpMsgAtCliq,
  dwnMsgsAtCliq,
  drawPose2DMC!,
  mcmcPose2D!,
  # drawUpMCMCPose2D!,
  # drawDwnMCMCPose2D!,
  drawLbl,
  predCurrFactorBeliefs,
  drawFactorBeliefs,
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
  drawFeatTrackers,
  saveImgSeq,
  stbPrtLineLayers!,
  drawPoses,
  drawLandms,
  drawPosesLandms,
  drawSubmaps,
  investigatePoseKDE, # not sure, likely obsolete -- use plotPose instead
  plotPose,
  drawMarginalContour,
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


end
