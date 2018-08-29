module RoMEPlotting

using Graphs
using KernelDensityEstimate, KernelDensityEstimatePlotting
using IncrementalInference, RoME
using Gadfly
using Colors
using DocStringExtensions

import KernelDensityEstimatePlotting: plot, drawHorDens, plotKDE
import IncrementalInference: CliqGibbsMC, DebugCliqMCMC
import Graphs: plot
import Gadfly: plot

export
  # Associated with IncrementalInference
  investigateMultidimKDE,
  drawHorDens,
  drawHorBeliefsList,
  spyCliqMat,
  plotKDE,
  plotKDEofnc,
  plotKDEresiduals,
  plotMCMC,
  plotKDE,
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
  saveplot,
  animateVertexBelief,
  getColorsByLength,

  # Associated with RoME
  togglePrtStbLines,
  plotLsrScanFeats,
  drawFeatTrackers,
  saveImgSeq,
  stbPrtLineLayers!,
  drawPoses,
  drawLandms,
  drawPosesLandms,
  drawSubmaps,
  investigatePoseKDE,
  drawMarginalContour,
  accumulateMarginalContours,
  plotPose3Pairs,
  progressExamplePlot,
  plotTrckStep,
  plotPose2Vels


include("SolverVisualization.jl")

include("RobotViz.jl")



end
