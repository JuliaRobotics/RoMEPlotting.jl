module RoMEPlotting

using IncrementalInference, RoME
using Gadfly

export

  # Associated with IncrementalInference
  investigateMultidimKDE,
  drawHorDens,
  drawHorBeliefsList,
  spyCliqMat,
  plotKDEofnc,
  plotKDEresiduals,
  plotMCMC,
  plotUpMsgsAtCliq,
  plotPriorsAtCliq,
  investigateMultidimKDE,
  draw,
  whosWith,
  drawUpMsgAtCliq,
  dwnMsgsAtCliq,
  drawPose2DMC!,
  mcmcPose2D!,
  # drawUpMCMCPose2D!,
  # drawDwnMCMCPose2D!,
  drawLbl,
  predCurrFactorBeliefs,
  drawHorDens,
  drawHorBeliefsList,
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
  plotTrckStep


include("SolverVisualization.jl")

include("RobotViz.jl")



end
