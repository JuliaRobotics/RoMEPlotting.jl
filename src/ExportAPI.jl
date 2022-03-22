

## =====================================================================================
# Primary supported functionality
## =====================================================================================


export plotBelief
export plotPose
export plotVariable2D
export plotSLAM2D, plotSLAM2D_KeyAndRef, plotSLAM2DLandmarks, plotSLAM2DPoses
export getColorsByLength

# FIXME what to do about this?
export  plot

# needs better testing (might be broken following upgrade to Manifolds.jl)
export plotPose2Vels
export plotFactor
export plotFactorValues, plotFactorMeasurements, plotVariableGivenFactor
export plotPairVariables, plotPairPose2, plotPose3Pairs
export localProduct, plotLocalProductCylinder, plotLocalProduct
export saveplot
export plotMarginalContour, accumulateMarginalContours


# plot tree functionality (needs testing, possible broken with recent upgrades)
export plotUpMsgsAtCliq
export plotTreeProductUp, plotTreeProductDown, plotCliqUpMsgs, plotCliqDownMsgs


## =====================================================================================
## Known FIXME issues
## =====================================================================================

# export reportFactor
# export plotPriorsAtCliq


## =====================================================================================
## CURRENTLY UNMAINTAINED, utility functions that must be refactored, consolidated, standardized, or deprecated
## =====================================================================================


export  investigateMultidimKDE,
  plotKDEofnc,
  plotKDEresiduals,
  investigateMultidimKDE,
  plotLbl,
  animateVertexBelief,
  # Associated with RoME
  plotTrajectoryArrayPose2,
  saveImgSeq,
  plotProductVsKDE

  
## =====================================================================================
## Old Victoria Park dataset features that were good but unmaintained
## =====================================================================================

export plotLsrScanFeats, progressExamplePlot, plotTrckStep



#