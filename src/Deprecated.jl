
## =============================================================================
## Delete end v0.3.x
## =============================================================================

export drawPoses, drawLandms, drawPosesLandms
export drawSubmaps, drawMarginalContour, drawFeatTrackers, investigatePoseKDE
export drawHorBeliefsList, drawUpMsgAtCliq, drawFactorBeliefs
export drawPose2DMC!

@deprecate drawPoses(x...) plotSLAM2DPoses(x...)
@deprecate drawLandms(x...) plotSLAM2DLandmarks(x...)
@deprecate drawPosesLandms(x...) plotSLAM2D(x...)
@deprecate drawSubmaps(x...) plotSLAM2DSubmaps(x...)
@deprecate drawMarginalContour(x...) plotMarginalContour(x...)
@deprecate drawFeatTrackers(x...) plotFeatTrackers(x...)
@deprecate drawHorBeliefsList(x...) plotHorBeliefsList(x...)
@deprecate drawUpMsgAtCliq(x...) plotUpMsgAtCliq(x...)
@deprecate drawFactorBeliefs(x...) plotFactorBeliefs(x...)
@deprecate drawPose2DMC!(x...) plotPose2DMC!(x...)


# deprecated
function investigatePoseKDE(p::BallTreeDensity, p0::BallTreeDensity)
    return investigateMultidimKDE(p, p0)
end


function investigatePoseKDE(p::Array{BallTreeDensity,1})
    return investigateMultidimKDE(p)
end

function investigatePoseKDE(p::BallTreeDensity)
    return investigateMultidimKDE(p)
end

@deprecate investigatePoseKDE(x...) plotPose(x...)



# function vArrPotentials(potens::Dict{Symbol,EasyMessage})
#   vv = Array{Gadfly.Compose.Context,1}(length(potens))
#   i = 0
#   oned=false
#   for p in potens
#       i+=1
#       pb = kde!(p[2].pts, p[2].bws)
#       if size(p[2].pts,1) > 3
#         # vv[i] = plotKDE(pb)
#         error("can't handle higher dimensional plots here yet")
#       elseif size(p[2].pts,1) > 1
#         vv[i] = investigateMultidimKDE(pb)
#       else
#         vv[i] = plotKDE(pb)
#       end
#   end
#   return vv
# end


# function draw(em::EasyMessage;xlbl="X")
#   p = Union{}
#   if size(em.pts,1) == 1
#     p=plotKDE(kde!(em),xlbl=xlbl)
#   else
#     p=plotKDE(kde!(em))
#   end
#   return p
# end

#
