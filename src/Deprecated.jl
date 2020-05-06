
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

export drawOneMC!
export drawMCMCDebug
export drawTreeUpwardMsgs
export drawFrontalDens
export drawUpMCMCPose2D!
export drawDwnMCMCPose2D!
export drawLbl
export drawAnalysis
export drawAllPose2DBeliefs
export drawComicStripLM
export drawComicStrip

@deprecate drawOneMC!(x...) plotOneMC!(x...)
@deprecate drawMCMCDebug(x...) plotMCMCDebug(x...)
@deprecate drawTreeUpwardMsgs(x...) plotTreeUpwardMsgs(x...)
@deprecate drawFrontalDens(x...) plotFrontalDens(x...)
@deprecate drawUpMCMCPose2D!(x...) plotUpMCMCPose2D!(x...)
@deprecate drawDwnMCMCPose2D!(x...) plotDwnMCMCPose2D!(x...)
@deprecate drawLbl(x...) plotLbl(x...)
@deprecate drawAnalysis(x...) plotAnalysis(x...)
@deprecate drawAllPose2DBeliefs(x...) plotAllPose2DBeliefs(x...)
@deprecate drawComicStripLM(x...) plotComicStripLM(x...)
@deprecate drawComicStrip(x...) plotComicStrip(x...)



# function plotKDE(fgl::FactorGraph,
#                  vsym::Vector{Symbol};
#                  axis=nothing,
#                  dims=nothing,
#                  c=getColorsByLength(length(vsym)),
#                  levels::Int=4,
#                  title::Union{Nothing, T}=nothing,
#                  overlay=nothing  ) where {T <: AbstractString}
#   #
#   @warn "plotKDE for FactorGraph is deprecated, use DistributedFactorGraphs objects instead."
#   verts = map((x)->getKDE(getVariable(fgl, x)), vsym)
#   plotKDE(verts, dims=dims, c=c, axis=axis, levels=levels, title=title, overlay=overlay )
# end
#
# function plotKDE(fgl::FactorGraph,
#                  vsym::Symbol;
#                  axis=nothing,
#                  dims=nothing,
#                  c=nothing,
#                  levels=4,
#                  title::Union{Nothing, T}=nothing) where {T <: AbstractString}
#   #
#   @warn "plotKDE for FactorGraph is deprecated, use DistributedFactorGraphs objects instead."
#   plotKDE(fgl, Symbol[vsym;], dims=dims, c=c, axis=axis, levels=levels, title=title)
# end



#
