


# function predCurrFactorBeliefs(fgl::G,
#                                fc::Graphs.ExVertex  ) where G <: AbstractDFG
#   #
#   # TODO update to use ls and lsv functions
#   prjcurvals = Dict{String, Array{BallTreeDensity,1}}()
#   for v in getNeighbors(fgl, fc)
#     pred = kde!(evalFactor(fgl, fc, v.index))
#     curr = kde!(getVal(v))
#     prjcurvals[v.attributes["label"]] = [curr; pred]
#   end
#   return prjcurvals, collect(keys(prjcurvals))
# end

# function plotFactorBeliefs(fgl::G,
#                            flbl::Symbol ) where G <: AbstractDFG
#   #
#   if !haskey(fgl.fIDs, flbl)
#     println("no key $(flbl)")
#     return nothing
#   end
#   # for fc in fgl.f
#     # if fc[2].attributes["label"] == flbl

#     fc = fgl.g.vertices[fgl.fIDs[flbl]]  # fc = fgl.f[fgl.fIDs[flbl]]
#       prjcurvals, lbls = predCurrFactorBeliefs(fgl, fc)
#       if length(lbls) == 3
#         return vstack(
#         plotKDE(prjcurvals[lbls[1]]),
#         plotKDE(prjcurvals[lbls[2]]),
#         plotKDE(prjcurvals[lbls[3]]),
#         )
#       elseif length(lbls) == 2
#         return vstack(
#         plotKDE(prjcurvals[lbls[1]]),
#         plotKDE(prjcurvals[lbls[2]]),
#         )
#       elseif length(lbls) == 1
#         return plotKDE(prjcurvals[lbls[1]])
#       end

#     # end
#   # end
#   nothing
# end
# plotFactorBeliefs(fgl::G, flbl::T) where {G <: AbstractDFG, T <: AbstractString} = plotFactorBeliefs(fgl, Symbol(flbl))


# function whosWith(cliq::Graphs.ExVertex)
#   println("$(cliq.attributes["label"])")
#   for pot in cliq.attributes["data"].potentials
#       println("$(pot)")
#   end
#   nothing
# end

# function whosWith(bt::BayesTree, frt::String)
#     whosWith(whichCliq(bt,frt))
# end


# function plotUpMsgAtCliq(fg::G,
#                          cliq::Graphs.ExVertex  ) where G <: AbstractDFG
#   #
#   for id in keys(cliq.attributes["data"].debug.outmsg.p)
#       print("$(getVariable(fg,id).label), ") #fg.v[id].
#   end
#   println("")
#   sleep(0.1)
#   potens = getData(cliq).debug.outmsg.p
#   vArrPotentials(potens)
# end

# function plotUpMsgAtCliq(fg::G,
#                          bt::BayesTree,
#                          lbl::Union{String,Symbol}  ) where G <: AbstractDFG
#   #
#   plotUpMsgAtCliq(fg, whichCliq(bt, Symbol(lbl)) )
# end


# function dwnMsgsAtCliq(fg::G,
#                        cliq::Graphs.ExVertex  ) where G <: AbstractDFG
#   #
#   for id in keys(cliq.attributes["data"].debugDwn.outmsg.p)
#       print("$(getVariable(fg,id).label), ") # fg.v[id].
#   end
#   println("")
#   sleep(0.1)
#   potens = cliq.attributes["data"].debugDwn.outmsg.p
#   potens
# end

# function dwnMsgsAtCliq(fg::G,
#                        bt::BayesTree,
#                        lbl::Union{String,Symbol}  ) where G <: AbstractDFG
#   #
#   dwnMsgsAtCliq(fg, whichCliq(bt, Symbol(lbl)) )
# end






## =============================================================================
## Delete end v0.3.x
## =============================================================================

export drawPoses, drawLandms, drawPosesLandms
export drawSubmaps, drawMarginalContour, drawFeatTrackers, investigatePoseKDE
export drawHorBeliefsList, drawUpMsgAtCliq, drawFactorBeliefs
export drawPose2DMC!


# """
#     $(SIGNATURES)
#
# Standardize the length colors used by RoMEPlotting.
#
# Notes
# - Duplicated in KernelDensityEstimatePlotting
# """
# function getColorsByLength(len::Int=10)::Vector{String}
#   COLORS = String["red";"green";"blue";"black";"deepskyblue";"yellow";"magenta"]
#   morecyan = ["cyan" for i in (length(COLORS)+1):len]
#   retc = [COLORS; morecyan]
#   return retc[1:len]
# end


@deprecate drawPoses(x...; kwargs...) plotSLAM2DPoses(x...; kwargs...)
@deprecate drawLandms(x...; kwargs...) plotSLAM2DLandmarks(x...; kwargs...)
@deprecate drawPosesLandms(x...; kwargs...) plotSLAM2D(x...; kwargs...)
@deprecate drawSubmaps(x...; kwargs...) plotSLAM2DSubmaps(x...; kwargs...)
@deprecate drawMarginalContour(x...; kwargs...) plotMarginalContour(x...; kwargs...)
@deprecate drawFeatTrackers(x...; kwargs...) plotFeatTrackers(x...; kwargs...)
@deprecate drawHorBeliefsList(x...; kwargs...) plotHorBeliefsList(x...; kwargs...)
@deprecate drawUpMsgAtCliq(x...; kwargs...) plotUpMsgAtCliq(x...; kwargs...)
@deprecate drawFactorBeliefs(x...; kwargs...) plotFactorBeliefs(x...; kwargs...)
@deprecate drawPose2DMC!(x...; kwargs...) plotPose2DMC!(x...; kwargs...)


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

@deprecate investigatePoseKDE(x...; kwargs...) plotPose(x...; kwargs...)



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

@deprecate drawOneMC!(x...; kwargs...) plotOneMC!(x...; kwargs...)
@deprecate drawMCMCDebug(x...; kwargs...) plotMCMCDebug(x...; kwargs...)
@deprecate drawTreeUpwardMsgs(x...; kwargs...) plotTreeUpwardMsgs(x...; kwargs...)
@deprecate drawFrontalDens(x...; kwargs...) plotFrontalDens(x...; kwargs...)
@deprecate drawUpMCMCPose2D!(x...; kwargs...) plotUpMCMCPose2D!(x...; kwargs...)
@deprecate drawDwnMCMCPose2D!(x...; kwargs...) plotDwnMCMCPose2D!(x...; kwargs...)
@deprecate drawLbl(x...; kwargs...) plotLbl(x...; kwargs...)
@deprecate drawAnalysis(x...; kwargs...) plotAnalysis(x...; kwargs...)
@deprecate drawAllPose2DBeliefs(x...; kwargs...) plotAllPose2DBeliefs(x...; kwargs...)
@deprecate drawComicStripLM(x...; kwargs...) plotComicStripLM(x...; kwargs...)
@deprecate drawComicStrip(x...; kwargs...) plotComicStrip(x...; kwargs...)



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
#   verts = map((x)->getBelief(getVariable(fgl, x)), vsym)
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
## =============================================================================
## Old Functions -- Possibly Deprecated
## =============================================================================
function plotHorBeliefsList(fgl::G,
                            lbls::Array{Symbol,1};
                            nhor::Int=-1,
                            gt=nothing,
                            N::Int=200,
                            extend=0.1  ) where G <: AbstractDFG
  #
  len = length(lbls)
  pDens = BallTreeDensity[]
  for lb in lbls
    ptkde = getBelief(getVariable(fgl,lb))
    push!(pDens, ptkde )
  end

  if nhor<1
    nhor = round(Int,sqrt(len))
  end
  vlen = ceil(Int, len/nhor)
  vv = Array{Gadfly.Compose.Context,1}(undef, vlen)
  conslb = deepcopy(lbls)
  vidx = 0
  for i in 1:nhor:len
    pH = BallTreeDensity[]
    gtvals = Dict{Int,Array{Float64,2}}()
    labels = String[]
    for j in 0:(nhor-1)
      if i+j <= len
        push!(pH, pDens[i+j])
        push!(labels, string(lbls[i+j]))
        if gt != nothing gtvals[j+1] = gt[lbls[i+j]]  end
      end
    end
    vidx+=1
    if gt !=nothing
      vv[vidx] = KernelDensityEstimatePlotting.drawHorDens(pH, N=N, gt=gtvals, lbls=labels, extend=extend)
    else
      vv[vidx] = KernelDensityEstimatePlotting.drawHorDens(pH, N=N, lbls=labels, extend=extend)
    end
  end
  vv
end
