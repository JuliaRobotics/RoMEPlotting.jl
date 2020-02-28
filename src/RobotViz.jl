# new exports from this file, WIP

export covEllipseParameterized, plotCovEllipseLayer


## Src

global DISABLESTBPRTLINES = false

function togglePrtStbLines()
  global DISABLESTBPRTLINES
  DISABLESTBPRTLINES = !DISABLESTBPRTLINES
end

function plotLsrScanFeats(br::Array{Float64,2})
  Cart = zeros(size(br))
  Cart[:,1] = br[:,2].*cos(br[:,1])
  Cart[:,2] = br[:,2].*sin(br[:,1])
  plot(x=Cart[:,1],y=Cart[:,2],Geom.point,
  Guide.xticks(ticks=collect(-60:10:60)),
  Guide.yticks(ticks=collect(0:10:80)))
end

function drawFeatTrackers(trkrs::Dict{Int64,Feature}, bfts::Array{Float64,2})
  musX = Float64[]
  varX = Float64[]
  musY = Float64[]
  varY = Float64[]
  allPtsX = Float64[]
  allPtsY = Float64[]

  for ftr in trkrs
    pts = getPoints(ftr[2].bel)
    allPtsX = [allPtsX; vec(pts[1,:])]
    allPtsY = [allPtsY; vec(pts[2,:])]

    push!(musX, Statistics.mean(vec(pts[1,:])))
    push!(varX, Statistics.std(vec(pts[1,:])))
    push!(musY, Statistics.mean(vec(pts[2,:])))
    push!(varY, Statistics.std(vec(pts[2,:])))
  end

  X = Float64[]
  Y = Float64[]

  if size(bfts,2) > 0
    if bfts[1,1] != 0.0 && bfts[2,1] != 0.0 && bfts[3,1] != 0.0
      for i in 1:size(bfts,2)
          u, R = p2c(vec(bfts[:,i]))
          push!(X, u[1])
          push!(Y, u[2])
      end
    end
  end

  # Guide.yticks(ticks=collect(-60:10:60)),
  # Guide.xticks(ticks=collect(0:10:80))
  p = plot(layer(x=musX, y=musY, Geom.point, Theme(default_color=colorant"red")),
  layer(x=allPtsX, y=allPtsY, Geom.histogram2d),
  Guide.yticks(ticks=collect(-70:10:70)),
  Guide.xticks(ticks=collect(-40:10:80)))
  for i in 1:length(X)
    push!(p.layers, Gadfly.layer(x=[0.0;X[i]], y=[0.0;Y[i]], Geom.line, Gadfly.Theme(default_color=colorant"magenta"))[1])
  end
  p
end


function saveImgSeq(d::Dict{Int64,Array{Float64,2}}; from::Int=1,to::Int=10,step::Int=1)
  for i in from:step:to
    p = plotLsrScanFeats(lsrBR(d[i]));
    Gadfly.draw(PNG(string("imgs/img",i,".png"),25cm,25cm),p)
  end
  nothing
end


"""
    $SIGNATURES

Calculate the ellipse from covariance matrix and return as lambda function.

Notes
- https://cookierobotics.com/007/

Related

plotCovEllipseLayer
"""
function covEllipseParameterized(distr::MvNormal; meanOffset::Bool=true)
  println(round.(cov(distr), digits=3) )
  a = cov(distr)[1,1]
  b = cov(distr)[1,2] + cov(distr)[2,1]
  b *= 0.5
  c = cov(distr)[2,2]

  λ1 = 0.5*(a+c) + sqrt(0.25*(a-c)^2 + b^2)
  λ2 = 0.5*(a+c) - sqrt(0.25*(a-c)^2 + b^2)
  θ = if b ≈ 0 && a >= c
    0
  elseif b ≈ 0 && a < c
    pi/2
  else
    atan(λ1-a, b)
  end
  sλ1 = sqrt(λ1)
  sλ2 = sqrt(λ2)
  sθ =  sin(θ)
  cθ =  cos(θ)
  ox = meanOffset ? distr.μ[1] : 0
  oy = meanOffset ? distr.μ[2] : 0
  (t) -> [sλ1*cθ*cos(t)-sλ2*sθ*sin(t)+ox; sλ1*sθ*cos(t)+sλ2*cθ*sin(t)+oy]
end

covEllipseParameterized(pts::Array{Float64,2}; meanOffset::Bool=true) = covEllipseParameterized( fit(MvNormal, pts), meanOffset=meanOffset )
covEllipseParameterized(X::BallTreeDensity; meanOffset::Bool=true) = covEllipseParameterized( getPoints(X), meanOffset=meanOffset )
covEllipseParameterized(dfg::AbstractDFG, sym::Symbol; meanOffset::Bool=true) = covEllipseParameterized( getKDE(dfg, sym), meanOffset=meanOffset )


"""
    $SIGNATURES

Plotting tool to draw Gadfly layers of ellipses of 2D covariance fitted to the belief of factor graph variable nonparametric points.

Related

covEllipseParameterized, drawPosesLandms
"""
function plotCovEllipseLayer(dfg::AbstractDFG,
                             vsym::Symbol;
                             points::Bool=true,
                             ellipseColor::AbstractString="gray30",
                             pointsColor::AbstractString="gray30",
                             drawEllipse::Bool=true  )::Vector
  #
  PL = []

  # points to work from
  pp = getPoints(getKDE(dfg, vsym))

  if drawEllipse
    # get ellipse function
    eX = covEllipseParameterized(pp[1:2,:], meanOffset=false)
    vEl = eX.(0:0.02:2pi)
    el = [(x->x[1]).(vEl) (x->x[2]).(vEl)]
    # add suggested PPE mean offset
    el[:,1] .+= getVariablePPE(dfg, vsym).suggested[1]
    el[:,2] .+= getVariablePPE(dfg, vsym).suggested[2]

    # add the ellipse layers
    plelX2 = Gadfly.layer(x=el[:,1], y=el[:,2], Geom.path, Theme(default_color=parse(Colorant, ellipseColor)))
    push!(PL, plelX2[1])
  end

  # add the points layer if needed
  if points
    plelX1 = Gadfly.layer(x=pp[1,:], y=pp[2,:], Geom.point, Theme(default_color=parse(Colorant, pointsColor), point_size=1pt))
    push!(PL, plelX1[1])
  end

  return PL
end


"""
    $SIGNATURES

Plot trajectory of Array{,2} with rows as consecutive entries and columns as x,y,theta.
"""
function plotTrajectoryArrayPose2(arr::AbstractMatrix__{<:Real};
                                  spscale::Float64=0.5,
                                  triadStride::Int=50)
  #
  @assert size(arr,2)==3
  trajPlt = Gadfly.plot(x=arr[:,1], y=arr[:,2], Geom.path, Coord.cartesian(fixed=true, aspect_ratio=1))

  if triadStride != -1
    Xpp = arr[1:triadStride:end,1]
    Ypp = arr[1:triadStride:end,2]
    Thpp = arr[1:triadStride:end,3]
    addXYLineLayers!(trajPlt, Xpp, Ypp, Thpp, l=spscale)
  end
  return trajPlt
end



## TODO -- you were here with port starboard lines
function stbPrtLineLayers!(pl, Xpp, Ypp, Thpp; l::Float64=5.0)
    if DISABLESTBPRTLINES
      return nothing
    end


    lnstpr = [0.0;l;0.0]
    lnstpg = [0.0;-l;0.0]

    Rd  =SE2(lnstpr)
    Gr = SE2(lnstpg)

    for i in 1:length(Xpp)
      lnstt = [Xpp[i];Ypp[i];Thpp[i]]
      Ps = SE2(lnstt)
      lnr = se2vee(Ps*Rd)
      lng = se2vee(Ps*Gr)
      xsr = [Xpp[i];lnr[1]]
      ysr = [Ypp[i];lnr[2]]
      xsg = [Xpp[i];lng[1]]
      ysg = [Ypp[i];lng[2]]

      push!(pl.layers, layer(x=xsr, y=ysr, Geom.path(), Gadfly.Theme(default_color=colorant"red", line_width=1.5pt))[1] )
      push!(pl.layers, layer(x=xsg, y=ysg, Geom.path(), Gadfly.Theme(default_color=colorant"green", line_width=1.5pt))[1] )
    end
    nothing
end

# draw the reference frame as a red-green dyad
function addXYLineLayers!(pl, Xpp, Ypp, Thpp; l::Float64=1.0, manualColor::Union{Nothing, AbstractString}=nothing  )
    lnstpr = [l;0.0;0.0]
    lnstpg = [0.0;l;0.0]

    Rd  =SE2(lnstpr)
    Gr = SE2(lnstpg)

    for i in 1:length(Xpp)
      lnstt = [Xpp[i];Ypp[i];Thpp[i]]
      Ps = SE2(lnstt)
      lnr = se2vee(Ps*Rd)
      lng = se2vee(Ps*Gr)
      xsr = [Xpp[i];lnr[1]]
      ysr = [Ypp[i];lnr[2]]
      xsg = [Xpp[i];lng[1]]
      ysg = [Ypp[i];lng[2]]

      push!(pl.layers, layer(x=xsr, y=ysr, Geom.path(), Gadfly.Theme(default_color=parse(Colorant, manualColor == nothing ? "red" : manualColor), line_width=1.5pt))[1] )
      push!(pl.layers, layer(x=xsg, y=ysg, Geom.path(), Gadfly.Theme(default_color=parse(Colorant, manualColor == nothing ? "green" : manualColor), line_width=1.5pt))[1] )
    end
    nothing
end

# function lblsFromTo(from,to)
#   lbls=String[]
#   [push!(lbls, "$(i)") for i in from:to]
#   return lbls
# end

"""
    $(SIGNATURES)

2D plot of all poses, assuming poses are labeled from ``::Symbol` type `:x0, :x1, ..., :xn`.  Use `to` and `from` to limit the range of numbers `n` to be drawn.  The underlying histogram can be enabled or disabled, and the size of maximum-point belief estimate cursors can be controlled with `spscale`.

Future:
- Relax to user defined pose labeling scheme, for example `:p1, :p2, ...`
"""
function drawPoses(fg::G;
                   from::Int64=0,
                   to::Int64=99999999,
                   meanmax=:max,
                   lbls=true,
                   drawhist=true,
                   spscale::Float64=5.0,
                   drawTriads::Bool=true,
                   contour::Bool=true, levels::Int=1,
                   regexPoses=r"x\d",
                   line_width=1pt,
                   manualColor=nothing  ) where G <: AbstractDFG
    #
    @info "drawPoses always sets orientation to max, regardless of meanmax setting.  TODO modefit."

    ## Use PPE.suggested
    vsyms = getVariablesLabelsWithinRange(fg, regexPoses, from=from, to=to)

    Ppes = map(x->getVariablePPE(fg, x), vsyms)
    mask = Ppes .!= nothing
    vsyms = vsyms[mask]
    suggPpes = (x->x.suggested).(Ppes[mask])

    Xpp  = (x->x[1]).(suggPpes)
    Ypp  = (x->x[2]).(suggPpes)
    Thpp = (x->x[3]).(suggPpes)
    LBLS = string.(vsyms)

      ## Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; LBLS=String[];
    # uXpp,uYpp, uThpp, uLBLS = get2DPoseMeans(fg, from=from, to=to)
    # xXpp,xYpp, xThpp, xLBLS = get2DPoseMax(fg, from=from, to=to)
    # # TODO use PPE.suggested instead
    # Xpp  = meanmax == :mean ? uXpp : xXpp
    # Ypp  = meanmax == :mean ? uYpp : xYpp
    # Thpp = xThpp # always use max -- should be modefit
    # LBLS = meanmax == :mean ? uLBLS : xLBLS


    # lbls = lblsFromTo(1,length(Xpp))
    psplt = Union{}
    thm = manualColor == nothing ? Theme(line_width=1pt) : Theme(line_width=line_width, default_color=parse(Colorant, manualColor))
    if lbls
      psplt = Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp,label=LBLS,Geom.path(), thm, Geom.label),
        Coord.cartesian(fixed=true)
      )
    else
      psplt = Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp,Geom.path(), thm),Coord.cartesian(fixed=true),
        Coord.cartesian(fixed=true)
      )
    end
	# return psplt
    drawTriads && addXYLineLayers!(psplt, Xpp, Ypp, Thpp, l=spscale, manualColor=manualColor)
    if drawhist
      Xp,Yp = get2DPoseSamples(fg, from=from, to=to)
      push!(psplt.layers,  Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1] )#(xbincount=100, ybincount=100))
    end
    # add contours to pose estimates
    if contour
      varsyms = Symbol.(LBLS)
      for vsym in varsyms
        pln = plotKDE(fg, vsym, dims=[1;2], levels=levels, c=[(manualColor == nothing ? "gray90" : manualColor);])
        union!(psplt.layers, pln.layers)
      end
    end
    return psplt
end


"""
    $(SIGNATURES)

2D plot of landmarks, assuming `:l1, :l2, ... :ln`.  Use `from` and `to` to control the range of landmarks `n` to include.
"""
function drawLandms(fg::AbstractDFG;
                    from::Int64=0, to::Int64=99999999,
                    minnei::Int64=0,
                    meanmax=:max,
                    lbls=true,showmm=false,drawhist=true,
                    contour::Bool=true, levels::Int=1,
                    manualColor=nothing,
                    c= manualColor==nothing ? "red" : manualColor,
                    MM::Dict{Int,T}=Dict{Int,Int}(),
                    point_size=1pt,
                    regexLandmark::Regex=r"l",
                    resampleGaussianFit::Int=0  ) where T
    #

    ## Use PPE.suggested
    vsyms = getVariablesLabelsWithinRange(fg, regexLandmark, from=from, to=to)

    Ppes = map(x->getVariablePPE(fg, x), vsyms)
    mask = Ppes .!= nothing
    vsyms = vsyms[mask]
    suggPpes = (x->x.suggested).(Ppes[mask])

    Xpp  = (x->x[1]).(suggPpes)
    Ypp  = (x->x[2]).(suggPpes)
    lbltags = string.(vsyms)

    # Xp,Yp = get2DLandmSamples(fg, from=from, to=to)
    # Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; lblstags=String[];
    # # TODO transition to new PPE.suggested
    # if meanmax==:mean
    #   Xpp,Ypp, t, lbltags = get2DLandmMeans(fg, from=from, to=to, regexLandmark=regexLandmark)
    # elseif meanmax==:max
    #   Xpp,Ypp, t, lbltags = get2DLandmMax(fg, from=from, to=to,showmm=showmm,MM=MM, regexLandmark=regexLandmark)
    # end

    if lbls
      psplt = Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp, label=lbltags, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=point_size), Geom.label),
        Coord.cartesian(fixed=true)
        # ,Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)#(xbincount=100, ybincount=100)
      )
    else
      psplt = Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=1pt)),
        Coord.cartesian(fixed=true)
      )
    end

    if drawhist
      push!(psplt.layers, Gadfly.layer(x=Xpp, y=Ypp, Geom.histogram2d)[1])#(xbincount=100, ybincount=100)
    end


    if contour
      # make pretty near Gaussian contours?
      if resampleGaussianFit != 0
        #
      end
      varsyms = Symbol.(lbltags)
      for vsym in varsyms
        @show manualColor
        pln = plotKDE(fg, vsym, dims=[1;2], levels=levels, c=[(manualColor==nothing ? "gray90" : manualColor);])
        union!(psplt.layers, pln.layers)
      end
    end

    psplt
end

"""
    $(SIGNATURES)

2D plot of both poses and landmarks contained in factor graph.  Assuming poses and landmarks are labeled `:x1, :x2, ...` and `:l0, :l1, ...`, respectively.  The rnage of numbers to include can be controlled with `from` and `to` along with other keyword functionality for manipulating the plot.

Notes
- assumes `:l1`, `:l2`, ... for landmarks -- not using `tags=[:LANDMARK]` here yet (TODO).
"""
function drawPosesLandms(fgl::AbstractDFG;
                         from::Int64=0, to::Int64=99999999, minnei::Int64=0,
                         meanmax=:max, lbls=true,
                         drawTriads::Bool=true,
                         spscale::Float64=5.0,
                         contour::Bool=true,
                         levels::Int=1,
                         drawhist=false, MM::Dict{Int,T}=Dict{Int,Int}(),
                         xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing,
                         showmm=true,
                         window::Union{Nothing, Tuple{Symbol, Real}}=nothing,
                         point_size=4pt,
                         line_width=1pt,
                         regexLandmark=r"l",
                         regexPoses=r"x",
                         manualColor=nothing  ) where {T}
  #
  xmin != nothing && xmax != nothing && xmin == xmax ? error("xmin must be less than xmax") : nothing
  ymin != nothing && ymax != nothing && ymin == ymax ? error("ymin must be less than ymax") : nothing
  ll = getVariableIds(fgl, regexLandmark)
  p = drawPoses(fgl, from=from,to=to,meanmax=meanmax,lbls=lbls,drawhist=drawhist, spscale=spscale, contour=contour, drawTriads=drawTriads, manualColor=manualColor, line_width=line_width)
  if length(ll) > 0
    pl = drawLandms(fgl, from=from, to=to, minnei=minnei,lbls=lbls,drawhist=drawhist, MM=MM, showmm=showmm, point_size=point_size, contour=contour, manualColor=manualColor)
    for l in pl.layers
      push!(p.layers, l)
    end
  end
  if window != nothing
    focusX = getKDEMax(getKDE(getVariable(fgl,window[1])))
    pwind = window[2]
    p.coord = Coord.cartesian(xmin=focusX[1]-pwind,xmax=focusX[1]+pwind,ymin=focusX[2]-pwind,ymax=focusX[2]+pwind)
  end
  co = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  p.coord = co
  return p
end

function drawSubmaps(fgl::G, fromto::Array{Int,2};
                     m1hist=false, m2hist=false, m3hist=false,
                     showmm=false, MM::Dict{Int,T} = Dict{Int,Any}(),
                     xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing ) where {G <: AbstractDFG, T}
  #
  p = drawLandms(fgl, from=fromto[1,1], to=fromto[1,2], drawhist=m1hist, showmm=showmm, MM=MM)
  if size(fromto,1) >1
    p2 = drawLandms(fgl, from=fromto[2,1], to=fromto[2,2], drawhist=m2hist,c="blue", showmm=showmm, MM=MM)
    for l in p2.layers
      push!(p.layers, l)
    end
  end
  if size(fromto,1) >2
    p3 = drawLandms(fgl, from=fromto[3,1], to=fromto[3,2], drawhist=m3hist,c="magenta", showmm=showmm, MM=MM)
    for l in p3.layers
      push!(p.layers, l)
    end
  end
  co = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  p.coord = co
  return p
end

function drawSubmaps(fgl::G, fromto::Array{Int,1}; spread::Int=25,
                     m1hist=false, m2hist=false, m3hist=false,
                     showmm=false, MM::Dict{Int,T}=Dict{Int,Any}(),
                     xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing ) where {G <: AbstractDFG, T}
  #
  ft = zeros(Int,length(fromto),2)
  for i in 1:length(fromto)
    ft[i,1] = fromto[i]-spread; ft[i,2] = fromto[i]+spread;
  end
  drawSubmaps(fgl, ft, m1hist=m1hist, m2hist=m2hist, m3hist=m3hist, showmm=showmm, MM=MM, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
end

# function getKDEMax(p::BallTreeDensity;N=200)
#   m = zeros(p.bt.dims)
#   for i in 1:p.bt.dims
#     mm = marginal(p,[i])
#     rangeV = getKDERange(mm)
#     X = linspace(rangeV[1],rangeV[2],N)
#     yV = evaluateDualTree(mm,X)
#     m[i] = X[findfirst(yV,maximum(yV))]
#   end
#   return m
# end


# function plotPose(::Pose2, bels::Vector{BallTreeDensity}, title; levels::Int=5, c=nothing)
#   p1 = plotKDE(bels, dims=[1;2], levels=levels, c=c, title=title)
#   p2 = plotKDE(bels, dims=[3], c=c)
#
#
#   Gadfly.vstack(p1,p2)
# end
import KernelDensityEstimate: getKDERange

function getKDERange(bds::Vector{BallTreeDensity}; extend=0.15)

  dims = Ndim(bds[1])
  ran = getKDERange(bds[1],extend=extend)

  for bd in bds
    rr = getKDERange(bd,extend=extend)
    for i in 2:dims, j in 1:2
      ran[i,j] = maximum([rr[i,j]; ran[i,j]])
    end
  end
  return ran
end

# import RoMEPlotting: plotPose

"""
    $SIGNATURES

Plot pose belief as contour information on visually sensible manifolds.
"""
function plotPose(pt::Pose2,
                  pp::Vector{BallTreeDensity},
                  title="plotPose2";
                  levels=3,
                  c=nothing,
                  axis=nothing,
                  scale::Float64=0.2,
                  overlay=nothing,
                  hdl=[]  )
  #
  # ops = buildHybridManifoldCallbacks(pt.manifolds)
  # @show ran = getKDERange(p, addop=ops[1], diffop=ops[2])
  ran = axis == nothing ? getKDERange(pp) : axis

  p1 = plotKDE(pp, dims=[1;2], levels=levels, c=c, title=title, axis=ran )
  # p2 = plotKDE(bels, dims=[3], c=c)

  cc = c == nothing ? getColorsByLength(length(pp)) : c

  GG = BallTreeDensity[]
  for ppc in pp
    gg = marginal(ppc,[3])
    # gg = (x)->pc(reshape([x], :,1))[1]
    push!(GG, gg)
  end
  # p2 = AMP.plotCircBeliefs(GG, c=cc)
  p2 = AMP.plotKDECircular(GG, scale=scale, c=cc)

  # deal with overlay


  push!(hdl, p1)
  push!(hdl, p2)

  Gadfly.hstack(p1,p2)
end


function plotPose(pt::Pose2,
                  pp::BallTreeDensity,
                  title="plotPose2";
                  levels=3,
                  c=nothing,
                  axis=nothing,
                  scale::Float64=0.2,
                  overlay=nothing,
                  hdl=[]  )
  #
  plotPose(pt, [pp;],title,levels=levels,c=c,axis=axis,scale=scale, overlay=overlay, hdl=hdl )
end


function plotPose(::DynPose2,
                  bels::Vector{BallTreeDensity},
                  title;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  hdl=[],
                  scale::Float64=0.2 )
  #
  p1 = plotKDE(bels, dims=[1;2], levels=levels, c=c, title=title)
  p2 = plotKDE(bels, dims=[3], c=c)
  p3 = plotKDE(bels, dims=[4;5], levels=levels, c=c)

  push!(hdl, p1)
  push!(hdl, p2)
  push!(hdl, p3)

  Gadfly.vstack(p1,p2,p3)
end

# import RoMEPlotting: plotPose

function plotPose(::Pose3,
                  bels::Vector{BallTreeDensity},
                  title;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  hdl=[],
                  scale::Float64=0.2  )
  #
  @show title
  p1 = plotKDE(bels, dims=[1;2], levels=levels, c=c, title=title)
  p2 = plotKDE(bels, dims=[3], c=c)
  p3 = plotKDE(bels, dims=[4;5], levels=levels, c=c)
  p4 = plotKDE(bels, dims=[6], c=c)

  push!(hdl, p1)
  push!(hdl, p2)
  push!(hdl, p3)
  push!(hdl, p4)

  Gadfly.vstack(p1,p2,p3,p4)
end

"""
    $(SIGNATURES)

Example: pl = plotPose(fg, [:x1; :x2; :x3])
"""
function plotPose(fgl::G,
                  syms::Vector{Symbol};
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  scale::Float64=0.2,
                  show::Bool=false,
                  filepath::AS="/tmp/tempposeplot.svg",
                  app::AS="eog",
                  hdl=[]  ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  typ = getSolverData(getVariable(fgl, syms[1])).softtype
  pt = string(string.(syms)...)
  getvertsgg = (sym) -> getKDE(getVariable(fgl, sym))
  pl = plotPose(typ, getvertsgg.(syms), pt, levels=levels, c=c, axis=axis, scale=scale, hdl=hdl )

  if length(filepath) > 0
    ext = split(filepath, '.')[end]
    cmd = getfield(Gadfly,Symbol(uppercase(ext)))

    h = 1*7Gadfly.cm
    if typ == DynPose2
        h *= 1.5
    end
    Gadfly.draw(cmd(filepath,15Gadfly.cm,h),pl)


    @async !show ? nothing : run(`$app $filepath`)
  end
  return pl
end

function plotPose(fgl::G,
                  sym::Symbol;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  scale::Float64=0.2,
                  show::Bool=false,
                  filepath::AS="/tmp/tempposeplot.svg",
                  app::AS="eog",
                  hdl=[]  ) where {G <: AbstractDFG, AS <: AbstractString}
  #
  plotPose(fgl, [sym;], levels=levels, axis=axis, show=show, filepath=filepath, app=app, hdl=hdl )
end

# deprecated
function investigatePoseKDE(p::BallTreeDensity, p0::BallTreeDensity)
    # co = ["black"; "blue"]
    # h = Union{}
    # x = plotKDE([marginal(p,[1]); marginal(p0,[1])], c=co )
    # y = plotKDE([marginal(p,[2]); marginal(p0,[2])], c=co )
    # if p.bt.dims >= 3
    #   th = plotKDE([marginal(p,[3]); marginal(p0,[3])], c=co )
    #   h = hstack(x,y,th)
    # else
    #   h = hstack(x,y)
    # end
    #
    # return h
    return investigateMultidimKDE(p, p0)
end


function investigatePoseKDE(p::Array{BallTreeDensity,1})
    # co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2";
    # "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta";
    # "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"]
    # # compute all the marginals
    # Pm = Array{Array{BallTreeDensity,1},1}()
    # push!(Pm,stackMarginals(p,1)) #[marginal(p[1],[1]); marginal(p[2],[1])]
    # push!(Pm,stackMarginals(p,2)) #[marginal(p[1],[2]); marginal(p[2],[2])]
    #
    # h = Union{}
    # x = plotKDE(Pm[1], c=co )
    # y = plotKDE(Pm[2], c=co )
    # if p[1].bt.dims >= 3
    #   #Pm3 = [marginal(p[1],[3]); marginal(p[2],[3])]
    #   push!(Pm,stackMarginals(p,3)) # [marginal(p[1],[3]); marginal(p[2],[3])]
    #   th = plotKDE(Pm[3], c=co )
    #   h = hstack(x,y,th)
    # else
    #   h = hstack(x,y)
    # end
    # return h
    return investigateMultidimKDE(p)
end

function investigatePoseKDE(p::BallTreeDensity)
    # x = plotKDE(marginal(p,[1]) )
    # y = plotKDE(marginal(p,[2]) )
    # if p.bt.dims >= 3
    #   th = plotKDE(marginal(p,[3]) )
    #   return hstack(x,y,th)
    # end
    # return hstack(x,y)
    return investigateMultidimKDE(p)
end

# import RoMEPlotting: drawMarginalContour

function drawMarginalContour(fgl::G, lbl::String;
    xmin=-150,xmax=150,ymin=-150,ymax=150,n=200 ) where G <: AbstractDFG
  #
  p = getKDE(getVariable(fgl,Symbol(lbl)))  # p = getKDE(getVert(fgl,lbl))
  Gadfly.plot(z=(x,y)->evaluateDualTree(p,vectoarr2([x,y]))[1],
    x=collect(range(xmin,stop=xmax,length=n)),
    y=collect(range(ymin,stop=ymax,length=n)),
    Geom.contour,
    Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    Guide.title(lbl)
  )
end

function accumulateMarginalContours(fgl, order;
    xmin=-150,xmax=150,ymin=-150,ymax=150,n=200 )
  #
  pl = drawMarginalContour(fgl, order[1],xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n)
  pl2 = nothing
  PL = []
  for or in order[1:end]
    pl2 = drawMarginalContour(fgl, or, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n)
    push!(PL, pl2)
    push!(pl.layers, pl2.layers[1])
  end
  return pl, PL
end




function plotPose3Pairs(fgl::G, sym::Symbol; fill::Bool=true) where {G <: AbstractDFG}
  p1= plotKDE(fgl, :x1, dims=[1;2], fill=fill)
  p2 = plotKDE(fgl, :x1, dims=[6;3], fill=fill)
  p3 = plotKDE(fgl, :x1, dims=[4;5], fill=fill)
  Gadfly.draw(PDF("/tmp/RoMEvstackPose3.pdf",15cm, 20cm), vstack(p1,p2,p3) )
  @async run(`evince /tmp/RoMEvstackPose3.pdf`)
  nothing
end





# Victoria Park Plotting functions


function progressExamplePlot(dOdo, lsrFeats; toT=Inf)
    len = length(dOdo)
    pose = SE2(zeros(3))
    lastpose = zeros(3)
    idx = 1
    T = dOdo[idx][4]
    lstlaseridx = 1
    WFTSX = Array{Float64,1}()
    WFTSY = Array{Float64,1}()
    WLBLS = ASCIIString[]

    lastX = Array{Float64,1}()
    lastY = Array{Float64,1}()

    while T < toT && idx <= len

      lastX = Array{Float64,1}()
      lastY = Array{Float64,1}()
      pose = pose*SE2(dOdo[idx][1:3]) # todo -- replace with inferred latest pose
      #@show idx, T, vec(pose[1:2,3])
      lastpose = vec(se2vee(pose))

      # lstlaseridx, Ta = getFeatsAtT(lsrFeats, T, prev=lstlaseridx)
      # bfts = lsrFeats[lstlaseridx].feats
      fe = lsrFeats[idx]
      if length(lsrFeats[idx]) > 0
        bfts = zeros(3,length(fe))
        lbls = ASCIIString[]
        k = collect(keys(fe))
        for i in 1:length(fe)
          bfts[1:length(fe[k[i]]),i] = fe[k[i]]
          push!(lbls, "l$(k[i])")
        end


        if bfts[1,1] != 0.0 && bfts[2,1] != 0.0 && bfts[3,1] != 0.0
          wfts = rotateFeatsToWorld(bfts, pose)
          for i in 1:size(wfts,2)
              push!(WFTSX, wfts[1,i])
              push!(WFTSY, wfts[2,i])
              push!(WLBLS, lbls[i])
              push!(lastX, wfts[1,i])
              push!(lastY, wfts[2,i])
          end
        end
      end
      idx += 1
      if idx <= len
        T = dOdo[idx][4]
      end
    end

    p = plotPoseDict(dOdo,to=idx-1)
    if length(WFTSX) > 0
      l = Gadfly.layer(x=WFTSX, y=WFTSY, label=WLBLS, Geom.label, Geom.point, Gadfly.Theme(default_color=colorant"red"))
      push!(p.layers, l[1])
      l2 = Gadfly.layer(x=WFTSX, y=WFTSY, Geom.point, Gadfly.Theme(default_color=colorant"red"))
      push!(p.layers, l2[1])
      for i in 1:length(lastX)
        push!(p.layers, Gadfly.layer(x=[lastpose[1];lastX[i]], y=[lastpose[2];lastY[i]], Geom.line, Gadfly.Theme(default_color=colorant"magenta"))[1])
      end
    end
    p
end


function plotTrckStep(DBG, i, fid, m)
  @show keys(DBG[i])
  pf = DBG[i][fid]
  arr = Array{BallTreeDensity,1}()
  for j in 1:3
    push!(arr, marginal(pf[j],[m]))
  end
  plotKDE(arr, c=["red";"green";"black"])
end



function plotPose3Pairs(fgl::FactorGraph, sym::Symbol; fill::Bool=true)
  p1= plotKDE(fgl, sym, dims=[1;2], fill=fill)
  p2 = plotKDE(fgl, sym, dims=[6;3], fill=fill)
  p3 = plotKDE(fgl, sym, dims=[4;5], fill=fill)
  Gadfly.draw(PDF("/tmp/RoMEvstackPose3.pdf",15cm, 20cm), vstack(p1,p2,p3) )
  @async run(`evince /tmp/RoMEvstackPose3.pdf`)
  nothing
end


function plotKDE(fgl::FactorGraph,
                 vsym::Vector{Symbol};
                 axis=nothing,
                 dims=nothing,
                 c=getColorsByLength(length(vsym)),
                 levels::Int=4,
                 title::Union{Nothing, T}=nothing,
                 overlay=nothing  ) where {T <: AbstractString}
  #
  @warn "plotKDE for FactorGraph is deprecated, use DistributedFactorGraphs objects instead."
  verts = map((x)->getKDE(getVariable(fgl, x)), vsym)
  plotKDE(verts, dims=dims, c=c, axis=axis, levels=levels, title=title, overlay=overlay )
end

function plotKDE(fgl::FactorGraph, vsym::Symbol; axis=nothing, dims=nothing, c=nothing, levels=4, title::Union{Nothing, T}=nothing) where {T <: AbstractString}
  @warn "plotKDE for FactorGraph is deprecated, use DistributedFactorGraphs objects instead."
  plotKDE(fgl, Symbol[vsym;], dims=dims, c=c, axis=axis, levels=levels, title=title)
end


"""
    $SIGNATURES

Convenience function to plot one Point2 or Pose2 location along with reference data if desired.
"""
function plotVariable2D(dfg::AbstractDFG,
                        varsym::Symbol;
                        refs::Vector=[],
                        levels::Int=10 )
  #
  # make sure variable is in the right family
  var = getVariable(dfg, varsym)
  @assert isa(getSofttype(var), Union{Pose2, Point2})
  pl = plotKDE(dfg, varsym, levels=levels)
  if 0 < length(refs)
    XX, YY = zeros(0), zeros(0)
    for dict in refs
        push!(XX,dict[varsym][1])
        push!(YY,dict[varsym][2])
    end
    p2 = Gadfly.plot(x=XX,
                     y=YY,
                     Geom.point,
                     Guide.Theme(default_color=colorant"red", point_size=5pt))
    union!(p2.layers, pl.layers)
    pl = p2
  end
  return pl
end
# pl = plotKDE(dfg, varsym, levels=levels)
# pl = Gadfly.plot(x=[landmarks_design[:l1][1]; landmarks_real[:l1][1]],
# y=[landmarks_design[:l1][2]; landmarks_real[:l1][2]],
# Geom.point,
# Guide.Theme(default_color=colorant"red", point_size=5pt))
# p2 = plotKDE(fg, :l1, levels=20)
# union!(pl.layers, p2.layers)



function plotTrailingPoses(pt::Pose2,
                           pp::Vector{BallTreeDensity},
                           title="";
                           levels=2,
                           c=nothing,
                           axis=nothing,
                           scale::Float64=0.2,
                           circlen::Int=5)

ran = axis == nothing ? getKDERange(pp) : axis

cc=["red"; ["pink" for i in 1:100]]

p1 = plotKDE(pp, dims=[1;2], levels=levels, c=cc, title=title, axis=ran )

GG = BallTreeDensity[]
for ppc in pp
  gg = marginal(ppc,[3])
  # gg = (x)->pc(reshape([x], :,1))[1]
  push!(GG, gg)
end
p2 = AMP.plotKDECircular(GG[(end-circlen):end], scale=scale, c=cc)

p2,p1
end



function plotTrailingPoses(fg::G,
                           pp::Vector{Symbol},
                           title="";
                           levels=2,
                           c=nothing,
                           axis=nothing,
                           scale::Float64=0.2,
                           circlen::Int=5) where G <: AbstractDFG
  #
  plotTrailingPoses(Pose2(), map(x->getKDE(fg,x),pp), scale=scale, title=title, circlen=circlen)
end

# gg = (x)->plotTrailingPoses(fg, [Symbol("x$i") for i in (x+60):-5:x],circlen=3)
#
# for i in 5:5:290
#  g1,g2 = gg(i)
#
#  g1 |> SVG("/tmp/trailingimgs/g1_$(i).svg")
#  g1 |> SVG("/tmp/trailingimgs/g1_$(i+1).svg")
#  g1 |> SVG("/tmp/trailingimgs/g1_$(i+2).svg")
#  g1 |> SVG("/tmp/trailingimgs/g1_$(i+3).svg")
#  g1 |> SVG("/tmp/trailingimgs/g1_$(i+4).svg")
#
#  g2 |> SVG("/tmp/trailingimgs/g2_$(i).svg")
#  g2 |> SVG("/tmp/trailingimgs/g2_$(i+1).svg")
#  g2 |> SVG("/tmp/trailingimgs/g2_$(i+2).svg")
#  g2 |> SVG("/tmp/trailingimgs/g2_$(i+3).svg")
#  g2 |> SVG("/tmp/trailingimgs/g2_$(i+4).svg")
# end

#
