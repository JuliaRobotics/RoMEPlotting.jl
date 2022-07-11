# new exports from this file, WIP

export calcDyadScaleAdaptive
export covEllipseParameterized, plotCovEllipseLayer
export plotSLAM2D_KeyAndRef

## Src



# draw the reference frame as a red-green dyad
function addXYLineLayers!(pl, Xpp, Ypp, Thpp; l::Real=1.0, manualColor::Union{Nothing, AbstractString}=nothing  )
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

      push!(pl.layers, layer(x=xsr, y=ysr, Geom.path(), Gadfly.Theme(default_color=parse(Colorant, manualColor === nothing ? "red" : manualColor), line_width=1.5pt))[1] )
      push!(pl.layers, layer(x=xsg, y=ysg, Geom.path(), Gadfly.Theme(default_color=parse(Colorant, manualColor === nothing ? "green" : manualColor), line_width=1.5pt))[1] )
    end
    nothing
end



"""
    $SIGNATURES
Adaptively size the `dyadScale` value for plotSLAM2DPose.
"""
function calcDyadScaleAdaptive( Xpp::AbstractVector{<:Real},
                                Ypp::AbstractVector{<:Real};
                                scaling::Real=0.25)

  dists = diff(Xpp).^2 + diff(Ypp).^2 |> vec .|> sqrt
  scaling*Statistics.mean(dists)
end



"""
    $(SIGNATURES)

2D plot of all poses, assuming poses are labeled from ``::Symbol` type `:x0, :x1, ..., :xn`.  Use `to` and `from` to limit the range of numbers `n` to be drawn.  The underlying histogram can be enabled or disabled, and the size of maximum-point belief estimate cursors can be controlled with `spscale`.

Future:
- Relax to user defined pose labeling scheme, for example `:p1, :p2, ...`
"""
function plotSLAM2DPoses( fg::AbstractDFG;
                          solveKey::Symbol=:default,
                          regexPoses=r"x\d",
                          from::Int=0,
                          to::Int=(2^(Sys.WORD_SIZE-1)-1),
                          variableList::AbstractVector{Symbol}=listVariablesLabelsWithinRange(fg, regexPoses, from=from, to=to),
                          meanmax=:null,
                          ppe=:suggested,
                          recalcPPEs::Bool=false,
                          lbls=true,
                          scale::Real=1,
                          x_off::Real=0,
                          y_off::Real=0,
                          drawhist=false,
                          spscale::Union{Nothing, <:Real}=nothing,
                          dyadScale::Union{Nothing, <:Real}=nothing,
                          drawTriads::Bool=true,
                          drawContour::Bool=true, levels::Int=1,
                          contour::Union{Nothing, Bool}=nothing,
                          line_width=1pt,
                          drawPoints::Bool=true,
                          pointsColor::AbstractString="gray30",
                          drawEllipse::Bool=false,
                          ellipseColor::AbstractString="gray30",
                          manualColor=nothing  )
    #
    # deprecations
    if meanmax != :null
      @warn "plotSLAM2DPoses meanmax keyword is deprecated, use ppe=:suggested instead."
      ppe = meanmax
    end
    !(spscale isa Nothing) ? (dyadScale=spscale; @warn("keyword spscale is deprecated, use dyadScale instead")) : nothing
    !(contour isa Nothing) ? (drawContour=contour; @warn("keyword contour is being deprecated, use drawContour instead")) : nothing

    ## Use PPE.suggested

    Ppes = if recalcPPEs
      map(x->calcVariablePPE(fg, x, solveKey=solveKey), variableList)
    else
      getPPE.(fg, variableList, solveKey)
    end
    mask = Ppes .!== nothing
    variableList = variableList[mask]
    suggPpes = (x->getfield(x,ppe)).(Ppes[mask])

    Xpp::Vector{Float64}  = Float64.( (x->x[1]).(suggPpes) )
    Ypp::Vector{Float64}  = Float64.( (x->x[2]).(suggPpes) )
    Thpp::Vector{Float64} = Float64.( (x->x[3]).(suggPpes) )
    LBLS = string.(variableList)

    Xpp .+= x_off
    Ypp .+= y_off

    Xpp .*= scale
    Ypp .*= scale

    # adaptively scale dyad size
    dyadScale = dyadScale isa Nothing ? calcDyadScaleAdaptive(Xpp, Ypp) : dyadScale

    # lbls = lblsFromTo(1,length(Xpp))
    psplt = Union{}
    thm = manualColor === nothing ? Theme(line_width=1pt) : Theme(line_width=line_width, default_color=parse(Colorant, manualColor))
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
    drawTriads && addXYLineLayers!(psplt, Xpp, Ypp, Thpp, l=dyadScale, manualColor=manualColor)
    if drawhist
      Xp,Yp = get2DPoseSamples(fg, from=from, to=to)
      push!(psplt.layers,  Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)[1] )#(xbincount=100, ybincount=100))
    end
    # add contours to pose estimates
      # varsyms = Symbol.(LBLS)
    if drawContour
      for vsym in variableList
        pln = plotKDE(fg, vsym, solveKey=solveKey, dims=[1;2], levels=levels, c=[(manualColor === nothing ? "gray90" : manualColor);])
        union!(psplt.layers, pln.layers)
      end
    end
    # drawEllipse
    for vsym in variableList
      pln = plotCovEllipseLayer(fg, vsym, solveKey=solveKey, drawEllipse=drawEllipse,drawPoints=drawPoints,ellipseColor=ellipseColor,pointsColor=pointsColor)
      union!(psplt.layers, pln)
    end
    return psplt
end



"""
    $(SIGNATURES)

2D plot of landmarks, assuming `:l1, :l2, ... :ln`.  Use `from` and `to` to control the range of landmarks `n` to include.
"""
function plotSLAM2DLandmarks( fg::AbstractDFG;
                              solveKey::Symbol=:default,
                              regexLandmark::Regex=r"l",
                              from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1),
                              minnei::Int=0,
                              variableList::AbstractVector{Symbol}=listVariablesLabelsWithinRange(fg, regexLandmark, from=from, to=to),
                              meanmax=:null,
                              ppe::Symbol=:suggested,
                              recalcPPEs::Bool=false,
                              lbls=true,showmm=false,
                              scale::Real=1,
                              x_off::Real=0,
                              y_off::Real=0,
                              drawhist=false,
                              drawContour::Bool=true, levels::Int=1,
                              contour::Union{Nothing, Bool}=nothing,
                              manualColor=nothing,
                              c= manualColor===nothing ? "red" : manualColor,
                              MM::Dict{Int,T}=Dict{Int,Int}(),
                              point_size=1pt,
                              drawPoints::Bool=true,
                              pointsColor::AbstractString="gray30",
                              drawEllipse::Bool=false,
                              ellipseColor::AbstractString="gray30",
                              resampleGaussianFit::Int=0  ) where T
    #
    if meanmax != :null
      @error "plotSLAM2DPoses meanmax keyword is deprecated, use ppe instead."
      ppe = meanmax
    end

    # remove before v0.10
    !(contour isa Nothing) ? (drawContour=contour; @warn("keyword contour is being deprecated, use drawContour instead")) : nothing

    ## Use PPE.suggested

    Ppes = if recalcPPEs
      map(x->calcVariablePPE(fg, x, solveKey=solveKey), variableList)
    else
      getPPE.(fg, variableList, solveKey)
    end
    mask = Ppes .!== nothing
    variableList = variableList[mask]
    suggPpes = (x->getfield(x,ppe)).(Ppes[mask])

    Xpp  = (x->x[1]).(suggPpes)
    Ypp  = (x->x[2]).(suggPpes)
    lbltags = string.(variableList)

    Xpp .+= x_off
    Ypp .+= y_off

    Xpp .*= scale
    Ypp .*= scale

    # Xp,Yp = get2DLandmSamples(fg, from=from, to=to)
    # Xpp = Float64[]; Ypp=Float64[]; Thpp=Float64[]; lblstags=String[];
    # # TODO transition to new PPE.suggested
    # if meanmax==:mean
    #   Xpp,Ypp, t, lbltags = get2DLandmMeans(fg, from=from, to=to, regexLandmark=regexLandmark)
    # elseif meanmax==:max
    #   Xpp,Ypp, t, lbltags = get2DLandmMax(fg, from=from, to=to,showmm=showmm,MM=MM, regexLandmark=regexLandmark)
    # end

    psplt = if lbls
      Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp, label=lbltags, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=point_size), Geom.label),
        Coord.cartesian(fixed=true)
        # ,Gadfly.layer(x=Xp, y=Yp, Geom.histogram2d)#(xbincount=100, ybincount=100)
      )
    else
      Gadfly.plot(
        Gadfly.layer(x=Xpp,y=Ypp, Geom.point, Theme(line_width=1pt, default_color=parse(Colorant,c), point_size=point_size)),
        Coord.cartesian(fixed=true)
      )
    end

    if drawhist
      push!(psplt.layers, Gadfly.layer(x=Xpp, y=Ypp, Geom.histogram2d)[1]) #(xbincount=100, ybincount=100)
    end


    if drawContour
      # make pretty near Gaussian contours?
      if resampleGaussianFit != 0
        #
      end
      varsyms = Symbol.(lbltags)
      for vsym in varsyms
        pln = plotKDE(fg, vsym, solveKey=solveKey, dims=[1;2], levels=levels, c=[(manualColor===nothing ? "gray90" : manualColor);])
        union!(psplt.layers, pln.layers)
      end
    end

    # if drawEllipse
    for vsym in variableList
      pln = plotCovEllipseLayer(fg, vsym, solveKey=solveKey, drawEllipse=drawEllipse,drawPoints=drawPoints,ellipseColor=ellipseColor,pointsColor=pointsColor)
      union!(psplt.layers, pln)
    end

    psplt
end


"""
    $(SIGNATURES)

2D plot of both poses and landmarks contained in factor graph.  Assuming poses and landmarks 
are labeled `:x1, :x2, ...` and `:l0, :l1, ...`, respectively.  The range of numbers to 
include can be controlled with `from` and `to` along with other keyword functionality for 
manipulating the plot.

Notes
- Assumes `:l1`, `:l2`, ... for landmarks -- 
- Can increase default Gadfly plot size (for JSSVG in browser): `Gadfly.set_default_plot_size(35cm,20cm)`.
- Enable or disable features such as the covariance ellipse with keyword `drawEllipse=true`.

DevNotes
- TODO update to use e.g. `tags=[:LANDMARK]`,
- TODO fix `drawHist`,
- TODO deprecate, `showmm`, `spscale`.

Examples:
```julia
fg = generateGraph_Hexagonal()
plotSLAM2D(fg)
plotSLAM2D(fg, drawPoints=false)
plotSLAM2D(fg, contour=false, drawEllipse=true)
plotSLAM2D(fg, contour=false, title="SLAM result 1")

# or load a factor graph
fg_ = loadDFG("somewhere.tar.gz")
plotSLAM2D(fg_)
```

Related

[`plotSLAM2DPoses`](@ref), [`plotSLAM2DLandmarks`](@ref), [`plotPose`](@ref), [`plotBelief`](@ref) 
"""
function plotSLAM2D(fgl::AbstractDFG;
                    solveKey::Symbol=:default,
                    from::Int=0, to::Int=(2^(Sys.WORD_SIZE-1)-1), 
                    minnei::Int=0,
                    meanmax=:null,
                    posesPPE=:suggested,
                    landmsPPE=:suggested,
                    recalcPPEs::Bool=false,
                    lbls=true,
                    scale::Real=1,
                    x_off::Real=0,
                    y_off::Real=0,
                    drawTriads::Bool=true,
                    spscale::Union{Nothing, <:Real}=nothing,
                    dyadScale::Union{Nothing, <:Real}=nothing,
                    levels::Int=1,
                    drawhist=false, MM::Dict{Int,T}=Dict{Int,Int}(),
                    xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing,
                    showmm=true,
                    window::Union{Nothing, Tuple{Symbol, Real}}=nothing,
                    point_size=4pt,
                    line_width=1pt,
                    regexLandmark=r"l\d+",
                    regexPoses=r"x\d+",
                    variableList::AbstractVector{Symbol}=union(listVariablesLabelsWithinRange(fgl, regexPoses, from=from, to=to), listVariablesLabelsWithinRange(fgl, regexLandmark, from=from, to=to)),
                    manualColor=nothing,
                    drawPoints::Bool=true,
                    pointsColor::AbstractString="gray30",
                    drawContour::Bool=true,
                    drawEllipse::Bool=false,
                    ellipseColor::AbstractString="gray30",
                    contour::Union{Nothing,Bool}=nothing,
                    title::AbstractString="",
                    aspect_ratio::Real=1  ) where T
  #
  # deprecations
  if meanmax != :null
    @warn "plotSLAM2DPosesLandms meanmax keyword is deprecated, use posesPPE or landmsPPE instead."
    posesPPE = meanmax
    landmsPPE = meanmax
  end
  !(spscale isa Nothing) ? (dyadScale=spscale; @error("keyword spscale is being deprecated, use dyadScale instead")) : nothing
  !(contour isa Nothing) ? (drawContour=contour; @error("keyword contour is being deprecated, use drawContour instead")) : nothing

  #
  xmin !== nothing && xmax !== nothing && xmin == xmax ? error("xmin must be less than xmax") : nothing
  ymin !== nothing && ymax !== nothing && ymin == ymax ? error("ymin must be less than ymax") : nothing
  
  p = plotSLAM2DPoses(fgl;
                      solveKey,
                      from,
                      to,
                      regexPoses,
                      variableList=intersect(variableList, listVariablesLabelsWithinRange(fgl, regexPoses; from, to)),
                      ppe=posesPPE,
                      lbls,
                      scale,
                      x_off,
                      y_off,
                      drawhist,
                      dyadScale,
                      drawContour,
                      drawTriads,
                      manualColor,
                      line_width,
                      drawPoints,
                      ellipseColor,
                      pointsColor,
                      drawEllipse,
                      recalcPPEs  )
  #

  ll = listVariables(fgl, regexLandmark)
  if length(ll) > 0
    pl = plotSLAM2DLandmarks( fgl;
                              solveKey,
                              from,
                              to,
                              regexLandmark,
                              variableList=intersect(variableList, listVariablesLabelsWithinRange(fgl, regexLandmark; from, to)),
                              ppe=landmsPPE,
                              minnei,
                              lbls,
                              scale,
                              x_off,
                              y_off,
                              drawhist,
                              MM,
                              showmm,
                              point_size,
                              drawContour,
                              manualColor,
                              drawPoints,
                              ellipseColor,
                              pointsColor,
                              drawEllipse,
                              recalcPPEs  )
    #

    for l in pl.layers
      push!(p.layers, l)
    end
  end
  if window !== nothing
    focusX = getPPE(fgl, window[1], solveKey).max # getKDEMax( getBelief(getVariable(fgl,window[1]),solveKey) )
    pwind = window[2]
    p.coord = Coord.cartesian(xmin=focusX[1]-pwind,xmax=focusX[1]+pwind,ymin=focusX[2]-pwind,ymax=focusX[2]+pwind)
  end
  co = Coord.Cartesian(;xmin,xmax,ymin,ymax,aspect_ratio)
  p.coord = co
  if title != ""
    push!(p.guides, Guide.title(title))
  end
  return p
end


function plotSLAM2DSubmaps( fgl::AbstractDFG, fromto::Array{Int,2};
                            m1hist=false, m2hist=false, m3hist=false,
                            showmm=false, MM::Dict{Int,T} = Dict{Int,Any}(),
                            xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing ) where T
  #
  p = plotSLAM2DLandmarks(fgl, from=fromto[1,1], to=fromto[1,2], drawhist=m1hist, showmm=showmm, MM=MM)
  if size(fromto,1) >1
    p2 = plotSLAM2DLandmarks(fgl, from=fromto[2,1], to=fromto[2,2], drawhist=m2hist,c="blue", showmm=showmm, MM=MM)
    for l in p2.layers
      push!(p.layers, l)
    end
  end
  if size(fromto,1) >2
    p3 = plotSLAM2DLandmarks(fgl, from=fromto[3,1], to=fromto[3,2], drawhist=m3hist,c="magenta", showmm=showmm, MM=MM)
    for l in p3.layers
      push!(p.layers, l)
    end
  end
  co = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  p.coord = co
  return p
end

function plotSLAM2DSubmaps( fgl::G, fromto::Array{Int,1}; spread::Int=25,
                            m1hist=false, m2hist=false, m3hist=false,
                            showmm=false, MM::Dict{Int,T}=Dict{Int,Any}(),
                            xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing ) where {G <: AbstractDFG, T}
  #
  ft = zeros(Int,length(fromto),2)
  for i in 1:length(fromto)
    ft[i,1] = fromto[i]-spread; ft[i,2] = fromto[i]+spread;
  end
  plotSLAM2DSubmaps(fgl, ft, m1hist=m1hist, m2hist=m2hist, m3hist=m3hist, showmm=showmm, MM=MM, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
end

function  plotSLAM2D_KeyAndRef(fg, solveKey=:default; refkey=:simulated, recalcPPEs=true)
  #
  pl_ = plotSLAM2D(fg, solveKey=refkey, drawContour=false, drawPoints=false, drawEllipse=false, manualColor="black", drawTriads=false);
  pl1 = plotSLAM2D(fg, solveKey=solveKey, recalcPPEs=recalcPPEs, drawContour=false, drawPoints=false, drawEllipse=false);

  union!(pl1.layers, pl_.layers);

  pl1
end



# import RoMEPlotting: drawMarginalContour

function plotMarginalContour( fgl::AbstractDFG, lbl::String;
                              solveKey::Symbol=:default,
                              xmin=-150,xmax=150,ymin=-150,ymax=150,
                              n::Int=200 )
  #
  p = getBelief(getVariable(fgl,Symbol(lbl)), solveKey)  # p = getBelief(getVert(fgl,lbl))
  Gadfly.plot(z=(x,y)->evaluateDualTree(p,vectoarr2([x,y]))[1],
    x=collect(range(xmin,stop=xmax,length=n)),
    y=collect(range(ymin,stop=ymax,length=n)),
    Geom.contour,
    Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    Guide.title(lbl)
  )
end
# ┌ Warning: `Scale.color_none` to be deprecated. Instead use e.g. `plot(..., Geom.contour, color=[colorant"black"])`
# └ @ Gadfly.Scale ~/.julia/packages/Gadfly/B5yQc/src/scale.jl:446


function accumulateMarginalContours(fgl, order;
                                    solveKey::Symbol=:default,
                                    xmin=-150,xmax=150,ymin=-150,ymax=150,n=200 )
  #
  pl = plotMarginalContour(fgl, order[1],solveKey=solveKey, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n)
  pl2 = nothing
  PL = []
  for or in order[1:end]
    pl2 = plotMarginalContour(fgl, or, solveKey=solveKey, xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,n=n)
    push!(PL, pl2)
    push!(pl.layers, pl2.layers[1])
  end
  return pl, PL
end



#
