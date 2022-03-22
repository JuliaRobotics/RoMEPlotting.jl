# Currently unmaintained



function plotLsrScanFeats(br::Array{Float64,2})
  Cart = zeros(size(br))
  Cart[:,1] = br[:,2].*cos(br[:,1])
  Cart[:,2] = br[:,2].*sin(br[:,1])
  plot(x=Cart[:,1],y=Cart[:,2],Geom.point,
  Guide.xticks(ticks=collect(-60:10:60)),
  Guide.yticks(ticks=collect(0:10:80)))
end

function plotFeatTrackers(trkrs::Dict{Int64,Feature}, bfts::Array{Float64,2})
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

Plot trajectory of Array{,2} with rows as consecutive entries and columns as x,y,theta.
"""
function plotTrajectoryArrayPose2(arr::AbstractMatrix__{<:Real};
                                  spscale::Real=0.5,
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
    WLBLS = String[]

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
        lbls = String[]
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



function plotPose3Pairs(fgl::AbstractDFG, sym::Symbol; fill::Bool=true)
  p1= plotKDE(fgl, sym, dims=[1;2], fill=fill)
  p2 = plotKDE(fgl, sym, dims=[6;3], fill=fill)
  p3 = plotKDE(fgl, sym, dims=[4;5], fill=fill)
  Gadfly.draw(PDF("/tmp/RoMEvstackPose3.pdf",15cm, 20cm), vstack(p1,p2,p3) )
  @async run(`evince /tmp/RoMEvstackPose3.pdf`)
  nothing
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
                           scale::Real=0.2,
                           circlen::Int=5)

ran = axis === nothing ? getKDERange(pp) : axis

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



function plotTrailingPoses(fg::AbstractDFG,
                           pp::Vector{Symbol},
                           title="";
                           solveKey::Symbol=:default,
                           levels=2,
                           c=nothing,
                           axis=nothing,
                           scale::Real=0.2,
                           circlen::Int=5)
  #
  plotTrailingPoses(Pose2(), map(x->getBelief(fg,x, solveKey),pp), scale=scale, title=title, circlen=circlen)
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





function investigateMultidimKDE(p::BallTreeDensity, p0::BallTreeDensity)
  co = ["black"; "blue"]
  h = Union{}
  x = plotKDE([marginal(p,[1]); marginal(p0,[1])], c=co )
  y = plotKDE([marginal(p,[2]); marginal(p0,[2])], c=co )
  if p.bt.dims >= 3
    th = plotKDE([marginal(p,[3]); marginal(p0,[3])], c=co )
    h = hstack(x,y,th)
  else
    h = hstack(x,y)
  end

  return h
end


function investigateMultidimKDE(p::Array{BallTreeDensity,1})
  co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2";
  "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"; "magenta";
  "cyan"; "cyan1"; "cyan2"; "magenta"; "cyan"; "cyan1"; "cyan2"]
  # compute all the marginals
  Pm = Array{Array{BallTreeDensity,1},1}()
  push!(Pm,stackMarginals(p,1)) #[marginal(p[1],[1]); marginal(p[2],[1])]
  push!(Pm,stackMarginals(p,2)) #[marginal(p[1],[2]); marginal(p[2],[2])]

  h = Union{}
  x = plotKDE(Pm[1], c=co )
  y = plotKDE(Pm[2], c=co )
  if p[1].bt.dims >= 3
    #Pm3 = [marginal(p[1],[3]); marginal(p[2],[3])]
    push!(Pm,stackMarginals(p,3)) # [marginal(p[1],[3]); marginal(p[2],[3])]
    th = plotKDE(Pm[3], c=co )
    h = hstack(x,y,th)
  else
    h = hstack(x,y)
  end
  return h
end

function investigateMultidimKDE(p::BallTreeDensity)
  x = plotKDE(marginal(p,[1]) )
  y = plotKDE(marginal(p,[2]) )
  if p.bt.dims >= 3
    th = plotKDE(marginal(p,[3]) )
    return hstack(x,y,th)
  end
  return hstack(x,y)
end


function plotLbl(fgl::G, lbl::Symbol) where G <: AbstractDFG
  plotKDE(getBelief(getVariable(fgl,lbl)))
end
drawLbl(fgl::G, lbl::T) where {G <: AbstractDFG, T <: AbstractString} = drawLbl(fgl, Symbol(lbl))




function plotFrontalDens( fg::AbstractDFG,
                          bt::AbstractBayesTree;
                          N=300,
                          gt=Union{} )
    #
    len = length(bt.cliques)
    vv = Array{Gadfly.Compose.Context,1}(len)
    i = 0
    for cliq in bt.cliques
        #@show cliq[2].attributes["label"]
        lenfr = length(cliq[2].attributes["data"].frontalIDs)

        p = Array{BallTreeDensity,1}(lenfr)
        j=0
        #pvals = Array{Array{Float64,2},1}(lenfr)
        gtvals = Dict{Int,Array{Float64,2}}()
        lbls = String[]

        for frid in cliq[2].attributes["data"].frontalIDs
            j+=1
            p[j] = getBelief(getVariable(fg, frid)) # getBelief(fg.v[frid])
            # p[j] = kde!(fg.v[frid].attributes["val"])

            #pvals[j] = fg.v[frid].attributes["val"]

            if gt!=Union{}
              gtvals[j] = gt[getVariable(fg,frid).label] # fg.v[frid].
              #push!(gtvals, gt[fg.v[frid].attributes["label"]][1])
              #push!(gtvals, gt[fg.v[frid].attributes["label"]][2])
            end
            push!(lbls, getVariable(fg,frid).label) # fg.v[frid].

        end

        #r = Array{RemoteRef,1}(lenfr)
        #[r[j] = @spawn kde!(pvals[j]) for j in 1:lenfr]
        #[p[j] = fetch(r[j]) for j in 1:lenfr]

        i+=1
        if length(gtvals) > 0
          #gtvals = reshape(gtvals,2,round(Int,length(gtvals)/2))'
          vv[i] = drawHorDens(p, N, gt=gtvals, lbls=lbls)
        else
          vv[i] = drawHorDens(p, N,lbls=lbls)
        end
    end

    #r = Array{RemoteRef,1}(lenfr)
    #[r[j] = @spawn kde!(pvals[j]) for j in 1:lenfr]
    #[p[j] = fetch(r[j]) for j in 1:lenfr]

    i+=1
    if length(gtvals) > 0
      #gtvals = reshape(gtvals,2,round(Int,length(gtvals)/2))'
      vv[i] = drawHorDens(p, N, gt=gtvals, lbls=lbls)
    else
      vv[i] = drawHorDens(p, N,lbls=lbls)
    end
  #
  return vv
end





function saveplot(pl;name="pl",frt=:png,w=25cm,h=25cm,nw=false,fill=true)
  if frt==:png
    Gadfly.draw(PNG(string(name,".png"),w,h),pl)
    # if fill run(`composite $(name).png plB.png $(name).png`) end
    if !nw run(`eog $(name).png`) end
  end
  if frt==:pdf
    Gadfly.draw(PDF(string(name,".pdf"),w,h),pl)
    if !nw run(`evince $(name).pdf`) end
  end
  nothing
end

function animateVertexBelief(FGL::Array{<:AbstractDFG,1}, lbl; nw=false)
  len = length(FGL)
  [saveplot(plotLocalProduct(FG[i],lbl),h=15cm,w=30cm,name="gifs/pl$(i)",nw=true) for i=1:len];
  run(`convert -delay 100 gifs/pl'*'.png result.gif`)
  if !nw run(`eog result.gif`) end
  nothing
end

function fixRotWrapErr!(RT::Array{Float64,1})

  for i in 1:length(RT)
    if RT[i] > pi
      RT[i] = abs(RT[i]-2.0*pi)
    end
  end
  nothing
end

function asyncUniComp(fgl::G, isamdict::Dict{Int,Array{Float64,1}}) where G <: AbstractDFG
  r,rt,lb = computeGraphResiduals(fgl,isamdict);
  fixRotWrapErr!(rt)
  return [sqrt(mean(r.^2));maximum(abs(r));sqrt(mean(rt.^2));maximum(rt)]
end

function unimodalCompare(FGL::Array{<:AbstractDFG,1},isamdict::Dict{Int,Array{Float64,1}})
  len = length(FGL)
  RMS = Float64[]
  MAX = Float64[]
  RMSth = Float64[]
  MAXth = Float64[]

  rr = Future[] #RemoteRef[]

  for fgl in FGL
    push!(rr, remotecall(uppA(),asyncUniComp, fgl, isamdict))
  end

  for r in rr
    err = fetch(r)
    push!(RMS, err[1])
    push!(MAX, err[2])
    push!(RMSth, err[3])
    push!(MAXth, err[4])
  end

  x=0:(len-1)
  df1 = DataFrame(x=x, y=RMS, label="rms")
  df2 = DataFrame(x=x, y=MAX, label="max")
  df3 = DataFrame(x=x, y=RMSth*180.0/pi, label="rmsth")
  df4 = DataFrame(x=x, y=MAXth*180.0/pi, label="maxth")
  df = vcat(df1, df2)
  dfth = vcat(df3,df4)

  return df,dfth
end

function asyncAnalyzeSolution(fgl::G, sym::Symbol) where G <: AbstractDFG
  lbl = string(sym)
  pp, arr, partials = IncrementalInference.localProduct(fgl, lbl)
  lpm = getKDEMax(pp)
  em = getKDEMax(getBelief(getVariable(fgl,lbl)))
  err1 = norm(lpm[1:2]-em[1:2])
  err2 = 0.0
  if lbl[1]=='x'
    err2 = abs(lpm[3]-em[3])
  end
  return [err1;err2]
end

function analyzeSolution(FGL::Array{<: AbstractDFG,1},fggt=Union{})
  len = length(FGL)
  RMS = Float64[]
  MAX = Float64[]
  RMSth = Float64[]
  MAXth = Float64[]
  for fgl in FGL
    xLB, lLB = ls(fgl)
    ERR = Float64[]
    ERRth = Float64[]
    ALB = [xLB;lLB]
    rr = Future[] #RemoteRef[]
    for lbl in ALB
      push!(rr, remotecall(uppA(),asyncAnalyzeSolution, fgl, lbl))
      # err
      # push!(ERR, err[1])
      # if lbl[1]=='x'
      #   push!(ERRth, err[2])
      # end
    end

    idx = 1
    for r in rr
      err = fetch(r)
      push!(ERR, err[1])
      if ALB[idx][1]=='x'
        push!(ERRth, err[2])
      end
      idx += 1
    end
    push!(RMS, sqrt(mean(ERR.^2)))
    push!(MAX, maximum(abs(ERR)))
    push!(RMSth, sqrt(mean(ERRth.^2)))
    push!(MAXth, maximum(ERRth))
  end

  x=0:(len-1)
  df1 = DataFrame(x=x, y=RMS, label="rms")
  df2 = DataFrame(x=x, y=MAX, label="max")
  df3 = DataFrame(x=x, y=RMSth*180.0/pi, label="rmsth")
  df4 = DataFrame(x=x, y=MAXth*180.0/pi, label="maxth")
  df = vcat(df1, df2)
  dfth = vcat(df3,df4)
  return df,dfth
end
# discrete_color_manual(colors...; levels=nothing,order=nothing) is deprecated, use color_discrete_manual(colors...; levels=levels,order=order) instead.

function plotAnalysis(df,dfth)
  return vstack(
  Gadfly.plot(df, x="x", y="y", color="label", Geom.line,
         Scale.discrete_color_manual("red","black")),
  Gadfly.plot(dfth, x="x", y="y", color="label", Geom.line,
        Scale.discrete_color_manual("red","black"))
        )
end

function getAllFGsKDEs(fgD::Array{<: AbstractDFG,1}, vertsym::Symbol)
  ret = Array{BallTreeDensity,1}()
  for i in 1:length(fgD)
    push!(ret, getBelief(getVariable(fgD[i],vertsym)) )
  end
  return ret
end

function plotAllPose2DBeliefs(plots::Array{Gadfly.Compose.Context,1}, fgD::Array{<: AbstractDFG,1})
    ids = sort(collect(keys(fgD[1].v)))
    co = ["black"; "blue"; "green"; "red"; "magenta"; "cyan"; "cyan1"; "cyan2"]
    println(co[1:length(fgD)])
    for i in ids
        getVariable(fgD[1],i).label
        kdes = getAllFGsKDEs(fgD, i)
        push!(plots, plotKDE(  kdes  )) # [kde!(getVal(V)); kde!(getVal(V0))]
    end
    vstackedPlots(plots)
end

# legacy function -- use the array version instead
function plotAllPose2DBeliefs(plots::Array{Gadfly.Compose.Context,1}, fgD::G, fgD0=Union{}) where G <: AbstractDFG
  println("WARNING: drawAllPose2DBeliefs -- legacy function -- use the array version instead.")
  if fgD0 != Union{}
    drawAllPose2DBeliefs(plots, [fgD;fgD0])
  else
    drawAllPose2DBeliefs(plots, [fgD])
  end
end

function plotComicStripLM(fgD::Array{<: AbstractDFG,1})
    comicA = Array{Gadfly.Plot,1}()
    for fgd in fgD
        cv = drawPosesLandms(fgd)
        # cv = drawPoses(fgd)
        push!(comicA,cv)
    end
    hstack(comicA)
end

function plotComicStrip(fgD::Array{<: AbstractDFG,1})
    comicA = Array{Gadfly.Plot,1}()
    for fgd in fgD
        cv = drawPoses(fgd)
        push!(comicA,cv)
    end
    hstack(comicA)
end


function compositeComic(fnc::Function, fgGT, fgA::Array{<: AbstractDFG,1})
    v = Union{}
    @show length(fgA)
    if length(fgA) == 2
        Gadfly.set_default_plot_size(25cm, 10cm)
        v = fnc([fgA[1:2];fgGT])
    elseif length(fgA) == 3
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = vstack(fnc(fgA[1:2])
        ,fnc([fgA[3];fgGT])    )
    elseif length(fgA) == 4
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc([fgA[4];fgGT])    )
    elseif length(fgA) == 7
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc([fgA[7];fgGT])    )
    elseif length(fgA) == 10
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc([fgA[10];fgGT])    )
    elseif length(fgA) == 13
        Gadfly.set_default_plot_size(25cm, 30cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc(fgA[10:12])
        ,fnc([fgA[13];fgGT])    )
    end
    v
end


function compositeComic(fnc::Function, fgA::Array{<: AbstractDFG,1})
    v = Union{}
    @show length(fgA)
    if length(fgA) == 2
        Gadfly.set_default_plot_size(25cm, 10cm)
        v = fnc(fgA[1:2])
    elseif length(fgA) == 3
        Gadfly.set_default_plot_size(25cm, 20cm)
        v = fnc(fgA[1:3])
    elseif length(fgA) == 4
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:2])
        ,fnc(fgA[3:4]))
    elseif length(fgA) == 6
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])  )
    elseif length(fgA) == 9
        Gadfly.set_default_plot_size(25cm, 25cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])    )
    elseif length(fgA) == 12
        Gadfly.set_default_plot_size(25cm, 30cm)
        v = vstack(fnc(fgA[1:3])
        ,fnc(fgA[4:6])
        ,fnc(fgA[7:9])
        ,fnc(fgA[10:12])    )
    end
    v
end

#
#
#
# function spyCliqMat(cliq::Graphs.ExVertex; showmsg=true)
#   mat = deepcopy(getCliqMat(cliq, showmsg=showmsg))
#   # TODO -- add improved visualization here, iter vs skip
#   mat = map(Float64, mat)*2.0.-1.0
#   numlcl = size(IIF.getCliqAssocMat(cliq),1)
#   mat[(numlcl+1):end,:] *= 0.9
#   mat[(numlcl+1):end,:] .-= 0.1
#   numfrtl1 = floor(Int,length(getData(cliq).frontalIDs) + 1)
#   mat[:,numfrtl1:end] *= 0.9
#   mat[:,numfrtl1:end] .-= 0.1
#   @show getData(cliq).itervarIDs
#   @show getData(cliq).directvarIDs
#   @show getData(cliq).msgskipIDs
#   @show getData(cliq).directFrtlMsgIDs
#   @show getData(cliq).directPriorMsgIDs
#   sp = Gadfly.spy(mat)
#   push!(sp.guides, Gadfly.Guide.title("$(cliq.attributes["label"]) || $(cliq.attributes["data"].frontalIDs) :$(cliq.attributes["data"].conditIDs)"))
#   push!(sp.guides, Gadfly.Guide.xlabel("fmcmcs $(cliq.attributes["data"].itervarIDs)"))
#   push!(sp.guides, Gadfly.Guide.ylabel("lcl=$(numlcl) || msg=$(size(getCliqMsgMat(cliq),1))" ))
#   return sp
# end
# function spyCliqMat(bt::BayesTree, lbl::Symbol; showmsg=true)
#   spyCliqMat(whichCliq(bt,lbl), showmsg=showmsg)
# end





function plotPose2Vels(fgl::G, sym::Symbol; coord=nothing) where G <: AbstractDFG
  X = getBelief(getVariable(fgl, sym))
  px = plotKDE(X, dims=[4], title="Velx")
  coord != nothing ? (px.coord = coord) : nothing
  py = plotKDE(X, dims=[5], title="Vely")
  coord != nothing ? (py.coord = coord) : nothing
  hstack(px, py)
end


"""
    $(SIGNATURES)

Analysis function to compare KDE plots between the factor graph centric product of a variable with
current value stored in the factor graph object.
"""
function plotProductVsKDE(fgl::G,
                          sym::Symbol;
                          levels::Int=3,
                          c::Vector{String}=["red";"black"] ) where  G <: AbstractDFG
    #
    plotKDE([IIF.localProduct(fgl, sym)[1], getBelief(getVariable(fgl, sym))], levels=3, c=c)
end





#