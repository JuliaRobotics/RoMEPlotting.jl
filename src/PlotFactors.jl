

function plotFactorMeasurements(dfg::AbstractDFG, fctsym::Symbol, fct::FunctorInferenceType; hdl=[])
  @error "plotFactorMeasurements not implemented yet for $(typeof(fct))."
end

function plotFactorMeasurements(dfg::AbstractDFG, fctsym::Symbol, fct::Pose2Pose2; hdl=[])
  #
  me, me0 = solveFactorMeasurements(dfg, fctsym)

  pt = plotKDE([manikde!(me[1:2,:], Point2);manikde!(me0[1:2,:], Point2)], c=["red";"blue"], legend=["pred";"meas"], levels=3)
  pc = plotKDECircular([manikde!(me[3:3,:], Sphere1);manikde!(me0[3:3,:], Sphere1)], c=["red";"blue"], legend=["pred";"meas"], title="inv. solve, Pose2Pose2, $fctsym")

  push!(hdl, pt)
  push!(hdl, pc)

  hstack(pt, pc)
end


function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol,
                                fct::Pose2Point2BearingRange;
                                hdl=[] )
  #
  me, me0 = solveFactorMeasurements(dfg, fctsym)

  pc = plotKDECircular([manikde!(me[1:1,:], Sphere1);manikde!(me0[1:1,:], Sphere1)], c=["red";"blue"], legend=["pred";"meas"], title="inv. solve, Pose2Point2BearingRange, $fctsym")
  pcl = plotKDE([manikde!(me[1:1,:], ContinuousScalar);manikde!(me0[1:1,:], ContinuousScalar)], c=["red";"blue"], legend=["pred";"meas"], title="unwrapped rotation, should be [-pi,pi)")
  pl = plotKDE([manikde!(me[2:2,:], ContinuousScalar);manikde!(me0[2:2,:], ContinuousScalar)], c=["red";"blue"], legend=["pred";"meas"])

  push!(hdl, pc)
  push!(hdl, pcl)
  push!(hdl, pl)

  hstack(vstack(pc,pcl), pl)
end


"""
    $SIGNATURES

Calculate the "inverse" SLAM solution to compare measured and predicted noise model samples.
"""
function plotFactorMeasurements(dfg::AbstractDFG, fctsym::Symbol;hdl=[])
  fct = getFactorType(dfg, fctsym)
  plotFactorMeasurements(dfg, fctsym, fct, hdl=hdl)
end



function plotFactor(dfg::AbstractDFG, fctsym::Symbol, fct::Pose2Point2Range; hdl=[])
  # variables
  vars = ls(dfg, fctsym)

  # the pose
  pose = intersect(vars, ls(dfg, Pose2))[1]
  poin = intersect(vars, ls(dfg, Point2))[1]

  # calculate predicted range (inverse solve)
  ptss = map(x->getPoints(getKDE(dfg, x)), vars)
  dpts = (ptss[1][1:2,:] - ptss[2][1:2,:]).^2
  pred = sqrt.(sum(dpts, dims=1))
  plpr = Gadfly.plot(x=pred, Geom.histogram(density=true), Theme(default_color=colorant"deepskyblue"), Guide.title("predicted"))

  # measured
  # meas = getSamples(fct, length(pred))
  smps = rand(fct.Z, length(pred))
  plme = Gadfly.plot(x=smps, Geom.histogram(density=true), Theme(default_color=colorant"magenta"), Guide.title("measured"))

  modl = Gadfly.plot(
  Gadfly.layer(x=smps, Geom.histogram(density=true), Theme(default_color=colorant"magenta")),
    Gadfly.layer(x=pred, Geom.histogram(density=true), Theme(default_color=colorant"deepskyblue")),
    Guide.manual_color_key("Legend", ["samples";"predicted"], ["magenta";"deepskyblue"]),
    Guide.title("Pose2Point2Range: Predicted and measurement probability model"),
    Guide.xlabel("range")
  )

  plxy = plotKDE(map(x->getKDE(dfg, x), vars), dims=[1;2], legend=string.(vars), title="Pose2Point2Range: XY plane")
  # plx = plotKDE([getKDE(dfg, vars[1]);], dims=[1;2], legend=[string(vars[1]);], title="Pose2Point2Range: XY plane")
  # ply = plotKDE([getKDE(dfg, vars[2]);], dims=[1;2], legend=[string(vars[2]);], title="Pose2Point2Range: XY plane")

  # mock projection
  tfg = initfg()
  addVariable!(tfg, pose, Pose2)
  addVariable!(tfg, poin, Point2)
  addFactor!(tfg, [pose;poin], fct, autoinit=false)
  manualinit!(tfg, pose, getKDE(dfg,pose))
  manualinit!(tfg, poin, getKDE(dfg,poin))
  apts=Array{Float64,2}(undef, 2, 0)
  loop =true
  while loop
    try
      apts = approxConv(tfg, ls(tfg,pose)[1], poin)
      loop = false
    catch
      @warn "plotFactor Pose2Point2BearingRange non-linear solve fail on approxConv, retrying"
    end
  end
  plhist2 = Gadfly.plot(x=apts[1,:], y=apts[2,:], Geom.histogram2d)
  spc = mean(smps)
  plt = drawPosesLandms(tfg, point_size=5pt, spscale=0.2*spc)
  union!(plhist2.layers, plt.layers)
  pla = plotKDE(manikde!(apts, Point2().manifolds),levels=3,c=["gray50"])

  # plot pose by itself
  # poseplot = plotPose(dfg, pose)
  posepl1 = plotKDE(dfg, pose, dims=[1;2], c=["green"])
  posepl2 = plotKDECircular(marginal(getKDE(dfg, pose), [3]))
  landmpl = plotKDE(dfg, poin, c=["red"])

  # plot handles
  push!(hdl, landmpl)
  push!(hdl, posepl1)
  push!(hdl, posepl2)
  push!(hdl, modl)
  push!(hdl, plme)
  push!(hdl, plpr)
  push!(hdl, plhist2)

  # put plot together
  botr = vstack(modl, plme, plpr) # 4,5,6
  # bot = hstack(plhist2, botr)     # 7,(4,5,6)
  farl = vstack(landmpl, posepl1, posepl2) # 1,2,3

  hstack(farl, botr, plhist2)     # (1,2,3),(4,5,6),7
end


function plotFactor(dfg::AbstractDFG, fctsym::Symbol, fct::Pose2Point2Bearing; hdl=[])
  #

  # variables
  vars = ls(dfg, fctsym)

  # the pose
  pose = intersect(vars, ls(dfg, Pose2))[1]
  poin = intersect(vars, ls(dfg, Point2))[1]

  # basic max estimates
  xp = getKDEMax(getKDE(dfg,pose))
  lp = getKDEMax(getKDE(dfg,poin))

  # convolve the yaw angle with bearing rotation model
  pX = marginal(getKDE(dfg, pose), [3])
  pts = approxConvCircular(pX, fct.bearing)

  # draw plots
  measest = manikde!(pts, Sphere1)

  # inverse solve for predicted bearing
  dx = getPoints(getKDE(dfg, poin))[1,:] - getPoints(getKDE(dfg, pose))[1,:]
  dy = getPoints(getKDE(dfg, poin))[2,:] - getPoints(getKDE(dfg, pose))[2,:]
  pred = reshape(atan.(dy,dx), 1,:)

  ppX = manikde!(pred, Sphere1)

  scal = 0.2*norm(xp[1:2]-lp)
  plcl = plotKDECircular( [measest; ppX], logpdf=true, legend=["Meas. Est.";"Predicted"], radix=scal, scale=0.2*scal, rVo=[xp[1:2];0.0] )

  # plot pose and point by itself
  posepl1 = plotKDE(dfg, pose, dims=[1;2], c=["green"])
  posepl2 = plotKDECircular(marginal(getKDE(dfg, pose), [3]))
  landmpl = plotKDE(dfg, poin, c=["red"])

  tfg = initfg()
  addVariable!(tfg, pose, Pose2)
  addVariable!(tfg, poin, Point2)
  addFactor!(tfg, [pose;poin], fct, autoinit=false)
  manualinit!(tfg, pose, getKDE(dfg,pose))
  manualinit!(tfg, poin, getKDE(dfg,poin))
  plt = drawPosesLandms(tfg, point_size=5pt)

  union!(plt.layers, plcl.layers)
  plt.coord = Coord.cartesian(aspect_ratio=1.0)

  # draw line from pose to point
  push!(plt.layers, Gadfly.layer(x=[xp[1];lp[1]], y=[xp[2];lp[2]], Geom.line)[1])

  # plot handles
  push!(hdl, landmpl)
  push!(hdl, posepl1)
  push!(hdl, posepl2)
  push!(hdl, plcl)
  push!(hdl, plt)

  hstack(vstack(landmpl,posepl1, posepl2), vstack(plcl, plt))
end


function plotFactor(dfg::AbstractDFG, fctsym::Symbol, fct::Pose2Point2BearingRange; hdl=[])
    #
    hdlb = []

    # variables
    vars = ls(dfg, fctsym)

    # the pose
    pose = intersect(vars, ls(dfg, Pose2))[1]
    poin = intersect(vars, ls(dfg, Point2))[1]

    # plot current pose & point
    # pl_poin = plotKDE(dfg, poin, levels=5, c=["red";])
    # pl_pose = plotPose(dfg, pose, c=["black"])

    # project landmark

    # do range model separately too
    pr = Pose2Point2Range(fct.range)
    plotFactor(dfg, fctsym, pr, hdl=hdl)

    br = Pose2Point2Bearing(fct.bearing)
    plotFactor(dfg, fctsym, br, hdl=hdlb)

    newguid = Gadfly.GuideElement[]
    for guid in hdl[4].guides
      if typeof(guid) == Gadfly.Guide.Title
        push!(newguid, Gadfly.Guide.Title("Pose2Point2BearingRange: Predicted and measurement probability model"))
      else
        push!(newguid, guid)
      end
    end
    hdl[4].guides = newguid


    union!(hdl, hdlb)

    # plot a prediction of landmark
    predpoints = approxConv(dfg, fctsym, poin)
    PP = manikde!(predpoints, Point2)
    predLandm = plotKDE([getKDE(dfg, poin); PP], levels=2, c=["red"; "deepskyblue"], legend=["landmark";"predicted"])

    push!(hdl, predLandm)

    # measurement solution
    plotFactorMeasurements(dfg, fctsym, hdl=hdl)

    return hstack(vstack(hdl[1],hdl[2],hdl[3]),vstack(hdl[4],hdl[5],hdl[6]),vstack(predLandm, hdlb[5],hdl[7]), vstack(hdl[14],hdl[15],hdl[16]))
end


function plotFactor(dfg::AbstractDFG, fctsym::Symbol, fct::Pose2Pose2; hdl=[])

  # variables
  fct = getFactor(dfg, fctsym)
  vars = fct._variableOrderSymbols

  pv1 = plotPose(dfg, vars[1], hdl=hdl)
  pv2 = plotPose(dfg, vars[2], hdl=hdl)

  pv12 = plotFactorMeasurements(dfg, fctsym, hdl=hdl)

  vstack(pv1, pv2, pv12)
end


"""
    $SIGNATURES

Return plot of each specific factor.
"""
function plotFactor(dfg::AbstractDFG, fctsym::Symbol; hdl=[])
  plotFactor(dfg, fctsym, getFactorType(dfg, fctsym), hdl=hdl)
end





## Reports


function reportFactors(dfg::AbstractDFG,
                       T::Union{Type{Pose2Pose2}, Type{Pose2Point2BearingRange}},
                       fcts::Vector{Symbol}=ls(dfg, T);
                       path="/tmp/caesar/random/report",
                       show::Bool=true  )
  #
  mkpath(path)

  files = String[]
  for fc in fcts
    file = joinpath(path,"$fc.pdf")
    plotFactor(dfg, fc) |> PDF(file)
    push!(files, file)
  end
  outpath = joinpath(path,"$T.pdf")
  push!(files, outpath)

  run(`pdfunite $files`)
  !show ? nothing : (@async run(`evince $outpath`))
  return outpath
end






#
