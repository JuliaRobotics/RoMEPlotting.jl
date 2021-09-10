
# import RoMEPlotting: plotFactor, reportFactors, plotFactorMeasurements, getTimeEasy

"""
    $SIGNATURES

Plot Pose2Pose2 measurement and predicted values in factor centric manner -- i.e. used with inverse solve.
"""
function plotFactorValues(asMeasured::AbstractMatrix{<:Real},
                          asPredicted::AbstractMatrix{<:Real},
                          fct::Type{T}=Pose2Pose2;
                          title="inv. solve",
                          fctsym="",
                          hdl=[],
                          dist::Vector{Float64}=zeros(1) ) where {T <: FunctorInferenceType}
  #
  PPg = manikde!(asMeasured[1:2,:], Point2)
  PP  = manikde!(asPredicted[1:2,:], Point2)
  # dist[1] = minimum(abs.([kld(PPg, PP)[1]; kld(PP, PPg)[1]]))
  mmd!(dist, asPredicted, asMeasured, SE2_Manifold)
  pt = plotKDE([PPg;PP], c=["blue";"red"], legend=["meas";"pred";], levels=3, title="$title $fctsym,\nmmd(..)=$(round(dist[1],digits=6))")

  pc = plotKDECircular([manikde!(asMeasured[3:3,:], Sphere1);manikde!(asPredicted[3:3,:], Sphere1);], c=["blue";"red";], legend=["meas";"pred";], title="$title $fctsym,\n$(T)")

  push!(hdl, pt)
  push!(hdl, pc)

  hstack(pt, pc)
end


"""
    $SIGNATURES

Calculate the "inverse" SLAM solution to compare measured and predicted noise model samples.
"""
function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol,
                                fct::FunctorInferenceType;
                                hdl=[],
                                dist::Vector{Float64}=zeros(1) )
  #
  @error "plotFactorMeasurements not implemented yet for $(typeof(fct))."
end



function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol,
                                fct::Pose2Point2BearingRange;
                                hdl=[],
                                dist::Vector{Float64}=[0.0;] )
  #
  asMeas, asPred = solveFactorMeasurements(dfg, fctsym)

  pc = plotKDECircular([manikde!(asPred[1:1,:], Sphere1);manikde!(asMeas[1:1,:], Sphere1)], c=["blue";"red";], legend=["meas";"pred";], title="inv. solve, $fctsym,\nPose2Point2BearingRange")
  pcl = plotKDE([manikde!(asPred[1:1,:], ContinuousScalar);manikde!(asMeas[1:1,:], ContinuousScalar)], c=["blue";"red";], legend=["meas";"pred";], title="unwrapped rotation, should be [-pi,pi)")
  pl = plotKDE([manikde!(asPred[2:2,:], ContinuousScalar);manikde!(asMeas[2:2,:], ContinuousScalar)], c=["blue";"red";], legend=["meas";"pred";])

  push!(hdl, pc)
  push!(hdl, pcl)
  push!(hdl, pl)

  hstack(vstack(pc,pcl), pl)
end


function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol;
                                hdl=[],
                                dist::Vector{Float64}=[0.0;]  )
  #
  fct = getFactorType(dfg, fctsym)
  plotFactorMeasurements(dfg, fctsym, fct, hdl=hdl, dist=dist)
end

"""
    $SIGNATURES

Return plot of each specific factor.
"""
function plotFactor(dfg::AbstractDFG,
                    fctsym::Symbol,
                    fct::Pose2Point2Range;
                    hdl=[],
                    dist::Vector{Float64}=Float64[0.0;]  )
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

  # measure kl divergence
  PP = manikde!(smps, ContinuousScalar)
  PPg = manikde!(pred, ContinuousScalar)
  dist[1] = minimum(abs.([kld(PPg, PP)[1]; kld(PP, PPg)[1]]))

  modl = Gadfly.plot(
  Gadfly.layer(x=smps, Geom.histogram(density=true), Theme(default_color=colorant"magenta")),
    Gadfly.layer(x=pred, Geom.histogram(density=true), Theme(default_color=colorant"deepskyblue")),
    Guide.manual_color_key("Legend", ["samples";"predicted"], ["magenta";"deepskyblue"]),
    Guide.title("Pose2Point2Range $fctsym\nmin(|kld...|)=$(round(dist[1],digits=3))"),
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


function plotFactor(dfg::AbstractDFG,
                    fctsym::Symbol,
                    fct::Pose2Point2Bearing;
                    hdl=[],
                    dist::Vector{Float64}=Float64[0.0;]  )
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
  pX = marginal(getBelief(dfg, pose), [3])
  pts = approxConvCircular(pX, fct.bearing) # TODO fix by using approxConvBelief instead

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


function plotFactor(dfg::AbstractDFG,
                    fctsym::Symbol,
                    fct::Pose2Point2BearingRange;
                    hdl=[],
                    dist::Vector{Float64}=Float64[0.0;]  )
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
        push!(newguid, Gadfly.Guide.Title("Pose2Point2BearingRange $fctsym"))
      else
        push!(newguid, guid)
      end
    end
    hdl[4].guides = newguid


    union!(hdl, hdlb)

    # plot a prediction of landmark
    predpoints = approxConv(dfg, fctsym, poin)
    PP = manikde!(predpoints, Point2)
    PPg= getKDE(dfg, poin)
    dist[1] = minimum(abs.([kld(PPg, PP)[1]; kld(PP, PPg)[1]]))
    predLandm = plotKDE([PPg; PP], levels=2, c=["red"; "deepskyblue"], legend=["landmark";"predicted"], title="min(|kld...|)=$(round(dist[1],digits=3))")

    push!(hdl, predLandm)

    # AMP.mmd!(dist, PPg, predpoints)

    # measurement solution
    plotFactorMeasurements(dfg, fctsym, hdl=hdl)

    return hstack(vstack(hdl[1],hdl[2],hdl[3]),vstack(hdl[4],hdl[5],hdl[6]),vstack(predLandm, hdlb[5],hdl[7]), vstack(hdl[14],hdl[15],hdl[16]))
end


function plotFactor(dfg::AbstractDFG,
                    fctsym::Symbol;
                    hdl=[],
                    dist::Vector{Float64}=Float64[0.0;]  )
  #
  plotFactor(dfg, fctsym, getFactorType(dfg, fctsym), hdl=hdl, dist=dist)
end





## Reports

getTimeEasy() = split(split("$(now())", 'T')[end],'.')[1]

# import RoMEPlotting: getTimeEasy, reportFactors
# import ApproxManifoldProducts: mmd!

## Also in IIF >v0.7.9
import IncrementalInference: isMultihypo
isMultihypo(fct) = isa(solverData(fct).fnc.hypotheses, Distribution)


"""
$SIGNATURES

Methods to generate plots of all requested factor user supplied vs. current predicted measurement -- i.e. inverse solve.
"""
reportFactors() = @info "Empty function to support both automated documentation and conditional features (@require)."





#
