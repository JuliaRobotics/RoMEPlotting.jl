
"""
    $SIGNATURES

Return plot of each specific factor.
"""
function plotFactor(dfg::G, fctsym::Symbol; hdl=[]) where G <: AbstractDFG
  plotFactor(dfg, fctsym, getFactorType(dfg, fctsym), hdl=hdl)
end

function plotFactor(dfg::G, fctsym::Symbol, fct::Pose2Point2Range; hdl=[]) where G <: AbstractDFG
  # variables
  vars = ls(dfg, fctsym)

  # the pose
  pose = intersect(vars, ls(dfg, Pose2))[1]
  poin = intersect(vars, ls(dfg, Point2))[1]

  # calculate predicted range
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
  apts = approxConv(tfg, ls(tfg,pose)[1], poin)
  plhist2 = Gadfly.plot(x=apts[1,:], y=apts[2,:], Geom.histogram2d)
  plt = drawPosesLandms(tfg, point_size=5pt)
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
  botr = vstack(modl, plme, plpr)
  bot = hstack(plhist2, botr) # vstack(plx,ply)
  farl = vstack(landmpl, posepl1, posepl2)
  # leftcol = vstack(bearpl, plhist2)
  # bothalf = hstack(farl, plhist2)

  hstack(farl, botr, plhist2)
end
