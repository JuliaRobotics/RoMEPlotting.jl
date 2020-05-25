# functions that might need to be loaded multiple times (dont add docs here)


function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol,
                                fct::ExtendedPose2Pose2Types;
                                hdl=[],
                                dist::Vector{Float64}=[0.0;]  )
  #
  me, me0 = solveFactorMeasurements(dfg, fctsym)

  PP = manikde!(me[1:2,:], Point2)
  PPg = manikde!(me0[1:2,:], Point2)
  dist[1] = minimum(abs.([kld(PPg, PP)[1]; kld(PP, PPg)[1]]))
  pt = plotKDE([PP;PPg], c=["red";"blue"], legend=["pred";"meas"], levels=3, title="inv. solve, $fctsym,\nmin(|kld(..)|)=$(round(dist[1],digits=3))")

  pc = plotKDECircular([manikde!(me[3:3,:], Sphere1);manikde!(me0[3:3,:], Sphere1)], c=["red";"blue"], legend=["pred";"meas"], title="inv. solve, $fctsym,\n$(typeof(fct))")


  push!(hdl, pt)
  push!(hdl, pc)

  hstack(pt, pc)
end


function plotFactor(dfg::AbstractDFG,
                    fctsym::Symbol,
                    fct::ExtendedPose2Pose2Types;
                    hdl=[],
                    dist::Vector{Float64}=Float64[0.0;]  )

  # variables
  fct = getFactor(dfg, fctsym)
  vars = fct._variableOrderSymbols

  pv1 = plotPose(dfg, vars[1], hdl=hdl)
  pv2 = plotPose(dfg, vars[2], hdl=hdl)

  pv12 = plotFactorMeasurements(dfg, fctsym, hdl=hdl, dist=dist)

  vstack(pv1, pv2, pv12)
end
