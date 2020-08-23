# functions that might need to be loaded multiple times (dont add docs here)


function plotFactorMeasurements(dfg::AbstractDFG,
                                fctsym::Symbol,
                                fct::ExtendedPose2Pose2Types;
                                hdl=[],
                                dist::Vector{Float64}=[0.0;]  )
  #
  asMeas, asPred = solveFactorMeasurements(dfg, fctsym)

  plotFactorValues(asMeas, asPred, Pose2Pose2, fctsym=fctsym,hdl=hdl,dist=dist)
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



function reportFactors(dfg::AbstractDFG,
                       T::PlotTypesPose2,
                       fcts::Vector{Symbol}=ls(dfg, T);
                       prefix::AbstractString="",
                       filepath=joinpath(getSolverParams(dfg).logpath, "$prefix"*getTimeEasy()*"_$T.pdf"),
                       show::Bool=true,
                       pdfWidth=20cm,
                       pdfHeight=30cm)
  #
  ss = split(filepath, '/')
  path = joinpath("/", joinpath(ss[1:(end-1)]...), "tmp")
  mkpath(path)
  alldists= Vector{Float64}()

  files = String[]
  ndist = Float64[0.0;]
  for fc in fcts
    if isMultihypo(getFactor(dfg, fc))
      # skip this factor
      continue
    end
    file = joinpath(path,"$fc.pdf")
    ndist[1] = 0.0
    plotFactor(dfg, fc, dist=ndist) |> PDF(file, pdfWidth, pdfHeight)
    push!(files, file)
    push!(alldists, ndist[1])
  end
  fileord = sortperm(alldists, rev=true)
  files = files[fileord]
  push!(files, filepath)

  2 < length(files) ? run(`pdfunite $files`) : nothing
  !show ? nothing : (@async run(`evince $filepath`))
  return filepath
end


#
