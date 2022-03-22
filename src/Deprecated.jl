
## =======================================================================================
## Remove before v0.11
## =======================================================================================

# see JuliaRobotics/Caesar.jl#811
@deprecate plotKDE(w...;kw...) plotBelief(w...;kw...)


## =======================================================================================
## Remove before v0.10
## =======================================================================================


export togglePrtStbLines
export stbPrtLineLayers!
# export drawHorDens # from KDEPlotting
# export plotMCMC, plotPose2DMC!, mcmcPose2D!


# for some reason we still need msgPlots of right size in the global for these functions to work.
# precall drawTreeUpwardMsgs or drawFrontalDens to make this work properly TODO
function vstackedDensities(msgPlots)
  #msgPlots = f(fg, bt) # drawTreeUpwardMsgs
  evalstr = ""
  for i in 1:length(msgPlots)
      evalstr = string(evalstr, ",msgPlots[$(i)]")
  end
  eval(parse(string("vstack(",evalstr[2:end],")")))
end

global DISABLESTBPRTLINES = false

function togglePrtStbLines()
  global DISABLESTBPRTLINES
#   DISABLESTBPRTLINES = !DISABLESTBPRTLINES
# end

## TODO -- you were here with port starboard lines
function stbPrtLineLayers!(pl, Xpp, Ypp, Thpp; l::Real=5.0)
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




#