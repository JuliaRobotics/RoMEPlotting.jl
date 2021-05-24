# new file for plotVariable related functionality

export plotSLAM2DSolveKeys

# function plotSolveKey(dfg::AbstractDFG,
#                       refSym::Symbol,
#                       refKey::Symbol,
#                       tstSym::Symbol,
#                       tstKey::Symbol;
#                       bw::AbstractVector{<:Real}=[0.001;],
#                       plotVarHack::Function=plotPose  )
#   #
#   pts, fctT = deconvSolveKey(dfg, refSym, refKey, tstSym, tstKey)
#   Xref = manikde!(pts[2],fctT)
#   Xtst = manikde!(pts[1],fctT)
#   mmdDist = mmd(Xref, Xtst, fctT, bw=bw)

#   # FIXME not all variables are Poses, so this is really hacky 
#   plotVarHack(getManifold(fctT)(), [Xref, Xtst], "$fctT\n$(refSym)_$(refKey) <-> $(tstSym)_$tstKey\nmmd=$(round(mmdDist,digits=6))", c=["red","blue"], legend=["Identity";"Delta"])
# end


"""
    $SIGNATURES

Analyzing repeat solves is easier if you can animate repeat solves.

Notes
- Repeat solves can be stored using `solveTree!(fg, storeOld=true)`. 
- Experimental, see DFG #641

Example

```julia
using Images, Caesar, RoMEPlotting

fg = generateCanonicalFG_Hexagonal()

for i in 1:10
  solveTree!(fg, storeOld=true)
end

# generate all the plots in RAM
plts = plotSLAM2DSolveKeys(fg)

# export all images to a video
imgs = convert.(Matrix{RGB},plts)
Caesar.writevideo("/tmp/test.avi", imgs)
@async run(`totem /tmp/test.avi`)
```

Related

writevideo, plotSolveKey, listSolveKeys
"""
function plotSLAM2DSolveKeys( dfg::AbstractDFG,
                              pattern::Regex=r"default_\d+";
                              xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing,
                              contour=false )
  #
  kys = listSolveKeys(dfg, filterSolveKeys=pattern) |> collect |> sortDFG
  kys .|> x -> plotSLAM2D(dfg, solveKey=x, contour=contour, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
end



#