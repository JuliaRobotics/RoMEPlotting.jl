# new file for plotVariable related functionality

export plotSolveKey

function plotSolveKey(dfg::AbstractDFG,
                      refSym::Symbol,
                      refKey::Symbol,
                      tstSym::Symbol,
                      tstKey::Symbol;
                      bw::AbstractVector{<:Real}=[0.001;],
                      plotVarHack::Function=plotPose  )
  #
  pts, fctT = deconvSolveKey(dfg, refSym, refKey, tstSym, tstKey)
  Xref = manikde!(pts[2],fctT)
  Xtst = manikde!(pts[1],fctT)
  mmdDist = mmd(Xref, Xtst, fctT, bw=bw)

  # FXIME not all variables are Poses, so this is really hacky 
  plotVarHack(getDomain(fctT)(), [Xref, Xtst], "$fctT\n$(refSym)_$(refKey) <-> $(tstSym)_$tstKey\nmmd=$(round(mmdDist,digits=6))", c=["red","blue"], legend=["Identity";"Delta"])
end




#