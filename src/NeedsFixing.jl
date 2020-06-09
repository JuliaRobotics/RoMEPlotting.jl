


"""
    $(SIGNATURES)

Plot absolute values of variables and measurement surrounding fsym factor.
"""
function plotKDEofnc(fgl::AbstractDFG,
                     fsym::Symbol;
                     solveKey::Symbol=:default,
                     marg=nothing,
                     N::Int=100 )
  #
  fnc = nothing
  if haskey(fgl.fIDs, fsym)
    fnc = getfnctype( fgl, fgl.fIDs[fsym] )
  else
    error("fIDs doesn't have $(fsym)")
  end
  p = kde!( getSample(fnc, N)[1]  )
  # mmarg = length(marg) > 0 ? marg : collect(1:Ndim(p))
  plotKDE(fgl, lsf(fgl, fsym), solveKey=solveKey, addt=[p], marg=marg)
end

"""
    $(SIGNATURES)

Try plot relative values of variables and measurement surrounding fsym factor.
"""
function plotKDEresiduals(fgl::AbstractDFG,
                          fsym::Symbol;
                          N::Int=100,
                          levels::Int=3,
                          dims=nothing,
                          fill=false  )
  #
  # COLORS = ["black";"red";"green";"blue";"cyan";"deepskyblue"]
  fnc = getfnctype( fgl, fgl.fIDs[fsym] )
  # @show sxi = lsf(fgl, fsym)[1]
  # @show sxj = lsf(fgl, fsym)[2]
  fct = getFactor(fgl, fsym)
  @show sxi = getData(fct).fncargvID[1]
  @show sxj = getData(fct).fncargvID[2]
  xi = getVal(fgl, sxi)
  xj = getVal(fgl, sxj)
  measM = getSample(fnc, N)
  meas = length(measM) == 1 ? (0*measM[1], ) : (0*measM[1], measM[2])
  @show size(meas[1])
  d = size(measM[1],1)
  RES = zeros(d,N)
  for i in 1:N
    res = zeros(d)
    fnc(res, i, meas, xi, xj)
    RES[:,i] = res
    if length(measM) > 1
      if measM[2][i] == 0
        RES[:,i] = 0.5*randn(d)
      end
    end
  end
  pR = kde!(RES)
  pM = kde!(measM[1])
  fncvar = getfnctype(fct)

  COLORS = getColorsByLength()
  plotKDE([pR;pM], c=COLORS[1:2], dims=dims,levels=3, legend=["residual";"model"], fill=fill, title=string(fsym, ", ", fncvar))
end

"""
    $(SIGNATURES)

Plot the product of priors and incoming upward messages for variable in clique.

plotPriorsAtCliq(tree, :x2, :x1[, marg=[1;2]])
"""
function plotPriorsAtCliq(tree::BayesTree,
                          lb::Symbol,
                          cllb::Symbol;
                          dims::Vector{Int}=Int[],
                          levels::Int=1,
                          fill::Bool=false  )
  #
  # COLORS = ["black";"red";"green";"blue";"cyan";"deepskyblue"]

  cliq = whichCliq(tree, cllb)
  cliqprs = cliq.attributes["debug"].priorprods[1]

  vidx = 1
  for lbel in cliqprs.lbls
    if lbel == lb
      break;
    else
      vidx += 1
    end
  end
  dims = length(dims)>0 ? dims : collect(1:size(cliqprs.prods[vidx].prev,1))

  arr = BallTreeDensity[]
  push!(arr, marginal(kde!(cliqprs.prods[vidx].prev), dims)  )
  push!(arr, marginal(kde!(cliqprs.prods[vidx].product), dims)  )
  len = length(cliqprs.prods[vidx].potentials)
  lg = String["p";"n"]
  i=0
  for pot in cliqprs.prods[vidx].potentials
    push!(arr, marginal(pot, dims)  )
    i+=1
    push!(lg, cliqprs.prods[vidx].potentialfac[i])
  end

  COLORS = getColorsByLength(len+2)
  # cc = COLORS[1:(len+2)]
  plotKDE(arr, c=COLORS, legend=lg, levels=levels, fill=fill );
end
