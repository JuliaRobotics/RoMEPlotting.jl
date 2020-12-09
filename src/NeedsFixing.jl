


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

# """
#     $(SIGNATURES)

# Plot the product of priors and incoming upward messages for variable in clique.

# plotPriorsAtCliq(tree, :x2, :x1[, marg=[1;2]])
# """
# function plotPriorsAtCliq(tree::BayesTree,
#                           lb::Symbol,
#                           cllb::Symbol;
#                           dims::Vector{Int}=Int[],
#                           levels::Int=1,
#                           fill::Bool=false  )
#   #
#   # COLORS = ["black";"red";"green";"blue";"cyan";"deepskyblue"]

#   cliq = whichCliq(tree, cllb)
#   cliqprs = cliq.attributes["debug"].priorprods[1]

#   vidx = 1
#   for lbel in cliqprs.lbls
#     if lbel == lb
#       break;
#     else
#       vidx += 1
#     end
#   end
#   dims = length(dims)>0 ? dims : collect(1:size(cliqprs.prods[vidx].prev,1))

#   arr = BallTreeDensity[]
#   push!(arr, marginal(kde!(cliqprs.prods[vidx].prev), dims)  )
#   push!(arr, marginal(kde!(cliqprs.prods[vidx].product), dims)  )
#   len = length(cliqprs.prods[vidx].potentials)
#   lg = String["p";"n"]
#   i=0
#   for pot in cliqprs.prods[vidx].potentials
#     push!(arr, marginal(pot, dims)  )
#     i+=1
#     push!(lg, cliqprs.prods[vidx].potentialfac[i])
#   end

#   COLORS = getColorsByLength(len+2)
#   # cc = COLORS[1:(len+2)]
#   plotKDE(arr, c=COLORS, legend=lg, levels=levels, fill=fill );
# end



# function plotMCMC(treel::BayesTree,
#                   lbll::Symbol;
#                   delay::Int=200,
#                   show::Bool=true,
#                   w=20cm, h=15cm,
#                   levels::Int=1,
#                   dims=nothing  )
#   #
#   cliq = whichCliq(treel, string(lbll))
#   cliqdbg = cliq.attributes["debug"]

#   vidx = 1
#   for lb in cliqdbg.mcmc[1].lbls
#     if lb == lbll
#       break;
#     else
#       vidx += 1
#     end
#   end

#   tmpfilepath = joinpath(dirname(@__FILE__),"tmpimgs")
#   ARR = BallTreeDensity[]
#   COLORS = getColorsByLength()
#   # COLORS = ["black";"red";"green";"blue";"cyan";"deepskyblue";"magenta"]
#   for i in 1:length(cliqdbg.mcmc)
#     ppr = kde!(cliqdbg.mcmc[i].prods[vidx].prev)
#     ppp = kde!(cliqdbg.mcmc[i].prods[vidx].product)
#     ARR = [ARR;ppr;ppr;cliqdbg.mcmc[i].prods[vidx].potentials]
#   end
#   rangeV = getKDERange(ARR)
#   ppp = nothing
#   for i in 1:length(cliqdbg.mcmc)
#     ppr = kde!(cliqdbg.mcmc[i].prods[vidx].prev)
#     ppp = kde!(cliqdbg.mcmc[i].prods[vidx].product)
#     arr = [ppr;ppp;cliqdbg.mcmc[i].prods[vidx].potentials]
#     len = length(cliqdbg.mcmc[i].prods[vidx].potentials)
#     lg = String["p";"n";cliqdbg.mcmc[i].prods[vidx].potentialfac] #map(string, 1:len)]
#     COLORS_l = getColorsByLength(len+2)
#     cc = plotKDE(arr, c=COLORS_l, legend=lg, levels=levels, fill=true, axis=rangeV, dims=dims );
#     Gadfly.draw(PNG(joinpath(tmpfilepath,"$(string(lbll))mcmc$(i).png"),w,h),cc)
#   end
#   # draw initial and final result
#   pp0 = kde!(cliqdbg.mcmc[1].prods[vidx].prev)
#   i = 0
#   cc = plotKDE([pp0], c=[COLORS[1]], legend=["0"], levels=levels, fill=true, axis=rangeV, dims=dims );
#   Gadfly.draw(PNG(joinpath(tmpfilepath,"$(string(lbll))mcmc$(i).png"),w,h),cc)

#   i = length(cliqdbg.mcmc)+1
#   cc = plotKDE([pp0;ppp], c=COLORS[1:2], legend=["0";"n"], levels=levels, fill=true, axis=rangeV, dims=dims );
#   Gadfly.draw(PNG(joinpath(tmpfilepath,"$(string(lbll))mcmc$(i).png"),w,h),cc)
#   # generate output
#   @warn "Maintenance required, run manually:  \"convert -delay $(delay) $(tmpfilepath)/$(string(lbll))mcmc*.png $(tmpfilepath)/$(string(lbll))mcmc.gif\""
#   # run(`convert -delay $(delay) $(tmpfilepath)/$(string(lbll))mcmc*.png $(tmpfilepath)/$(string(lbll))mcmc.gif`)
#   # !show ? nothing : (@async run(`eog $(tmpfilepath)/$(string(lbll))mcmc.gif`) )
#   # return "$(tmpfilepath)/$(string(lbll))mcmc.gif"
# end

# function plotOneMC!(cliqMC::CliqGibbsMC, minmax, mcmc=0; offs=2.0)

#     if mcmc>minmax[4]
#       minmax[4]=mcmc
#     end

#     i = 0.0
#     for prod in cliqMC.prods
#         prodVal = kde!(prod.product,"lcv")
#         plotKDEProd!([prodVal;prod.potentials],minmax, h=-i*offs, mcmc=mcmc)
#         i += 1.0
#     end

# end

# function plotMCMCDebug(cliq; offs=2.0)
#     println("$(cliq.attributes["label"])")
#     cliqDbg = cliq.attributes["data"].debug
#     MCd = 100.0/length(cliqDbg.mcmc)
#     MCo=0.0
#     minmax=[99999,-99999,-10,-99999.0]
#     for mc in cliqDbg.mcmc
#         drawOneMC!(mc, minmax, MCo, offs=offs)
#         MCo+=MCd
#     end

#     n = 3
#     y = linspace(minmax[1], minmax[2], n)
#     x = linspace(minmax[3],minmax[4],n)
#     xgrid = repmat(x',n,1)
#     ygrid = repmat(y,1,n)
#     z = zeros(n,n)

#     for i in 1:length(cliqDbg.mcmc[1].prods)
#       surf(xgrid,ygrid,z-(i-1)*offs,alpha=0.04, linewidth=0.0)
#     end
# end
