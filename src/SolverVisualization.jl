

"""
    $(SIGNATURES)

A peneric KDE plotting function that allows marginals of higher dimensional beliefs and various keyword options.

Example:
```julia
p = kde!(randn(3,100))

plotKDE(p)
plotKDE(p, dims=[1;2], levels=3)
plotKDE(p, dims=[1])

q = kde!(5*randn(3,100))
plotKDE([p;q])
plotKDE([p;q], dims=[1;2], levels=3)
plotKDE([p;q], dims=[1])
```
"""
function plotKDE( fgl::AbstractDFG,
                  sym::Symbol;
                  solveKey::Symbol=:default,
                  dims=nothing,
                  title="",
                  levels::Int=5,
                  fill::Bool=false,
                  layers::Bool=false,
                  c=nothing,
                  overlay=nothing  )
  #
  p = getKDE(getVariable(fgl,sym), solveKey)
  # mmarg = length(marg) > 0 ? marg : collect(1:Ndim(p))
  # mp = marginal(p,mmarg)
  plotKDE(p, levels=levels, dims=dims, title=string(sym, "  ", title), fill=fill, layers=layers, c=c, overlay=overlay )
end
function plotKDE( fgl::AbstractDFG,
                  syms::Vector{Symbol};
                  solveKey::Symbol=:default,
                  addt::Vector{BallTreeDensity}=BallTreeDensity[],
                  dims=nothing,
                  title=nothing,
                  levels=3,
                  layers::Bool=false,
                  c=getColorsByLength(length(addt)),
                  overlay=nothing  )
  #
  # TODO -- consider automated rotisary of color
  # colors = ["black";"red";"green";"blue";"cyan";"deepskyblue"; "yellow"]
  # COLORS = repmat(colors, 10)
  # COLORS = getColorsByLength(length(syms))
  MP = BallTreeDensity[]
  LEG = String[]
  # mmarg = Int[]
  for sym in syms
    p = getBelief(getVariable(fgl,sym), solveKey)
    # mmarg = length(marg) > 0 ? marg : collect(1:Ndim(p))
    # mp = marginal(p,mmarg)
    push!(MP, p)
    push!(LEG, string(sym))
  end
  for p in addt
    # mp = marginal(p,mmarg)
    push!(MP, p)
    push!(LEG, "add")
  end
  plotKDE(MP, c=c, levels=levels, dims=dims, legend=LEG, title=title, layers=layers, overlay=overlay)
end



"""
    $(SIGNATURES)

Draw the upward belief from clique `cllb::Symbol` message on variable `lb::Symbol`.

Example:

```julia
plotUpMsgsAtCliq(tree, :x2, :x1)
```

Related

plotKDE, getUpMsgs
"""
function plotUpMsgsAtCliq(treel::AbstractBayesTree,
                          cllb::Symbol,
                          lb::Symbol;
                          show::Bool=true,
                          w=20cm, h=15cm,
                          levels::Int=1,
                          dims::Union{Vector{Int}, Nothing}=nothing )
  #

  cliq = getClique(treel, cllb)
  cliqoutmsg = getUpMsgs(cliq)
  bel = convert(BallTreeDensity, cliqoutmsg.belief[lb])
  plotKDE(bel, dims=dims)
end


function plotTreeUpwardMsgs(fgl::G,
                            bt::AbstractBayesTree;
                            N=300 ) where G <: AbstractDFG
  #
  len = length(bt.cliques)-1
  vv = Array{Gadfly.Compose.Context,1}(undef, len)
  #r = Array{RemoteRef,1}(len)
  i = 0
  for cliq in bt.cliques
      if cliq[1] == 1 println("No upward msgs from root."); continue; end
      @show cliq[2].attributes["label"]
      i+=1
      vv[i] = drawHorDens(fgl, cliq[2].attributes["debug"].outmsg.p, N)
  end
  vv
end


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



"""
    $(SIGNATURES)

Plot the proposal belief from neighboring factors to `lbl` in the factor graph (ignoring Bayes tree representation),
and show with new product approximation for reference.

DevNotes
- TODO, standardize around ::MIME="image/svg", see JuliaRobotics/DistributedFactorGraphs.jl#640
"""
function plotLocalProduct(fgl::AbstractDFG,
                          lbl::Symbol;
                          solveKey::Symbol=:default,
                          N::Int=100,
                          dims::Vector{Int}=Int[],
                          levels::Int=1,
                          show::Bool=false,
                          dirpath="/tmp/",
                          mimetype::AbstractString="svg",
                          sidelength=20cm,
                          title::String="Local product ($solveKey), "  )
  #
  @warn "not showing partial constraints, but included in the product"
  arr = Array{BallTreeDensity,1}()
  lbls = String[]
  push!(arr, getBelief(getVariable(fgl, lbl), solveKey))
  push!(lbls, "curr")
  pl = nothing
  pp, parr, partials, lb = IncrementalInference.localProduct(fgl, lbl, N=N, solveKey=solveKey)

  # helper functions
  function plotDirectProducts()
      if pp != parr[1]
        push!(arr,pp)
        push!(lbls, "prod")
        for a in parr
          push!(arr, a)
        end
        @show lb, lbls
        lbls = union(lbls, string.(lb))
      end
      dims = length(dims) > 0 ? dims : collect(1:Ndim(pp))
      colors = getColorsByLength(length(arr))
      plotKDE(arr, dims=dims, levels=levels, c=colors, title=string(title,lbl), legend=string.(lbls)) #
  end
  function plotPartialProducts()
      # stack 1d plots to accomodate all the partials
      PL = []
      lbls = String["prod";"curr";string.(lb)]
      pdims = sort(collect(keys(partials)))
      for dimn in pdims
        vals = partials[dimn]
        proddim = marginal(pp, [dimn])
        colors = getColorsByLength(length(vals)+2)
        pl = plotKDE([proddim;getBelief(getVariable(fgl, lbl),solveKey);vals], dims=[1;], levels=levels, c=colors, title=string("Local product, dim=$(dimn), ",lbl))
        push!(PL, pl)
      end
      Gadfly.vstack(PL...)
  end

  if length(parr) > 0 && length(partials) == 0
    pl = plotDirectProducts()
  elseif length(parr) == 0 && length(partials) > 0
    pl = plotPartialProducts()
  else
    return error("plotLocalProduct not built for lengths parr, partials = $(length(parr)), $(length(partials)) yet.")
  end

  # now let's export:
  backend = getfield(Gadfly, Symbol(uppercase(mimetype)))
  Gadfly.draw(backend(dirpath*"test_$(lbl).$(mimetype)",sidelength,0.75*sidelength), pl)
  driver = mimetype in ["pdf"] ? "evince" : "eog"
  show ? (@async run(`$driver $(dirpath)test_$(lbl).$(mimetype)`)) : nothing

  return pl
end


"""
    $(SIGNATURES)

Plot---for the cylindrical manifold only---the proposal belief from neighboring factors to `lbl` in the factor graph (ignoring Bayes tree representation),
and show with new product approximation for reference.  The linear and circular belief products are returned as a tuple.
"""
function plotLocalProductCylinder(fgl::G,
                                  lbl::Symbol;
                                  N::Int=100,
                                  levels::Int=1,
                                  show=false,
                                  dirpath="/tmp/",
                                  mimetype::AbstractString="svg",
                                  sidelength=30cm,
                                  scale::Float64=0.2 ) where G <: AbstractDFG
  #
  @warn "not showing partial constraints, but included in the product"
  arr = Array{BallTreeDensity,1}()
  carr = Array{BallTreeDensity,1}()
  lbls = String[]
  push!(arr, getKDE(getVariable(fgl, lbl)))
  push!(lbls, "curr")
  pl = nothing
  plc = nothing
  pp, parr, partials, lb = IncrementalInference.localProduct(fgl, lbl, N=N)
  if length(parr) > 0 && length(partials) == 0
    if pp != parr[1]
      push!(arr, marginal(pp,[1]))
      push!(lbls, "prod")
      push!(carr, marginal(pp, [2]))
      for a in parr
        push!(arr, marginal(a,[1]))
        push!(carr, marginal(a,[2]))
      end
      lbls = union(lbls, lb)
    end
    colors = getColorsByLength(length(arr))
    pll = plotKDE(arr, dims=[1;], levels=levels, c=colors, legend=lbls, title=string("Local product, ",lbl))
    plc = AMP.plotKDECircular(carr, c=colors, scale=scale)
    # (pll, plc)
  elseif length(parr) == 0 && length(partials) > 0
    # stack 1d plots to accomodate all the partials
    PL = []
    PLC = []
    lbls = ["prod";"curr";lb]
    pdims = sort(collect(keys(partials)))
    for dimn in pdims
      vals = partials[dimn]
      proddim = marginal(pp, [dimn])
      colors = getColorsByLength(length(vals)+2)
      if dimn == 1
        pll = plotKDE([proddim;getKDE(getVariable(fgl, lbl));vals], dims=[1;], levels=levels, c=colors, title=string("Local product, dim=$(dimn), ",lbl))
        push!(PL, pl)
      else
        plc = AMP.plotKDECircular([proddim; marginal(getKDE(getVariable(fgl, lbl)),[2]);vals], c=colors, scale=scale)
        push!(PLC, plc)
      end
    end
    pl = Gadfly.vstack(PL..., PLC...)
  else
    return error("plotLocalProduct not built for lengths parr, partials = $(length(parr)), $(length(partials)) yet.")
  end

  # now let's export:
  # backend = getfield(Gadfly, Symbol(uppercase(mimetype)))
  # Gadfly.draw(backend(dirpath*"test_$(lbl).$(mimetype)",sidelength,sidelength), pl)
  # driver = mimetype in ["pdf"] ? "evince" : "eog"
  # show ? (@async run(`$driver $(dirpath)test_$(lbl).$(mimetype)`)) : nothing

  return pll, plc
end


"""
    $SIGNATURES

Plot the proposal belief from neighboring factors to `lbl` in the factor graph (ignoring Bayes tree representation),
and show with new product approximation for reference. String version is obsolete and will be deprecated.
"""
plotLocalProduct(fgl::G, lbl::T; N::Int=100, dims::Vector{Int}=Int[]) where {G <: AbstractDFG, T <: AbstractString} = plotLocalProduct(fgl, Symbol(lbl), N=N, dims=dims)

"""
    $SIGNATURES

Project (convolve) to and take product of variable in Bayes/Junction tree.

Notes
- assume cliq and var sym the same, unless both specified.
- `cliqsym` defines a frontal variable of a clique.
"""
function plotTreeProductUp(fgl::G,
                           treel::AbstractBayesTree,
                           cliqsym::Symbol,
                           varsym::Symbol=cliqsym;
                           levels::Int=1,
                           dims::Vector{Int}=Int[]  ) where G <: AbstractDFG
  #
  # build a subgraph copy of clique
  cliq = getClique(treel, cliqsym)
  syms = getCliqAllVarIds(cliq)
  subfg = buildSubgraph(fgl, syms)

  # add upward messages to subgraph
  msgs = fetchMsgsUpChildren(treel,cliq, TreeBelief)
  # @show typeof(msgs)
  addMsgFactors!.(subfg, msgs, IIF.UpwardPass)

  # predictBelief
  # stuff = treeProductUp(fgl, treel, cliqsym, varsym)
  # plotKDE(manikde!(stuff[1], getManifolds(fgl, varsym)))
  cllbl = cliq.attributes["label"]
  return plotLocalProduct(subfg, varsym, title="Tree Up $(cllbl) | ", levels=levels, dims=dims)
end


function plotTreeProductDown( fgl::AbstractDFG,
                              treel::AbstractBayesTree,
                              cliqsym::Symbol,
                              varsym::Symbol=cliqsym;
                              levels::Int=1  )
  #
  # build a subgraph copy of clique
  cliq = whichCliq(treel, cliqsym)
  syms = getCliqAllVarIds(cliq)
  subfg = buildSubgraphFromLabels!(fgl,syms)

  # add upward messages to subgraph
  msgs = getCliqParentMsgDown(treel,cliq)
  addMsgFactors!(subfg, msgs)

  # predictBelief
  # stuff = treeProductUp(fgl, treel, cliqsym, varsym)
  # plotKDE(manikde!(stuff[1], getManifolds(fgl, varsym)))
  cllbl = cliq.attributes["label"]
  return plotLocalProduct(subfg, varsym, title="Tree Dwn $(cllbl) | ", levels=levels)
end


"""
    $SIGNATURES

Overlay plot all upward messages from cliques.
"""
function plotCliqUpMsgs(fg::G,
                        tree::AbstractBayesTree,
                        sym::Symbol;
                        show::Bool=true,
                        dims::Vector{Int}=Int[],
                        levels::Int=1,
                        c=nothing,
                        title="up msgs on $(sym)"    ) where G <: AbstractDFG
  #
  # get all msgs
  allmsgs = getTreeCliqUpMsgsAll(tree)
  # stack messages by variable
  sckmsgs = stackCliqUpMsgsByVariable(tree, allmsgs)

  if !haskey(sckmsgs, sym)
    @warn "plotCliqUpMsgs -- tree does not have up messages for $sym."
    return nothing
  end

  # prepend current estimate too
  Xs = getBelief(fg, sym)
  # vectorize beliefs
  beliefs = BallTreeDensity[Xs;]
  lbls = String["curr,-1";]
  for msg in sckmsgs[sym]
    push!(beliefs, msg.belief)
    push!(lbls, "$(msg.cliqId),$(msg.depth)")
  end

  # ignoring legend and color information
  cc = c === nothing ? getColorsByLength(length(beliefs)) : c

  # plot and return
  plotKDE(beliefs, levels=levels, c=cc, title=title, legend=lbls)
end


"""
    $SIGNATURES

Development function to plot the same variable from both factor graphs together.
"""
function plotPairVariables(dfg1::G1,
                           dfg2::G2,
                           sym::Symbol;
                           dims=nothing,
                           levels::Int=3,
                           title::String="") where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  X1 = getBelief(dfg1, sym)
  X2 = getBelief(dfg2, sym)

  plotKDE([X1;X2], c=["black";"red"], legend=["dfg1";"dfg2"], dims=dims, levels=levels, title=title*" $sym")
end


function plotPairPose2(dfg1::G1,
                           dfg2::G2,
                           sym::Symbol;
                           dims=nothing,
                           levels::Int=3,
                           title::String="") where {G1 <: AbstractDFG, G2 <: AbstractDFG}
  #
  X1 = getKDE(dfg1, sym)
  X2 = getKDE(dfg2, sym)

  plotPose(Pose2(), [X1;X2], title*" $sym", c=["black";"red"], levels=levels)
  # , legend=["dfg1";"dfg2"]
end

"""
    $SIGNATURES

Plot new proposal (convolution) for factor x of variable y given all other -- i.e. y|X

Notes
- plot new result on `towards` as first color.
- Plot other variables in factor on 2nd to m colors.
"""
function plotVariableGivenFactor(dfg::G,
                                 fct::Symbol,
                                 towards::Symbol;
                                 levels::Int=2,
                                 dims=nothing  ) where G <: AbstractDFG
  #
  pts = approxConv(dfg,fct,towards)
  mani = getManifolds(dfg, towards)
  res = manikde!(pts,mani)

  lie = ls(dfg, fct)
  setdiff!(lie, [towards])

  otr = map(x->getKDE(dfg,x),lie)
  lbls = string.([towards;lie])

  pl = plotKDE([res;otr],dims=dims,levels=levels,legend=lbls)

  return pl
end


"""
    $SIGNATURES

Plot the downward messages currently stored in a clique.
"""
function plotCliqDownMsgs(tree::AbstractBayesTree,
                          frnt::Symbol;
                          show::Bool=true,
                          levels::Int=2,
                          dims=nothing,
                          existing=nothing  )
  #
  cliq = getClique(tree,frnt)
  msgs = IIF.getMessageBuffer(cliq).downRx.belief

  PL = []

  for (key, beldim) in msgs
    npl = plotKDE(manikde!(beldim.val, beldim.softtype), levels=levels, title="dwn msg $key", dims=dims)
    existing === nothing ? nothing : union!(npl.layers, existing.layers)
    push!(PL, npl)
  end

  existing === nothing ? nothing : push!(PL, existing)
  pl = vstack(PL...)

  folderpath = "/tmp/caesar/random/"
  filepath = folderpath*"downmsgs_cliq$(cliq.index).pdf"
  Base.mkpath(folderpath)
  pl |> PDF(filepath, 20cm, length(PL)*12cm)
  @async run(`evince $filepath`)
  pl
end
