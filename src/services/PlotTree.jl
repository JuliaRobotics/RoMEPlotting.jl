


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
    npl = plotKDE(manikde!(beldim.varType, beldim.val), levels=levels, title="dwn msg $key", dims=dims)
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



#