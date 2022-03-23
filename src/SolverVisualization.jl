

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
                          title::String="Local product ($solveKey), ",
                          xmin=nothing, xmax=nothing, ymin=nothing, ymax=nothing  )
  #
  @warn "not showing partial constraints, but included in the product"
  arr = Array{BallTreeDensity,1}()
  lbls = String[]
  push!(arr, getBelief(getVariable(fgl, lbl), solveKey))
  push!(lbls, "curr")
  pl = nothing
  pp, parr, lb, sinfd = IIF.localProduct(fgl, lbl, N=N, solveKey=solveKey)

  partials = []

  # another sanity check
  xmin !== nothing && xmax !== nothing && xmin == xmax ? error("xmin must be less than xmax") : nothing
  ymin !== nothing && ymax !== nothing && ymin == ymax ? error("ymin must be less than ymax") : nothing

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
        newbel_ = getBelief(getVariable(fgl, lbl),solveKey)
        newbel = newbel_ isa ManifoldKernelDensity ? newbel_.belief : newbel_
        pl = plotKDE( [proddim;newbel;vals], dims=[1;], levels=levels, c=colors, title=string("Local product, dim=$(dimn), ",lbl))
        push!(PL, pl)
      end
      Gadfly.vstack(PL...)
  end

  @show length(parr), length(partials)
  if length(parr) > 0 && length(partials) == 0
    pl = plotDirectProducts()
  elseif length(parr) == 0 && length(partials) > 0
    pl = plotPartialProducts()
  else
    return error("plotLocalProduct not built for lengths parr, partials = $(length(parr)), $(length(partials)) yet.")
  end

  # # set coordinates accordingly
  # co = Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  # pl.coord = co

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
  push!(arr, getBelief(getVariable(fgl, lbl)))
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
        pll = plotKDE([proddim;getBelief(getVariable(fgl, lbl));vals], dims=[1;], levels=levels, c=colors, title=string("Local product, dim=$(dimn), ",lbl))
        push!(PL, pl)
      else
        plc = AMP.plotKDECircular([proddim; marginal(getBelief(getVariable(fgl, lbl)),[2]);vals], c=colors, scale=scale)
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
plotLocalProduct(fgl::AbstractDFG, lbl::AbstractString; N::Int=100, dims::Vector{Int}=Int[] ) = plotLocalProduct(fgl, Symbol(lbl), N=N, dims=dims)



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
  X1 = getBelief(dfg1, sym)
  X2 = getBelief(dfg2, sym)

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
  mani = getManifold(dfg, towards)
  res = manikde!(mani, pts)

  lie = ls(dfg, fct)
  setdiff!(lie, [towards])

  otr = map(x->getBelief(dfg,x),lie)
  lbls = string.([towards;lie])

  pl = plotKDE([res;otr],dims=dims,levels=levels,legend=lbls)

  return pl
end


#