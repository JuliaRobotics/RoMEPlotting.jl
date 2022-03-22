


plotBelief(mkd::ManifoldKernelDensity,w...;kw...) = plotKDE(mkd.belief, w...;kw...)
plotBelief(arr::AbstractVector{<:ManifoldKernelDensity},w...;kw...) = plotKDE((x->x.belief).(arr), w...;kw...)


"""
    $(SIGNATURES)

A peneric KDE plotting function that allows marginals of higher dimensional beliefs and various keyword options.

Example for `Position2`:
-----------------------
```julia

p = manikde!(Position2, [randn(2) for _ in 1:100])
q = manikde!(Position2, [randn(2).+[5;0] for _ in 1:100])

plotBelief(p)
plotBelief(p, dims=[1;2], levels=3)
plotBelief(p, dims=[1])

plotBelief([p;q])
plotBelief([p;q], dims=[1;2], levels=3)
plotBelief([p;q], dims=[1])
```

Example for `Pose2`:
--------------------
```julia
# TODO
```
"""
function plotBelief(fgl::AbstractDFG,
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
  p = getBelief(getVariable(fgl,sym), solveKey)
  # mmarg = length(marg) > 0 ? marg : collect(1:Ndim(p))
  # mp = marginal(p,mmarg)
  bel = p isa BallTreeDensity ? p : p.belief
  plotKDE(bel, levels=levels, dims=dims, title=string(sym, "  ", title), fill=fill, layers=layers, c=c, overlay=overlay )
end

function plotBelief(fgl::AbstractDFG,
                    syms::Vector{Symbol};
                    solveKey::Symbol=:default,
                    addt::Union{<:AbstractVector{<:BallTreeDensity},AbstractVector{<:ManifoldKernelDensity}}=BallTreeDensity[],
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
    if p isa BallTreeDensity
      push!(MP, p)
    else
      push!(MP, p.belief)
    end
    push!(LEG, string(sym))
  end
  for p in addt
    # mp = marginal(p,mmarg)
    push!(MP, p)
    push!(LEG, "add")
  end
  plotKDE(MP, c=c, levels=levels, dims=dims, legend=LEG, title=title, layers=layers, overlay=overlay)
end


#