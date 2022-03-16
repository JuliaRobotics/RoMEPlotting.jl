# plotPose functions


"""
    $SIGNATURES

Plot pose belief as contour information on visually sensible manifolds.

Example:

```julia
fg = generateCanonicalFG_Hexagonal()
solveTree!(fg);
plotPose(fg, :x6)
```

Related

[`plotSLAM2D`](@ref), [`plotSLAM2DPoses`](@ref), [`plotKDE`](@ref), `plotKDECircular`
"""
function plotPose(::IIF.InstanceType{Pose2},
                  pp::AbstractVector{<:BallTreeDensity},
                  title="plotPose2";
                  levels=3,
                  c=nothing,
                  legend=nothing,
                  axis=nothing,
                  scale::Real=0.2,
                  overlay=nothing,
                  hdl=[]  )
  #
  # ops = buildHybridManifoldCallbacks(pt.manifolds)
  # @show ran = getKDERange(p, addop=ops[1], diffop=ops[2])
  ran = axis === nothing ? getKDERange(pp) : axis

  p1 = plotKDE(pp, dims=[1;2], levels=levels, c=c, axis=ran )
  # p2 = plotKDE(bels, dims=[3], c=c)

  cc = c === nothing ? getColorsByLength(length(pp)) : c

  GG = BallTreeDensity[]
  for ppc in pp
    gg = marginal(ppc,[3])
    # gg = (x)->pc(reshape([x], :,1))[1]
    push!(GG, gg)
  end
  # p2 = AMP.plotCircBeliefs(GG, c=cc)
  p2 = AMP.plotKDECircular(GG, scale=scale, c=cc, legend=legend, title=title)

  # deal with overlay

  push!(hdl, p1)
  push!(hdl, p2)

  Gadfly.hstack(p1,p2)
end

plotPose(pt::IIF.InstanceType{Pose2}, bds::AbstractVector{<:ManifoldKernelDensity}, w...;kw...) = plotPose(pt, (s->s.belief).(bds), w...; kw...)

function plotPose(pt::IIF.InstanceType{Pose2},
                  pp::Union{<:BallTreeDensity,<:ManifoldKernelDensity},
                  title="plotPose2";
                  levels=3,
                  c=nothing,
                  axis=nothing,
                  scale::Real=0.2,
                  overlay=nothing,
                  hdl=[]  )
  #
  plotPose(pt, [pp;],title,levels=levels,c=c,axis=axis,scale=scale, overlay=overlay, hdl=hdl )
end


function plotPose(::DynPose2,
                  bels::Union{<:AbstractVector{<:BallTreeDensity},<:AbstractVector{<:ManifoldKernelDensity}},
                  title;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  hdl=[],
                  scale::Real=0.2 )
  #
  p1 = plotKDE(bels, dims=[1;2], levels=levels, c=c, title=title)
  p2 = plotKDE(bels, dims=[3], c=c)
  p3 = plotKDE(bels, dims=[4;5], levels=levels, c=c)

  push!(hdl, p1)
  push!(hdl, p2)
  push!(hdl, p3)

  Gadfly.vstack(p1,p2,p3)
end

# import RoMEPlotting: plotPose

function plotPose(::Pose3,
                  bels::Union{<:AbstractVector{<:BallTreeDensity},<:AbstractVector{<:ManifoldKernelDensity}},
                  title;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  hdl=[],
                  scale::Real=0.2  )
  #
  @show title
  p1 = plotKDE(bels, dims=[1;2], levels=levels, c=c, title=title)
  p2 = plotKDE(bels, dims=[3], c=c)
  p3 = plotKDE(bels, dims=[4;5], levels=levels, c=c)
  p4 = plotKDE(bels, dims=[6], c=c)

  push!(hdl, p1)
  push!(hdl, p2)
  push!(hdl, p3)
  push!(hdl, p4)

  Gadfly.vstack(p1,p2,p3,p4)
end

"""
    $(SIGNATURES)

Example: pl = plotPose(fg, [:x1; :x2; :x3])
"""
function plotPose(fgl::AbstractDFG,
                  syms::Vector{Symbol};
                  solveKey::Symbol=:default,
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  scale::Real=0.2,
                  show::Bool=false,
                  filepath::AbstractString="/tmp/tempposeplot.svg",
                  app::AbstractString="eog",
                  hdl=[]  )
  #
  typ = getSolverData(getVariable(fgl, syms[1]), solveKey).variableType
  pt = string(string.(syms)...)
  getvertsgg = (sym) -> getBelief(getVariable(fgl, sym), solveKey)
  pl = plotPose(typ, getvertsgg.(syms), pt, levels=levels, c=c, axis=axis, scale=scale, hdl=hdl )

  if length(filepath) > 0
    ext = split(filepath, '.')[end]
    cmd = getfield(Gadfly,Symbol(uppercase(ext)))

    h = 1*7Gadfly.cm
    if typ == DynPose2
        h *= 1.5
    end
    Gadfly.draw(cmd(filepath,15Gadfly.cm,h),pl)


    @async !show ? nothing : run(`$app $filepath`)
  end
  return pl
end

function plotPose(fgl::AbstractDFG,
                  sym::Symbol;
                  levels::Int=5,
                  c=nothing,
                  axis=nothing,
                  scale::Real=0.2,
                  show::Bool=false,
                  filepath::AbstractString="/tmp/tempposeplot.svg",
                  app::AbstractString="eog",
                  hdl=[]  )
  #
  plotPose(fgl, [sym;], levels=levels, axis=axis, show=show, filepath=filepath, app=app, hdl=hdl )
end
