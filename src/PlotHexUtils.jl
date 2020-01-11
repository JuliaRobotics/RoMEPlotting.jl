
export layerLine, plotHex!, plotHex, plotHex_6, plotBeehive_6

# draw a line on plot
function layerLine(start::Vector{Float64}, stop::Vector{Float64}, color=colorant"green")
  layer(x=[start[1];stop[1]], y=[start[2];stop[2]], Geom.path, Theme(default_color=color))[1]
end

# draw a default hex shape on plot
function plotHex!(offsetxy::Vector{Float64}=zeros(2);
                  len::Float64=10.0,
                  PLT = [])
    #
    push!(PLT, layerLine([offsetxy[1]+0.0;offsetxy[2]+0.0], [offsetxy[1]+len;offsetxy[2]+0.0]))
    push!(PLT, layerLine([offsetxy[1]+len;offsetxy[2]+0.0], [offsetxy[1]+len*3/2;offsetxy[2]+len*sqrt(3)/2]))
    push!(PLT, layerLine([offsetxy[1]+len+5.0;offsetxy[2]+len*sqrt(3)/2], [offsetxy[1]+len;offsetxy[2]+len*sqrt(3)]))
    push!(PLT, layerLine([offsetxy[1]+0;offsetxy[2]+len*sqrt(3)], [offsetxy[1]+10;offsetxy[2]+len*sqrt(3)]))
    push!(PLT, layerLine([offsetxy[1]+0;offsetxy[2]+len*sqrt(3)], [offsetxy[1]-len/2; offsetxy[2]+len*sqrt(3)/2]))
    push!(PLT, layerLine([offsetxy[1]-len/2; offsetxy[2]+len*sqrt(3)/2], [offsetxy[1]+0.0; offsetxy[2]+0]))

    PLT
end

# draw a default hex shape on plot
function plotHex(offsetxy::Vector{Float64}=zeros(2);
                 len::Float64=10.0)
    #
    Gadfly.plot(plotHex!(offsetxy, len=len)...)
end

# draw a default ring of 6 hexes on plot
function plotHex_6(;PLT=[], len::Float64=10.0)
    plotHex!(zeros(2),PLT=PLT)
    plotHex!([3/2*len;-sqrt(3)/2*len],PLT=PLT)
    plotHex!([3/2*len;-sqrt(3)*3/2*len],PLT=PLT)
    plotHex!([0.0;-sqrt(3)*2*len],PLT=PLT)
    plotHex!([-3/2*len;-sqrt(3)*3/2*len],PLT=PLT)
    plotHex!([-3/2*len;-sqrt(3)/2*len],PLT=PLT)

    Gadfly.plot(PLT...)
end

"""
    $SIGNATURES

Plot first ring (ground truth) of beehive example.
"""
function plotBeehive_6(fgl::G; from::Int64=0, to::Int64=99999999, meanmax::Symbol=:max) where G <: AbstractDFG
  pl = plotHex_6()
  pl2 = drawPosesLandms(fgl, from=from, to=to, meanmax=meanmax)

  union!(pl2.layers, pl.layers)

  return pl2
end
