
@info "RoMEPlotting.jl is loading tools based on ImageMagick.jl"


using .ImageMagick

export makeImage

"""
$SIGNATURES

Use PNG from Cairo and Fontconfig to convert a Gadfly plot into a Julia Images style matrix.

Example
```julia
pl = Gadlfy.plot(y=randn(10), Geom.line)
makeImage(pl)
```

See also: [`toFormat`](@ref)
"""
function makeImage(pl::Gadfly.Plot)
  buf = IOBuffer()
  pl |> PNG(buf)
  plbytes = take!(buf)
  ImageMagick.readblob(plbytes)
end

# deprecate before RoME v0.18.0
import Base: convert
@deprecate convert(::Type{<:AbstractArray{<:Color, 2}}, pl::Gadfly.Plot) makeImage(pl)


#