
@info "RoMEPlotting is loading tools based on ImageMagick.jl"

import Base: convert

export convert

"""
    $SIGNATURES

Use PNG from Cairo and Fontconfig to convert a Gadfly plot into a Julia Images style matrix.

Example
```julia
pl = Gadlfy.plot(y=randn(10), Geom.line)
convert(Matrix{RGB}, pl)
```
"""
function convert(::Type{<:AbstractArray{<:Color, 2}}, pl::Gadfly.Plot)
  buf = IOBuffer()
  pl |> PNG(buf)
  plbytes = take!(buf)
  ImageMagick.readblob(plbytes)
end