


"""
    $SIGNATURES

Calculate the ellipse from covariance matrix and return as lambda function.

Notes
- https://cookierobotics.com/007/

Related

plotCovEllipseLayer
"""
function covEllipseParameterized(distr::MvNormal; meanOffset::Bool=true)
  # println(round.(cov(distr), digits=3) )
  a = cov(distr)[1,1]
  b = cov(distr)[1,2] + cov(distr)[2,1]
  b *= 0.5
  c = cov(distr)[2,2]

  λ1 = 0.5*(a+c) + sqrt(0.25*(a-c)^2 + b^2)
  λ2 = 0.5*(a+c) - sqrt(0.25*(a-c)^2 + b^2)
  θ = if b ≈ 0 && a >= c
    0
  elseif b ≈ 0 && a < c
    pi/2
  else
    atan(λ1-a, b)
  end
  sλ1 = sqrt(λ1)
  sλ2 = sqrt(λ2)
  sθ =  sin(θ)
  cθ =  cos(θ)
  ox = meanOffset ? distr.μ[1] : 0
  oy = meanOffset ? distr.μ[2] : 0
  (t) -> [sλ1*cθ*cos(t)-sλ2*sθ*sin(t)+ox; sλ1*sθ*cos(t)+sλ2*cθ*sin(t)+oy]
end

covEllipseParameterized(pts::Array{Float64,2}; meanOffset::Bool=true) = covEllipseParameterized( fit(MvNormal, pts), meanOffset=meanOffset )
covEllipseParameterized(X::Union{<:BallTreeDensity,<:ManifoldKernelDensity}; meanOffset::Bool=true) = covEllipseParameterized( getPoints(X), meanOffset=meanOffset )
covEllipseParameterized(dfg::AbstractDFG, sym::Symbol; meanOffset::Bool=true, solveKey::Symbol=:default) = covEllipseParameterized( getBelief(dfg, sym, solveKey), meanOffset=meanOffset, solveKey=solveKey )


"""
    $SIGNATURES

Plotting tool to draw Gadfly layers of ellipses of 2D covariance fitted to the belief of factor graph variable nonparametric points.

Related

covEllipseParameterized, plotSLAM2DPosesLandms
"""
function plotCovEllipseLayer( dfg::AbstractDFG,
                              vsym::Symbol;
                              solveKey::Symbol=:default,
                              drawPoints::Bool=true,
                              ellipseColor::AbstractString="gray30",
                              pointsColor::AbstractString="gray30",
                              drawEllipse::Bool=true  )
  #
  PL = []

  if !(drawEllipse || drawPoints)
    return PL
  end

  # points to work from
  bel = getBelief(dfg, vsym, solveKey)
  pp__ = getPoints(bel)
  pp_ = if bel.manifold isa SpecialEuclidean
    # assume Pose2
    (x->x.parts[1]).(pp__)
  else
    # assume Point2
    pp__
  end

  @cast pp[i,j] := pp_[j][i]

  if drawEllipse
    # get ellipse function
    eX = covEllipseParameterized(pp[1:2,:], meanOffset=false)
    vEl = eX.(0:0.02:2pi)
    el = [(x->x[1]).(vEl) (x->x[2]).(vEl)]
    # add suggested PPE mean offset
    el[:,1] .+= getVariablePPE(dfg, vsym, solveKey).suggested[1]
    el[:,2] .+= getVariablePPE(dfg, vsym, solveKey).suggested[2]

    # add the ellipse layers
    plelX2 = Gadfly.layer(x=el[:,1], y=el[:,2], Geom.path, Theme(default_color=parse(Colorant, ellipseColor)))
    push!(PL, plelX2[1])
  end

  # add the points layer if needed
  if drawPoints
    plelX1 = Gadfly.layer(x=pp[1,:],
                          y=pp[2,:],
                          Geom.point,
                          Theme(default_color=parse(Colorant, pointsColor), 
                          point_size=1pt))
    #
    push!(PL, plelX1[1])
  end

  return PL
end


