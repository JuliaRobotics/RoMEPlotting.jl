# features specifically related to Flux.jl
@info "RoMEPlotting is adding Flux related functionality."

@eval PlotTypesPose2 = Union{Type{Pose2Pose2}, Type{Pose2Point2BearingRange}, Type{Pose2Point2Range}, Type{Pose2Point2Bearing}, Type{MixtureFluxPose2Pose2}}

@eval ExtendedPose2Pose2Types = Union{Pose2Pose2, MixtureFluxPose2Pose2}

# rerun since a few functions need to be reloaded with Flux
include("PlotFactorsReload.jl")

#
