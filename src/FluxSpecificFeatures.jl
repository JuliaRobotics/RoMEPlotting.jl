# features specifically related to Flux.jl

@eval PlotTypesPose2 = Union{Type{Pose2Pose2}, Type{Pose2Point2BearingRange}, Type{Pose2Point2Range}, Type{Pose2Point2Bearing}, Type{FluxModelsPose2Pose2}}

@eval ExtendedPose2Pose2Types = Union{Pose2Pose2, FluxModelsPose2Pose2}


#
