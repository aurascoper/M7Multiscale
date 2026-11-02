# M7Multiscale/visualization/Timeline4D/Timeline4D.jl
# Workspace 5: 4D Spatiotemporal & Comparison World (3-clock timeline, Difference Mode, Conservation).

function extract_timeline4d_summary(frames::AbstractVector{M7SpatialFrame})
    !isempty(frames) || throw(ArgumentError("frames must be non-empty"))
    cons = audit_m7_spatiotemporal_conservation(frames)
    return (
        workspace = :spatiotemporal_world,
        frame_count = length(frames),
        physical_time_span_s = (first(frames).clocks.physical_time_s, last(frames).clocks.physical_time_s),
        biological_time_span_s = (first(frames).clocks.biological_time_s, last(frames).clocks.biological_time_s),
        sobolev_step_span = (first(frames).clocks.sobolev_relaxation_step, last(frames).clocks.sobolev_relaxation_step),
        conservation_audit = cons
    )
end
