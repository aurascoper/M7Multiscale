# M7Multiscale/contracts/TimeFrame.jl
# Three-clock multirate temporal contract (physical_time_s, biological_time_s, sobolev_relaxation_step).

using .M7VisualizationBridge:
    M7MultirateClock,
    validate_m7_temporal_interpolation,
    interpolate_m7_spatial_frames
