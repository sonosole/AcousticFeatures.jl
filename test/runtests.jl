using Plots: plot, heatmap, title!, gui
using ChirpSignal: chirp
using AcousticFeatures
using Test

include("fbank-dtype.jl")
include("fbank-wins.jl")
include("powspec.jl")

@test true
