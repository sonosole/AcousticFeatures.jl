using BenchmarkTools

fbank32 = MelSpec(winlen=256, stride=128, dtype=Float32);
fbank64 = MelSpec(winlen=256, stride=128, dtype=Float64);

x32 = randn(Float32, 16000*7, 1);
x64 = randn(Float64, 16000*7, 1);

@btime fbank32(x32);
@btime fbank64(x64);
