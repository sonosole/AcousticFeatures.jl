fsrate  = 16000
fbank32 = MelSpec(fs=fsrate, winlen=256, stride=128, dtype=Float32);
fbank64 = MelSpec(fs=fsrate, winlen=256, stride=128, dtype=Float64);

x32 = randn(Float32, fsrate*7, 1);
x64 = randn(Float64, fsrate*7, 1);

@time fbank32(x32);
@time fbank64(x64);
