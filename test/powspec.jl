# create chirp signals
fl = 1000.0  # start frequency
fh = 6000.0  # final frequency
Fs = 16000   # sampling frequency
T  = 5.0     # duration of signal
data1 = chirp(T, Fs, fl, fh, method="linear");
data2 = chirp(T, Fs, fl, fh, method="logarithmic");
datas = [data1 data2];   # 2-channel input

# single channel powspec
singlepow = PowerSpec(
            fs = Fs,
            fmin = 0,
            fmax = 8000,
            winlen = 256,
            dilation = 1,
            stride = 128,
            nffts = 256,
            winfunc = bohman,
            type = Vector{Float64}) # Vector represents single channel input

# multi-channel powspec
multipows = PowerSpec(
            fs = Fs,
            fmin = 0,
            fmax = 8000,
            winlen = 256,
            dilation = 1,
            stride = 128,
            nffts = 256,
            winfunc = bohman,
            type = Matrix{Float64}) # Matrix represents mutiple channels input

# get power spectrums
y1 = singlepow(data1);
y2 = singlepow(data2);
ys = multipows(datas);

@test all(ys[1,:,:] .== y1)
@test all(ys[2,:,:] .== y2)

plot(heatmap(ys[1,:,:]), heatmap(ys[2,:,:]), layout=(2,1))
gui();sleep(5)
