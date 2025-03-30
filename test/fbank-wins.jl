# different window functions
fnlist = [
          barthann,
          bartlett,
          blackman,
          blackmanharris,
          bohman,
          flattop,
          hamming,
          hanning,
          nuttall,
          parzen,
          rectangular,
          triangular];

# create a chirp signal
fl = 50.0
fh = 6000.0
Fs = 16000
T  = 2.0
data = chirp(T, Fs, fl, fh);

# compare via heatmap
plts = []
for winfn ∈ fnlist
    fbank = MelSpec(
          fs = Fs,
         eps = 1e-5,
       alpha = 0.97,
      winlen = 1024,
      stride = 128,
      nbanks = 64,
      nffts  = 1024,
     winfunc = winfn)

    feat = fbank(copy(data));
    temp = heatmap(feat, legend=nothing, framestyle=:box, ticks=nothing, titlefont=(13, "times"))
    title!("$winfn")
    push!(plts, temp)
end

plot(plts..., layout=(3,4), size=(1800,1000))
gui()
sleep(5)
