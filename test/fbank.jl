using Plots:plot, heatmap
using WAV:wavread

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
triangular]

plts = []
data, Fs = wavread("AMonoAudioFile.wav");
for winfn in fnlist
    fbank = MelSpec(
          fs = floor(Int,Fs),
       alpha = 0.5,
      winlen = 1024,
      stride = 128,
      nbanks = 128,
      nffts  = 1024,
     epsilon = 1e-6,
     winfunc = winfn)

    feat = fbank(copy(data));
    temp = heatmap(feat, legend=nothing, framestyle=:box, ticks=nothing, titlefont=(20, "times"))
    title!("$winfn")
    push!(plts, temp)
end

plot(plts...,layout=(4,3),size=(1600,1200))
