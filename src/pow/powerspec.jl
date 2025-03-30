"""
Operator for computing the power spectrum
# Constructor
    PowerSpec(;
              fs       :: Int = 16000,          # sampling frequency
              fmin     :: Real = 0,             # min frequency
              fmax     :: Real = fs/2,          # max frequency
              winlen   :: Int = 256,            # window length
              dilation :: Int = 1,              # down-sample ratio
              stride   :: Int = winlen/2,       # distance between two adjacent frames
              nffts    :: Int = winlen,         # the bins of FFT transform
              donorm   :: Bool = false,         # whether normalize the spectrum
              winfunc  :: Function = hanning,   # the windowing function
              type     :: DataType = Vector{Float32})

where [`fmin`, `fmax`] is the frequency region of interest, `type` is the input's type (Matrix or Vector). 
If the `donorm` is true, then the elements of window would be divided by the sqrt of window energy.
"""
mutable struct PowerSpec{T}
    fs       :: Int     # 采样率如16kHz
    fminidx  :: Int     # 最小频率下标
    fmaxidx  :: Int     # 最大频率下标
    winlen   :: Int     # 分帧参数-帧长
    dilation :: Int     # 膨胀系数/降采样因子
    stride   :: Int     # 分帧参数-帧移
    nffts    :: Int     # 傅立叶变换点数
    window   :: T       # 窗函数数值表示
    function PowerSpec(;
        fs       :: Int = 16000,
        fmin     :: Real = 0,
        fmax     :: Real = fs>>1,
        winlen   :: Int = 256,
        dilation :: Int = 1,
        stride   :: Int = winlen>>1,
        nffts    :: Int = winlen,
        winfunc  :: Function = hanning,
        donorm   :: Bool = false,
        type     :: DataType = Vector{Float32})

        @assert nffts ≥ winlen > 0 begin
            "nffts ≥ winlen, but got nffts=$nffts, winlen=$winlen"
        end
        @assert 0 ≤ fmin < fmax ≤ fs/2 begin
            "0 ≤ fmin < fmax ≤ fs/2, but got 0 ≤ $fmin < $fmax ≤ $(fs/2)"
        end
        # (idx-1)/(nffts-1) = (f - 0)/(fs - 0)
        minidx = floor(Int, (nffts-1) * fmin / fs + 1)
        maxidx = floor(Int, (nffts-1) * fmax / fs + 1)
        window = winfunc(winlen, type)
        if donorm
            coeff = sqrt(sum(window .^ 2))
            window .*= inv(coeff)
        end
        new{type}(fs, minidx, maxidx, winlen, dilation, stride, nffts, window)
    end
end


# pretty showing
function Base.show(io::IO, ::MIME"text/plain", P::PowerSpec{T}) where T
    fs = P.fs
    winlen = trunc(P.winlen / fs * 1000, digits=3)
    stride = trunc(P.stride / fs * 1000, digits=3)
    fmin = floor(Int, (P.fminidx - 1)/(P.nffts - 1) * fs)
    fmax =  ceil(Int, (P.fmaxidx - 1)/(P.nffts - 1) * fs)
    println(io, "════════════════════════════")
    println(io, " PowerSpec{$T}")
    println(io, "════════════════════════════")
    println(io, "dilation = $(P.dilation)")
    println(io, "  winlen = $(P.winlen) ($winlen ms)")
    println(io, "  stride = $(P.stride) ($stride ms)")
    println(io, "   nffts = $(P.nffts)")
    println(io, " minfreq = $(fmin)")
    println(io, " maxfreq = $(fmax)")
    println(io, "      fs = $(fs) Hz")
    println(io, "════════════════════════════")
end


# deal with single-channel input
function (P::PowerSpec{Vector{D}})(data::Vector{D}) where D <: AbstractFloat
    frames = split4fft(data, P.window, P.dilation, P.stride, P.nffts)
    L = P.fminidx            # Lower bound
    U = P.fmaxidx            # Upper bound
    T = size(frames, 2)      # Time steps
    C = fft(frames, 1)       # 时域到频域 𝐑ⁿ → 𝐂ⁿ,按帧所在维度计算
    X = abs2.(C[L:U,1:T])    # 功率谱,提取有用部分
    return X
end


# deal with single-channel input
function (P::PowerSpec{Vector{D}})(data::Matrix{D}) where D <: AbstractFloat
    @assert size(data, 2) == 1 "PowerSpec{Vector{$D}} only supports one channel input"
    frames = split4fft(vec(data), P.window, P.dilation, P.stride, P.nffts)
    L = P.fminidx            # Lower bound
    U = P.fmaxidx            # Upper bound
    T = size(frames, 2)      # Time steps
    C = fft(frames, 1)       # 时域到频域 𝐑ⁿ → 𝐂ⁿ,按帧所在维度计算
    X = abs2.(C[L:U,1:T])    # 功率谱,提取有用部分
    return X
end


# deal with multi-channel input
function (P::PowerSpec{Matrix{D}})(data::Matrix{D}) where D <: AbstractFloat
    frames = split4fft(data, P.window, P.dilation, P.stride, P.nffts)
    L = P.fminidx             # Lower bound
    U = P.fmaxidx             # Upper bound
    T = size(frames, 3)       # Time steps
    C = fft(frames, 2)        # 时域到频域 𝐑ⁿ → 𝐂ⁿ,按帧所在维度计算
    X = abs2.(C[:, L:U, 1:T]) # 功率谱,提取感兴趣部分
    return X
end

