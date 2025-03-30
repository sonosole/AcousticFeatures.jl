mutable struct OnlineMelSpec{T}
    fs      :: Int   # 采样率，一般 16kHz
    winlen  :: Int   # 分帧参数-帧长
    stride  :: Int   # 帧移参数
    nbanks  :: Int   # 梅尔滤波器个数
    nffts   :: Int   # 傅立叶变换点数
    maxfreq :: Int   # 频域最大有效频率下标
    keptone :: T     # 上一帧需要保留的采样点
    alpha   :: T     # 预加重系数，默认 0.97
    eps     :: T     # 防止下溢系数
    winfunc :: Function       # 窗计算函数，如 hanning
    window  :: AbstractArray  # 窗函数数值表示
    melbank :: AbstractArray  # mel滤波器矩阵

    function OnlineMelSpec{T}(;
        fs      :: Int = 16000,
        winlen  :: Int = 256,
        stride  :: Int = winlen,
        nffts   :: Int = winlen,
        nbanks  :: Int = 32,
        alpha   :: AbstractFloat = 0.97,
        eps     :: AbstractFloat = 1e-6,
        winfunc :: Function = hanning) where T <: AbstractFloat

        @assert nffts>=winlen "nffts ≥ winlen, but got nffts=$nffts, winlen=$winlen"
        maxfreq = nffts>>1
        window  = winfunc(winlen, T)
        melbank = AcousticFeatures.filterbanks(nbanks, nffts, fs, dtype=T)
        new{T}(fs, winlen, stride, nbanks, nffts, maxfreq, 0.0, alpha, eps, winfunc, window, melbank)
    end
end


function OnlineMelSpec(;
    fs     :: Int = 16000,
    winlen :: Int = 256,
    stride :: Int = winlen,
    nffts  :: Int = winlen,
    nbanks :: Int = 32,
    alpha  :: AbstractFloat = 0.97,
    eps    :: AbstractFloat = 1e-6,
    winfunc:: Function = hanning,
    dtype  :: DataType = Float32)

    return OnlineMelSpec{dtype}(
        fs      = fs,
        winlen  = winlen,
        stride  = stride,
        nffts   = nffts,
        nbanks  = nbanks,
        alpha   = alpha,
        eps     = eps,
        winfunc = winfunc)
end


function (filter::OnlineMelSpec{T})(wav::S, func::Union{Function,Nothing}=log) where {T,S <: AbstractArray}
    W = filter.melbank
    α = filter.alpha
    B = filter.eps
    F = filter.maxfreq

    winlen = filter.winlen
    stride = filter.stride
    window = filter.window
    nffts  = filter.nffts

    # 滤波
    keptone = wav[stride]
    for i = winlen:-1:2
        wav[i] = wav[i] - α * wav[i-1]
    end;wav[1] = wav[1] - α * filter.keptone
    filter.keptone = keptone

    # 加窗 & fft & 功率谱
    d = wav .* window
    if winlen < nffts
        x = zeros(eltype(d), nffts, 1)
        copyto!(x, d)
        C = fft(x)
        X = abs2.(C[1:F])
    else
        C = fft(d)
        X = abs2.(C[1:F])
    end

    if func !== nothing
        return func.(W * X .+ B)
    else
        return W * X
    end
end

function Base.show(io::IO, ::MIME"text/plain", f::OnlineMelSpec{T}) where T
    winlen = trunc(f.winlen / f.fs * 1000, digits=3)
    stride = trunc(f.stride / f.fs * 1000, digits=3)
    println(io, "═════════════════════════════")
    println(io, " OnlineMelSpec{dtype=$(T)}")
    println(io, "═════════════════════════════")
    println(io, " winfunc = $(f.winfunc)")
    println(io, "  winlen = $(f.winlen) ($winlen ms)")
    println(io, "  stride = $(f.stride) ($stride ms)")
    println(io, "   nffts = $(f.nffts)")
    println(io, " maxfreq = $(f.maxfreq)")
    println(io, "  nbanks = $(f.nbanks)")
    println(io, "   alpha = $(f.alpha)")
    println(io, "     eps = $(f.eps)")
    println(io, "      fs = $(f.fs) Hz")
    println(io, "═════════════════════════════")
end
