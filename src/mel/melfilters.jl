mutable struct MelSpec{T}
    fs      :: Int   # 采样率，一般 16kHz
    winlen  :: Int   # 分帧参数-帧长
    stride  :: Int   # 分帧参数-帧移
    nbanks  :: Int   # 梅尔滤波器个数
    nffts   :: Int   # 傅立叶变换点数
    maxfreq :: Int   # 频域最大有效频率下标
    alpha   :: T     # 预加重系数，默认0.97
    eps     :: T     # 防止下溢系数
    winfunc :: Function       # 窗计算函数，如 hanning
    window  :: AbstractArray  # 窗函数数值表示
    melbank :: AbstractArray  # mel滤波器矩阵

    function MelSpec{T}(;
        fs      :: Int = 16000,
        winlen  :: Int = 256,
        stride  :: Int = winlen>>1,
        nffts   :: Int = winlen,
        nbanks  :: Int = 32,
        alpha   :: AbstractFloat = 0.97,
        eps     :: AbstractFloat = 1e-6,
        winfunc :: Function = hanning) where T <: AbstractFloat

        @assert nffts>=winlen "nffts ≥ winlen, but got nffts=$nffts, winlen=$winlen"
        maxfreq = nffts>>1
        window  = winfunc(winlen, T)
        melbank = filterbanks(nbanks, nffts, fs, dtype=T)
        new{T}(fs, winlen, stride, nbanks, nffts, maxfreq, alpha, eps, winfunc, window, melbank)
    end
end


function MelSpec(;
    fs     :: Int = 16000,
    winlen :: Int = 256,
    stride :: Int = winlen>>1,
    nffts  :: Int = winlen,
    nbanks :: Int = 32,
    alpha  :: AbstractFloat = 0.97,
    eps    :: AbstractFloat = 1e-6,
    winfunc:: Function = hanning,
    dtype  :: DataType = Float32)

    return MelSpec{dtype}(
        fs      = fs,
        winlen  = winlen,
        stride  = stride,
        nffts   = nffts,
        nbanks  = nbanks,
        alpha   = alpha,
        eps     = eps,
        winfunc = winfunc)
end


function filterbanks(nbanks::Int, nffts::Int, fs::Int; dtype::DataType=Float32)
    MAX   = nffts>>1;                    # 正频率部分的下标最大值
    freq  = (0:(MAX-1))/MAX * fs/2;      # 下标映射到频率
    Fmel  = 2595*log10.(1 .+ freq/700);  # 频率映射到梅尔频率
    dFmel = Fmel[MAX]/(nbanks+1);        # 将 Mel 带平分
    bank  = zeros(dtype, nbanks, MAX);   # 滤波器的频域权重系数
    cFmel = 0.0;                         # 每个Mel频带的中心Mel频率
    for n = 1:nbanks
        cFmel = cFmel + dFmel
        for m = 1:MAX
            if (cFmel - dFmel) ≤ Fmel[m] ≤ (cFmel + dFmel)
                bank[n,m] = 1.0 - abs( Fmel[m] - cFmel )/dFmel
            end
        end
    end
    return bank
end


function filterwav(data, α)
    n = length(data)
    @. data[1:n-1] = data[2:n] - α * data[1:n-1]
    data[n : n] = data[n-1 : n-1]
    return data
end


function (filter::MelSpec{T})(wav::S, func::Union{Function,Nothing}=log) where {T,S <: AbstractArray}
    W = filter.melbank
    α = filter.alpha
    B = filter.eps
    F = filter.maxfreq
    winlen = filter.winlen
    stride = filter.stride
    window = filter.window
    nffts  = filter.nffts
    # 滤波 & 功率谱
    d = filterwav(wav, α)
    X = spectrum(d, window, winlen, stride, nffts, dtype=T)

    if func !== nothing
        return func.(W * X .+ B)
    else
        return W * X
    end
end

function Base.show(io::IO, ::MIME"text/plain", f::MelSpec{T}) where T
    winlen = trunc(f.winlen / f.fs * 1000, digits=3)
    stride = trunc(f.stride / f.fs * 1000, digits=3)
    println(io, "════════════════════════════")
    println(io, " MelSpec{dtype=$(T)}")
    println(io, "════════════════════════════")
    println(io, " winfunc = $(f.winfunc)")
    println(io, "  winlen = $(f.winlen) ($winlen ms)")
    println(io, "  stride = $(f.stride) ($stride ms)")
    println(io, "   nffts = $(f.nffts)")
    println(io, " maxfreq = $(f.maxfreq)")
    println(io, "  nbanks = $(f.nbanks)")
    println(io, "   alpha = $(f.alpha)")
    println(io, "     eps = $(f.eps)")
    println(io, "      fs = $(f.fs) Hz")
    println(io, "════════════════════════════")
end
