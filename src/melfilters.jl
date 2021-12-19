using FFTW:fft
export MelFilters


mutable struct MelFilters
    winlen::Int   # 分帧参数-帧长
    stride::Int   # 分帧参数-帧移
    nbanks::Int   # 梅尔滤波器个数
    nffts ::Int   # 傅立叶变换点数
    alpha::Real   # 预加重系数，默认0.97
    fs::Int       # 采样率，一般16kHz
    maxfreq::Int  # 频域最大有效频率下标
    eps::Real     # 防止下溢系数
    winfunc::Function       # 窗计算函数，如 hanning
    window ::AbstractArray  # 窗函数数值表示
    melbank::AbstractArray  # mel滤波器矩阵

    function MelFilters(;
        winlen::Int = 256,
        stride::Int = winlen>>1,
        nffts ::Int = winlen,
        nbanks::Int = 32,
        alpha::Real = 0.97,
        fs::Int     = 16000,
        eps::Real   = 1e-6,
        winfunc::Function = hanning)
        @assert nffts>=winlen "nffts ≥ winlen, but got nffts=$nffts, winlen=$winlen"
        maxfreq = nffts>>1
        window  = winfunc(winlen)
        melbank = filterbanks(nbanks, nffts, fs)
        new(winlen,stride, nbanks,nffts,alpha,fs, maxfreq,eps, winfunc, window,melbank)
    end
end


function filterbanks(nbanks::Int, nffts::Int, fs::Int)
    MAX   = nffts>>1;                    # 正频率部分的下标最大值
    freq  = (0:(MAX-1))/MAX * fs/2;      # 下标映射到频率
    Fmel  = 2595*log10.(1 .+ freq/700);  # 频率映射到梅尔频率
    dFmel = Fmel[MAX]/(nbanks+1);        # 将Mel带平分
    bank  = zeros(nbanks, MAX);          # 滤波器的频域权重系数
    cFmel = 0.0;                         # 每个Mel频带的中心Mel频率
    for n = 1:nbanks
        cFmel = cFmel + dFmel
        for m = 1:MAX
            if ( Fmel[m] >= cFmel-dFmel ) && ( Fmel[m] <= cFmel+dFmel )
                bank[n,m] = 1.0 - abs( Fmel[m] - cFmel )/dFmel
            end
        end
    end
    return bank
end


function (filter::MelFilters)(wav, func::Union{Function,Nothing}=log)
    W = filter.melbank
    α = filter.alpha
    B = filter.eps
    F = filter.maxfreq
    winlen = filter.winlen
    stride = filter.stride
    window = filter.window
    nffts  = filter.nffts

    d = filterwav(wav, α)                            # 滤波
    X = spectrum(d, window, winlen, stride, nffts)   # 功率谱

    if func !== nothing
        return func.(W * X .+ B)
    else
        return W * X
    end
end


function Base.show(io::IO, f::MelFilters)
    winlen = trunc(f.winlen / f.fs * 1000, digits=3)
    stride = trunc(f.stride / f.fs * 1000, digits=3)
    println(io, "———————————————————————")
    println(io, " winfunc = $(f.winfunc)")
    println(io, "  winlen = $(f.winlen) ($winlen ms)")
    println(io, "  stride = $(f.stride) ($stride ms)")
    println(io, "   nffts = $(f.nffts)")
    println(io, " maxfreq = $(f.maxfreq)")
    println(io, "  nbanks = $(f.nbanks)")
    println(io, "   alpha = $(f.alpha)")
    println(io, "     eps = $(f.eps)")
    println(io, "      fs = $(f.fs) Hz")
    println(io, "———————————————————————")
end
