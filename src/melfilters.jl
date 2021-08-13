using FFTW:fft
export MelFilters


mutable struct MelFilters
    winlen::Int   # 分帧参数-帧长
    stride::Int   # 分帧参数-帧移
    nbanks::Int   # 梅尔滤波器个数
    nffts::Int    # 傅立叶变换点数
    alpha::Real   # 预加重系数，默认0.97
    fs::Int       # 采样率，一般16kHz
    maxfreq::Int  # 频域最大有效频率下标
    epsilon::Real # 防止下溢系数
    winfunc::Function
    windown::AbstractArray
    melbank::AbstractArray

    function MelFilters(;
        winlen::Int       = 256,
        stride::Int       = winlen>>1,
        nffts::Int        = winlen,
        nbanks::Int       = 32,
        alpha::Real       = 0.97,
        fs::Int           = 16000,
        epsilon::Real     = 1e-6,
        winfunc::Function = hanning)

        maxfreq = nffts>>1
        windown = winfunc(winlen)
        melbank = filterbanks(nbanks, nffts, fs)
        new(winlen,stride, nbanks,nffts,alpha,fs, maxfreq,epsilon, winfunc, windown,melbank)
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
    B = filter.epsilon
    F = filter.maxfreq
    winlen  = filter.winlen
    stride  = filter.stride
    windown = filter.windown

    data = filterwav(wav, α)                             # 滤波
    frames, T = splitwav(data, windown, winlen, stride)  # 分帧
    C = fft(frames, 1)       # 时域到频域 𝐑ⁿ ↣ 𝐂ⁿ ,按列计算
    X = abs2.(C[1:F,1:T])    # 功率谱,提取有用部分
    if func !== nothing
        return func.(W * X .+ B)
    else
        return W * X
    end
end


function Base.show(io::IO, f::MelFilters)
    println(io, "———————————————————————")
    println(io, "  winlen = $(f.winfunc)")
    println(io, "  winlen = $(f.winlen)")
    println(io, "  stride = $(f.stride)")
    println(io, "  nbanks = $(f.nbanks)")
    println(io, "   nffts = $(f.nffts)")
    println(io, "   alpha = $(f.alpha)")
    println(io, "      fs = $(f.fs)")
    println(io, " maxfreq = $(f.maxfreq)")
    println(io, " epsilon = $(f.epsilon)")
    println(io, "———————————————————————")
end
