export splitwav
export spectrum


function splitwav(data::Vector{D}, win::Vector{D}, winlen::Int, stride::Int) where D <: AbstractFloat
    nframes  = div(length(data)-winlen, stride) + 1
    firstIds = (0:(nframes-1)) .* stride .+ 1   # 帧起始下标
    lasstIds = firstIds .+ (winlen - 1)         # 帧结束下标
    frames   = zeros(eltype(data), winlen, nframes)
    for i = 1:nframes
        s = firstIds[i]
        f = lasstIds[i]
        frames[:,i] = data[s:f] .* win
    end
    return frames, nframes
end


function spectrum(data, win, winlen::Int, stride::Int, nffts::Int; dtype=Float32)
    # @assert nffts>=winlen
    T = div(length(data)-winlen, stride) + 1   # nframs
    F = nffts>>1                               # max effective freq bin index
    firstIds = (0:(T-1)) .* stride .+ 1        # 帧起始下标
    lasstIds = firstIds .+ (winlen - 1)        # 帧结束下标
    frames   = zeros(dtype, nffts, T)          # 如果 nffts>winlen 则零填充帧
    for i = 1:T
        s = firstIds[i]
        f = lasstIds[i]
        frames[1:winlen,i] = data[s:f] .* win
    end
    C = fft(frames, 1)       # 时域到频域 𝐑ⁿ ↣ 𝐂ⁿ ,按列计算
    X = abs2.(C[1:F,1:T])    # 功率谱,提取有用部分
    return X
end



"""
    split4fft(X::Vector{D}, W::Vector{D}, dilation::Int, stride::Int, nffts::Int) -> Y::Array{D,2}

Split single-channel 1-D signal `X` into frame buffers `Y`. The following conditions are satisfied:
    nsamples = length(`X`)
    nffts, nframes = size(`Y`)
"""
function split4fft(X::Vector{D}, W::Vector{D}, dilation::Int, stride::Int, nffts::Int) where {D <: AbstractFloat}
    winlen = length(W)
    kernel = (winlen - 1) * dilation + 1        # 等效窗长
    T = div(length(X) - kernel, stride) + 1    # 分帧帧数
    starts = 1 .+ (0:(T-1)) .* stride         # 帧起始下标
    finals = starts .+ (kernel - 1)          # 帧结束下标
    frames = Array{D}(undef, nffts, T)      # 如nffts>winlen则零填充

    for t = 1:T
        s = starts[t]
        f = finals[t]
        frames[1:winlen, t:t] .= X[s:dilation:f] .* W
    end
    if nffts > winlen
        frames[winlen+1 : nffts, 1:T] .= zero(D)
    end
    return frames
end


"""
    split4fft(X::Array{D,2}, W::Array{D,2}, dilation::Int, stride::Int, nffts::Int) -> Y::Array{D,3}

Split multi-channel 1-D signal `X` into frame buffers `Y`. The following conditions are satisfied:
    nsamples, nchannels = size(`X`)
    nchannels, nffts, nframes = size(`Y`)
"""
function split4fft(X::Matrix{D}, W::Matrix{D}, dilation::Int, stride::Int, nffts::Int) where {D <: AbstractFloat}
    L, C = size(X)  # length & channels
    winlen = length(W)
    kernel = (winlen - 1) * dilation + 1        # 等效窗长
    T = (L - kernel) ÷ stride + 1               # 分帧帧数
    starts = 1 .+ (0:(T-1)) .* stride           # 帧起始下标
    finals = starts .+ (kernel - 1)             # 帧结束下标
    frames = Array{D,3}(undef, C, nffts, T)     # 如nffts>winlen则零填充

    Threads.@threads for c = 1:C
        for t = 1:T
            s = starts[t]
            f = finals[t]
            frames[c, 1:winlen, t:t] .= X[s:dilation:f, c] .* W
        end
    end
    if nffts > winlen
        frames[1:C, winlen+1 : nffts, 1:T] .= zero(D)
    end
    return frames
end
