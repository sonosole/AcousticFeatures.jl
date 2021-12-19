export filterwav
export splitwav
export spectrum


function filterwav(data, α)
    n = length(data)
    @. data[1:n-1] = data[2:n] - α * data[1:n-1]
    data[n] = data[n-1]
    return data
end



function splitwav(data, win, winlen::Int, stride::Int)
    nframes  = div(length(data)-winlen, stride) + 1
    firstIds = (0:(nframes-1)) .* stride .+ 1   # 帧起始下标
    lasstIds = firstIds .+ (winlen - 1)         # 帧结束下标
    frames   = zeros(eltype(data), winlen, nframes)
    for i = 1:nframes
        frames[:,i] = data[firstIds[i]:lasstIds[i]] .* win
    end
    return frames, nframes
end

function spectrum(data, win, winlen::Int, stride::Int, nffts::Int)
    # @assert nffts>=winlen
    T = div(length(data)-winlen, stride) + 1   # nframs
    F = nffts>>1                               # max effective freq bin index
    firstIds = (0:(T-1)) .* stride .+ 1        # 帧起始下标
    lasstIds = firstIds .+ (winlen - 1)        # 帧结束下标
    frames   = zeros(eltype(data), nffts, T)   # 如果 nffts>winlen 则零填充帧
    for i = 1:T
        frames[1:winlen,i] = data[firstIds[i]:lasstIds[i]] .* win
    end
    C = fft(frames, 1)       # 时域到频域 𝐑ⁿ ↣ 𝐂ⁿ ,按列计算
    X = abs2.(C[1:F,1:T])    # 功率谱,提取有用部分
    return X
end
