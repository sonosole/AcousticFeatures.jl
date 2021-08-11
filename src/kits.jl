export filterwav

function filterwav(data, α)
    n = length(data)
    @. data[1:n-1] = data[2:n] - α * data[1:n-1]
    data[n] = data[n-1]
    return data
end


export splitwav

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
