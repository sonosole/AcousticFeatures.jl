module AcousticFeatures

include("./kits.jl")

include("./winfuns.jl")
using .WindowFunctions
export barthann
export bartlett
export blackman
export blackmanharris
export bohman
export flattop
export hamming
export hanning
export nuttall
export parzen
export rectangular
export triangular


include("./melfilters.jl")


end # module
