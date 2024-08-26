module SPGBox

using TestItems: @testitem
using Compat: @compat

export spgbox!
export spgbox
export SPGBoxResult
@compat public VAux

include("./SPGBoxResult.jl")
include("./VAux.jl")
include("./auxiliary_functions.jl")
include("./spgbox_main.jl")

end
