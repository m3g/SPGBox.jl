module SPGBox

using TestItems: @testitem
using Compat: @compat

export spgbox!
export spgbox
export SPGBoxResult
export SPGBoxOpt
@compat public VAux

include("./SPGBoxResult.jl")
include("./VAux.jl")
include("./auxiliary_functions.jl")
include("./spgbox_main.jl")
include("./SPGBoxOpt.jl")

end
