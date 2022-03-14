#
# Structure to store auxiliary arrays
#
struct VAux{T,F}
    g::T
    xn::T
    gn::T
    fprev::Vector{F}
end
function VAux(x,fx;m=10)
    fprev = typeof(fx)[]
    for i in 1:m
        push!(fprev,zero(fx))
    end
    return VAux{typeof(x),typeof(fx)}(
        zero(similar(x)),
        zero(similar(x)),
        zero(similar(x)),
        fprev
    )
end