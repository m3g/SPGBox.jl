#
# Structure to store auxiliary arrays
#
struct VAux{T}
  g :: Vector{T}
  xn :: Vector{T}
  gn :: Vector{T}
  fprev :: Vector{T}
end
VAux(n,m) = VAux(Float64,n,m)

VAux(::Type{T},n,m) where T = VAux( zeros(T,n), #BigFloat needs to be allocated first
                                    zeros(T,n),
                                    zeros(T,n),
                                    zeros(T,m))
# The default value for m
VAux(n) = VAux(n,10)
