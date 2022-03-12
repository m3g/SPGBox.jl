#
# Structure to store auxiliary arrays
#
struct VAux{T,F}
  g :: Vector{T}
  xn :: Vector{T}
  gn :: Vector{T}
  fprev :: Vector{F}
end
VAux(::Type{T},::Type{F},n,m) where {T,F} = VAux{T,F}(zeros(T,n),zeros(T,n),zeros(T,n),zeros(F,m))
VAux(n,m) = VAux(Float64,Float64,n,m)
VAux(::Type{T},::Type{F},n) where {T,F} = VAux(T,F,n,10)
VAux(n) = VAux(Float64,Float64,n)
