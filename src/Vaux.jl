#
# Structure to store auxiliary arrays
#
struct VAux
  g :: Vector{Float64}
  xn :: Vector{Float64}
  gn :: Vector{Float64}
  fprev :: Vector{Float64}
end
VAux(n,m) = VAux( Vector{Float64}(undef,n), 
                  Vector{Float64}(undef,n), 
                  Vector{Float64}(undef,n), 
                  Vector{Float64}(undef,m) )
# The default value for m
VAux(n) = VAux(n,10)
