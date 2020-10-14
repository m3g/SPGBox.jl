#
# Computes the norm of the projected gradient
#
function pr_gradnorm(x,g,l,u)
  gnorm = 0.
  for i in 1:length(x)
    z = max(l[i], min(u[i], x[i]-g[i])) - x[i]
    gnorm = max(gnorm, abs(z))
  end
  return gnorm
end
