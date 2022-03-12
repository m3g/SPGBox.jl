#
# Computes the norm of the projected gradient
# One method is defined for each possible scenario 
# of bounds (all, only lower, only upper, none). This is
# needed such that the lower and upper bounds vectors can 
# be of type Nothing if those bounds are not defined
#
function pr_gradnorm(g, x, lower, upper)
    if isnothing(lower) && isnothing(upper)
        return pr_gradnorm_no_bounds(g)
    elseif isnothing(lower)
        pr_gradnorm_upper(g, x, upper)
    elseif isnothing(upper)
        return pr_gradnorm_lower(g, x, lower)
    else
        return pr_gradnorm_both_bounds(g, x, lower, upper)
    end
end
# Both bounds
function pr_gradnorm_both_bounds(g, x, lower, upper)
    gnorm = zero(eltype(x))
    for i in eachindex(x)
        z = max(lower[i], min(upper[i], x[i] - g[i])) - x[i]
        gnorm = max(gnorm, abs(z))
    end
    return gnorm
end
# Only lower
function pr_gradnorm_lower(g, x, lower)
    gnorm = zero(eltype(x))
    for i in eachindex(x)
        z = max(lower[i], x[i] - g[i]) - x[i]
        gnorm = max(gnorm, abs(z))
    end
    return gnorm
end
# Only upper
function pr_gradnorm_upper(g, x, upper)
    gnorm = zero(eltype(x))
    for i in eachindex(x)
        z = min(upper[i], x[i] - g[i]) - x[i]
        gnorm = max(gnorm, abs(z))
    end
    return gnorm
end
# No bounds
function pr_gradnorm_no_bounds(g)
    gnorm = zero(eltype(g))
    for i in eachindex(g)
        gnorm = max(gnorm, abs(g[i]))
    end
    return gnorm
end
