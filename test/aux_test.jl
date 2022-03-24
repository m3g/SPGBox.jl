using NLPModels
using CUTEst

function cutest2spg(probname)
    prob = CUTEstModel(probname)
    
    f(x) = obj(prob, x)
    function g!(g, x)
        g .= grad(prob, x)
    end
    x0 = copy(prob.meta.x0)
    lower = copy(prob.meta.lvar)
    upper = copy(prob.meta.uvar)
    return f, g!, x0, lower, upper, prob
end
