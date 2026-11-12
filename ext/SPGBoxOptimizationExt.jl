module SPGBoxOptimizationExt

import OptimizationBase
import SciMLBase
using SciMLBase: ReturnCode, OptimizationStats
using SPGBox: SPGBoxOpt, spgbox!

SciMLBase.allowscallback(::SPGBoxOpt) = true
SciMLBase.has_init(::SPGBoxOpt) = true
SciMLBase.allowsbounds(::SPGBoxOpt) = true
SciMLBase.requiresgradient(::SPGBoxOpt) = true

function SciMLBase.__solve(cache::OptimizationBase.OptimizationCache{<:SPGBoxOpt})
    opt = cache.opt
    x = copy(cache.u0)
    lower = isnothing(cache.lb) ? nothing : convert(typeof(x), cache.lb)
    upper = isnothing(cache.ub) ? nothing : convert(typeof(x), cache.ub)

    maxiters = OptimizationBase._check_and_convert_maxiters(cache.solver_args.maxiters)
    maxtime = cache.solver_args.maxtime
    limits = if isnothing(maxiters)
        (;)
    else
        (; nitmax=maxiters, nfevalmax=get(opt.kwargs, :nfevalmax, typemax(Int)))
    end

    ngeval = Ref(0)
    f = u -> cache.f(u, cache.p)
    g! = (g, u) -> begin
        ngeval[] += 1
        cache.f.grad(g, u)
    end

    # With MaxSense, OptimizationBase passes -f to the solver; report the user's objective.
    sgn = cache.sense === SciMLBase.MaxSense ? -1 : 1

    t0 = time()
    timed_out = Ref(false)
    callback = R -> begin
        objective = sgn * R.f
        state = OptimizationBase.OptimizationState(;
            iter=R.nit, u=R.x, objective, original=R
        )
        stop = cache.callback(state, objective)
        if !isnothing(maxtime) && time() - t0 > maxtime
            timed_out[] = true
            stop = true
        end
        return stop
    end

    res = spgbox!(f, g!, x; lower, upper, callback, merge(opt.kwargs, limits)...)

    retcode = if timed_out[]
        ReturnCode.MaxTime
    elseif res.return_from_callback
        ReturnCode.Terminated
    elseif res.ierr == 0
        ReturnCode.Success
    else
        ReturnCode.MaxIters
    end
    stats = OptimizationStats(;
        iterations=res.nit, time=time() - t0, fevals=res.nfeval, gevals=ngeval[]
    )
    return SciMLBase.build_solution(
        cache, opt, res.x, res.f; original=res, retcode, stats
    )
end

end
