import Metaheuristics
const mh = Metaheuristics


function apply_nested!(P, f, D_upper, D_lower, lower_bounds)
	return 0
end

function hbo(F::Function, f::Function, D_upper, D_lower, upper_bounds, lower_bounds)
    # general parameters
    searchType = :minimize
    η_max = 2.0
    K = 7
    α = 10
    D = D_upper + D_lower
    N = 2K*D
    bounds = [upper_bounds lower_bounds]
    τ = 10
    max_evals = 1000D
    #############################
    # auto conf
    T = round(Int, max_evals / N)
    τ_ratio = div(T, τ)
    #############################

    # penalized function
    φ(z) = F(z[1:D_upper], z[D_upper+1:end]) + α * f(z[1:D_upper], z[D_upper+1:end])

    a, b = bounds[1,:], bounds[2,:]

    # initialize population
    Population = mh.initializePop(φ, N, D, a, b, :uniform)

    # current evaluations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solution
    best = mh.getBest(Population, searchType)
  
    # start search
    while !stop
        I = randperm(N)

        for i in 1:N

            # current
            x = Population[i].x

            # generate U masses
            U = mh.getU(Population, K, I, i, N)
            
            # generate center of mass
            c, u_worst, u_best = mh.center(U, searchType)

            # stepsize
            η = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            
            # current-to-center/bin
            h = x + η * (c - u)
            
            h = mh.correct(h, a, b, true)

            sol = mh.generateChild(h, φ(h))
            nevals += 1


            # replace worst element
            if mh.Selection(Population[i], sol, searchType)
                Population[mh.getWorstInd(Population, searchType)] = sol

                if mh.Selection(best, sol, searchType)
                    best = sol
                end
            end
            
            stop = nevals >= max_evals
            if stop
                break
            end
        end

        if t != 0 && t % τ_ratio == 0
            apply_nested!(Population, f, D_upper, D_lower, lower_bounds)
        end

        t += 1
    end

    println("+----------------------------------+")
    println("|          HBO results             |")
    println("+----------------------------------+")
    mh.printResults(best, Population, t, nevals)
    println("+----------------------------------+")
    return best.x, best.f
end