include("tools.jl")

function apply_nested!(P, f, D_upper, D_lower, lower_bounds)
	return 0
end

function center(U::Array, mass::Vector{Float64})
    d = length(U[1].x)

    c = zeros(Float64, d)
    
    for i = 1:length(mass)
        c += mass[i] * U[i].x
    end

    return c / sum(mass)
end

function center(U::Array, V::Array, searchType::Symbol)
    n = length(U)

    mass = getMass(U, searchType)

    return center(U, mass), getWorstInd(U, searchType), getBestInd(U, searchType)
end

function fitnessToMass(fitness::Vector{Float64}, searchType::Symbol)
    m = minimum(fitness)
    
    if m < 0
        fitness = 2abs(m) + fitness
    end

    if searchType == :minimize
        fitness = 2maximum(fitness) - fitness
    end

    return fitness
end

function getMass(U, searchType::Symbol)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)
    
    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness, searchType)
end


function hbo(F::Function, f::Function, D_upper, D_lower, bounds_ul::Matrix, bounds_ll::Matrix)
    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 

    # general parameters
    searchType = :minimize
    η_max = 2.0
    κ = 7
    α = 10
    D = D_upper + D_lower
    N = 2κ*D
    τ = 10
    max_evals = 1000D
    #############################
    # auto conf
    T = round(Int, max_evals / N)
    τ_ratio = div(T, τ)
    #############################

    # initialize population
    Population = initializePop(F, f, N, D, bounds_ul, bounds_ll)

    # current evaluations
    nevals = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solution
    best = getBest(Population, searchType)
  
    # start search
    while !stop
        I_ul = randperm(N)
        I_ll = randperm(N)

        for i in 1:N

            # current
            x = Population[i].x
            y = Population[i].y

            # generate U masses
            U = mh.getU(Population, κ, I_ul, i, N)
            V = mh.getU(Population, κ, I_ll, i, N)
            
            # generate center of mass
            c_ul, c_ll, u_worst, v_worst = center(U, V, searchType)

            # stepsize
            η_ul = η_max * rand()
            η_ll = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            v = U[u_worst].y
            
            # current-to-center
            h_ul = x + η_ul * (c_ul - u)
            h_ll = y + η_ll * (c_ll - v)
            
            h_ul = mh.correct(h, a, b, true)
            h_ll = mh.correct(h, a, b, true)

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