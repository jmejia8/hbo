import Metaheuristics
const mh = Metaheuristics

include("structures.jl")
include("tools.jl")
include("operators.jl")
include("eca.jl")


function fitnessToMass(fitness::Array{Float64,2}, searchType::Symbol)
    m_ul = minimum(fitness[:,1])
    m_ll = minimum(fitness[:,2])
    
    m_ul < 0 && (fitness[:,1] = 2abs(m_ul) + fitness[:,1])
    m_ll < 0 && (fitness[:,1] = 2abs(m_ll) + fitness[:,2])
  

    if searchType == :minimize
        fitness[:,1] = 2maximum(fitness[:,1]) - fitness[:,1]
        fitness[:,2] = 2maximum(fitness[:,2]) - fitness[:,2]
    end

    return fitness
end


function getMass(U::Array, V::Array, α, β, searchType::Symbol)
    n = length(U)

    fitness = zeros(Float64, n, 2)
    
    for i = 1:n
        fitness[i, 1] = U[i].F + α*U[i].f
        fitness[i, 2] = β*V[i].F + V[i].f
    end

    return fitnessToMass(fitness, searchType)
end


function center(U::Array, V::Array, mass::Array{Float64,2})
    d_ul = length(U[1].x)
    d_ll = length(V[1].y)

    c_ul = zeros(Float64, d_ul)
    c_ll = zeros(Float64, d_ll)
    
    for i = 1:length(U)
        c_ul += mass[i,1] * U[i].x
        c_ll += mass[i,2] * V[i].y
    end

    return c_ul / sum(mass[:,1]), c_ll / sum(mass[:,2])
end

function center(U::Array, V::Array, α, β, searchType::Symbol)
    n = length(U)

    mass = getMass(U, V, α, β, searchType)

    a, b = center(U, V, mass)
    return a, b, getWorstInd(U, searchType), indmin(mass[:,2])
end


function apply_nested!(P, F, f, D_ul, D_ll, bounds)
	nevals_ul = 0
    nevals_ll = 0
    
    for i = 1:length(P)
		y, fy, nevals = eca(f, P[i], D_ll, bounds)

        if fy >= P[i].f
            continue
        end


        P[i].y = y
        oldf = P[i].f
        P[i].f = fy
        oldF = P[i].F
        P[i].F = F(P[i].x, P[i].y)

        nevals_ul += 1
        nevals_ll += nevals
        @printf("fnew = %e \t fold = %e \t | \t Fnew = %e \t Fold = %e  \n", P[i].f, oldf, P[i].F, oldF)
	end

    return nevals_ul, nevals_ll

end

function getWorstVals(Population::Array, searchType=:minimize)
    Fworst, fworst = Population[1].F, Population[1].f

    for i = 2:length(Population)
        Fworst <  Population[i].F && (Fworst =  Population[i].F)
        fworst <  Population[i].f && (fworst =  Population[i].f)
    end

    return Fworst, fworst
end

α = β = 0.1

function hbo(F::Function, f::Function, D_ul, D_ll, bounds_ul::Matrix, bounds_ll::Matrix; showResults = true)
    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 

    # general parameters
    searchType = :minimize
    D = D_ul + D_ll
    κ = 3
    η_max = 2.0
    N = D < 5 ? κ*5 : κ*D
    max_evals_ul = 3000D
    max_evals_ll = 3000D
    #############################
    # auto conf
    # T = round(Int, max_evals / N)
    #############################

    # initialize population
    Population = initializePop(F, f, N, bounds_ul, bounds_ll)


    Fworst, fworst = getWorstVals(Population, searchType)

    # println(Fworst)
    # println(fworst)
 
    ########################################################
    # the sequences
    α0 = β0 = 0.1
    global α = α0
    global β = β0
    # β(t::Int, β0::Float64 = 10.0) = (1/β0) * ( 1.0 - (t / max_evals_ll)^3 )
    # α(t::Int, α0::Float64 = 10.0) = β(t, α0)
    ########################################################

    # current evaluations
    nevals_ul = N
    nevals_ll = N

    # stop condition
    stop = false

    # current generation
    t = 0

    # best solution
    best = getBest(Population, searchType)

    convergence = [best]

  
    # start search
    while !stop
        I_ul = randperm(N)
        I_ll = randperm(N)

        # @printf("t = %d \t a = %e \t b = %e \n", t, α, β)
        global α = α0 * (1.0 - (nevals_ul / max_evals_ul)^3)
        global β = β0 * (1.0 - (nevals_ul / max_evals_ll)^9)
        for i in 1:N

            # current
            x = Population[i].x
            y = Population[i].y

            # generate U masses
            U = mh.getU(Population, κ, I_ul, i, N)
            V = mh.getU(Population, κ, I_ll, i, N)
            
            # generate center of mass
            c_ul, c_ll, u_worst, v_worst = center(U, V, α, β, searchType)

            # stepsize
            η_ul = η_max * rand()
            η_ll = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            v = V[v_worst].y
            
            # current-to-center
            p = x + η_ul * (c_ul - u)
            q = y + η_ll * (c_ll - v)
            
            p = correct(p, bounds_ul)
            q = correct(q, bounds_ll)

            sol = generateChild(p, q, F(p, q), f(p, q))
            nevals_ul += 1
            nevals_ll += 1


            # replace worst element
            if is_better_mass(sol, Population[i], searchType)
                Population[getWorstInd(Population, searchType)] = sol

                if is_better_mass(sol, best, searchType)
                    best = sol
                    push!(convergence, best)
                    # @printf("F = %e \t f = %e  \n", best.F, best.f)
                end
            end
            
            stop = nevals_ul >= max_evals_ul || nevals_ll >= max_evals_ll || (abs(best.f) < 1e-4 && abs(best.F) < 1e-4)
            if stop
                break
            end
        end

        # if t != 0 && t % τ_ratio == 0
        #     n1, n2 = apply_nested!(Population, F, f, D_ul, D_ll, bounds_ll)
        #     nevals_ul += n1
        #     nevals_ll += n2  
        # end

        t += 1
    end



    if showResults
        println("+----------------------------------+")
        println("|          HBO results             |")
        println("+----------------------------------+")
        printResults(best, Population, t, nevals_ul, nevals_ll)
        println("+----------------------------------+")
    end

    return best.x, best.y, best, nevals_ul, Population, convergence
end