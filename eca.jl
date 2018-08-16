
function initializePopY(fobj, best, N, D, a, b)
    L = 0.5*(b - a)
    r = 0.01

    # println(size(L))
    y0 = best.y

    a2 = y0 .- r*L
    b2 = y0 .+ r*L

    Y = a2' .+ (b2- a2)' .* rand(N, D)


    y = Y[1,:]
    child = mh.generateChild(y, fobj(y))
    individual = typeof(child)

    # population array
    population = Array{individual, 1}([])

    # first individual
    push!(population, child)
    push!(population, mh.generateChild(y0, best.f))

    for i in 3:N
        y = Y[i,:]
        child = mh.generateChild(y, fobj(y))
        push!(population, child)
    end

    return population
end

function eca(f::Function, best, D::Int, bounds)
    # general parameters
    searchType = :minimize
    η_max = 2.0
    K = 7
    α = 10
    N = div(K*D, 2)
    max_evals = 100D
    #############################

    x = best.x
    fobj(y) = f(x, y)

    a, b = bounds[1,:], bounds[2,:]

    # initialize population
    Population = initializePopY(fobj, best, N, D, a, b)

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
            y = Population[i].x

            # generate U masses
            U = mh.getU(Population, K, I, i, N)
            
            # generate center of mass
            c, u_worst, u_best = mh.center(U, searchType)

            # stepsize
            η = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            
            # current-to-center/bin
            h = y + η * (c - u)
            
            h = mh.correct(h, a, b, true)

            sol = mh.generateChild(h, fobj(h))
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


        t += 1
    end

    return best.x, best.f, nevals
end

