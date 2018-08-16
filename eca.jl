import Metaheuristics
const mh = Metaheuristics

function eca(f::Function, best, D::Int, bounds)
    # general parameters
    searchType = :minimize
    η_max = 2.0
    K = 7
    α = 10
    N = div(K*D, 2)
    bounds = [upper_bounds lower_bounds]
    max_evals = 100D
    #############################

    fobj(y) = f(best.x, y)

    a, b = bounds[1,:], bounds[2,:]

    # initialize population
    Population = initializePopY(fobj, N, D, a, b, :uniform)

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

    return best
end