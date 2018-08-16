###################################################
#      Solutions and population functions
#          for Matrix representation
###################################################
function correctSol(z::Vector{Float64}, bounds::Matrix)
    # Correct solution

    for i = 1:length(z)
        if !( bounds[1,i] <= z[i] <= bounds[2,i] )
            z[i] = bounds[1,i] + (bounds[2,i] - bounds[1,i])*rand()
        end
    end
    
    return z
end

correct(z::Vector{Float64}, bounds::Matrix) = correctSol(z, bounds)


function initializePop(N::Int, D::Int, a::Vector{Float64}, b::Vector{Float64}, initType::Symbol=:uniform)
    # a, b should be D × 1

    if initType == :cheb
        chebPts(x, a, b) = 0.5*(a + b) + 0.5*(b-a)*cos.( x )
        X = zeros(N, D)
        for j in 1:D
            X[:, j] = chebPts(2π*rand(N), a[j], b[j])
        end

        return X
    end
    return a'  .+ (b - a)' .* rand(N, D)
end

function initializePop(F::Function, f::Function, N::Int, bounds_ul::Matrix, bounds_ll::Matrix)
    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 
   
    a_ul, b_ul = bounds_ul[1,:], bounds_ul[2,:]
    a_ll, b_ll = bounds_ll[1,:], bounds_ll[2,:]
   
    X = initializePop(N, D_ul, a_ul, b_ul)
    Y = initializePop(N, D_ll, a_ll, b_ll)
    
    # infers datatype
    x = X[1,:]
    y = Y[1,:]

    child = generateChild(x, y, F(x, y), f(x, y))
    individual = typeof(child)

    # population array
    population = Array{individual, 1}([])

    # first individual
    push!(population, child)

    for i in 2:N
        x = X[i,:]
        y = Y[i,:]

        child = generateChild(x, y, F(x, y), f(x, y))
        push!(population, child)
    end

    return population
end

function getfValues(P::Array)
    F = zeros(length(P))
    f = zeros(length(P))

    for i = 1:length(P)
        F[i] = P[i].F
        f[i] = P[i].f
    end

    return F, f
end

function getPositions(P::Array)
    N, D_ul, D_ll = length(P), length(P[1].x), length(P[1].y)

    X = zeros(N, D_ul)
    Y = zeros(N, D_ll)

    for i = 1:N
        X[i,:] = P[i].x
        Y[i,:] = P[i].y
    end

    return X, Y
end


# COP functions
function violationsSum(g::Vector, h::Vector)
    sum_g = 0
    sum_h = 0

    for i = 1:length(g)
        if g[i] > 0
            sum_g += g[i]  end
    end

    for i = 1:length(h)
        if h[i] != 0.0
            sum_h += abs(h[i])  end
    end

    return sum_g + sum_h
end

# for Deb rules
function countViolations(g::Vector, h::Vector)
    sum_g = 0
    sum_h = 0

    for i = 1:length(g)
        if g[i] > 0
            sum_g += 1  end
    end

    for i = 1:length(h)
        if h[i] != 0.0
            sum_h += 1  end
    end

    return sum_g + sum_h
end

function printResults(best::xf_indiv, P, t, nevals_ul, nevals_ll)
    println("| Generations = $t")
    println("| Evals UL    = ", nevals_ul)
    println("| Evals LL    = ", nevals_ll)
    @printf("| best F.     = %e\n", best.F)
    @printf("| best f.     = %e\n", best.f)
    
end

function printResults(best::xfg_indiv, P, t, nevals_ul, nevals_ll)
    println("| Generations = $t")
    println("| Evals UL    = ", nevals_ul)
    println("| Evals LL    = ", nevals_ll)
    @printf("| best F.     = %e\n", best.F)
    @printf("| best f.     = %e\n", best.f)
    @printf("| No. vio. G. = %i\n", countViolations(best.G,[]))
    @printf("| No. vio. g. = %i\n", countViolations(best.g,[]))
    
end


function printResults(best::xfgh_indiv, P, t, nevals_ul, nevals_ll)
    println("| Generations = $t")
    println("| Evals UL    = ", nevals_ul)
    println("| Evals LL    = ", nevals_ll)
    @printf("| best f.     = %e\n", best.f)
    @printf("| # vios. UL  = %i\n", countViolations(best.G, best.H))
    @printf("| # vios. LL  = %i\n", countViolations(best.g, best.h))
    
end

function isfeasible(element::xf_indiv)
    return true
end

function isfeasible(element::xfg_indiv)
    return countViolations(element.g, []) == 0
end

function isfeasible(element::xfgh_indiv)
    return countViolations(element.g, element.h) == 0
end