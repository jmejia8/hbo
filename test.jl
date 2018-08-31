using BilevelBenchmark

include("hbo.jl")

# configures problem
function getBilevel(fnum::Int)

    if fnum < 0
       F_(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
       f_(x, y) = sum((x.^2 + y.^2 - 1.0).^2)
       D_ul = 1
       D_ll = 1

       bounds_ul = Array{Float64, 2}([-1  1])'
       
       bounds_ll = bounds_ul

       return f_, F_, D_ll, D_ul, bounds_ll, bounds_ul
    end
    
    D_ul = D_ll = 5

    bounds_ul, bounds_ll = bilevel_ranges(D_ul, D_ll, fnum)

    # leader
    F(x::Array{Float64}, y::Array{Float64}) =  bilevel_leader(x, y, fnum)

    # follower
    f(x::Array{Float64}, y::Array{Float64}) =bilevel_follower(x, y, fnum)
    
    return f, F, D_ll, D_ul, bounds_ll, bounds_ul
end

function test()

    for fnum = 1:8
      f, F, D_ll, D_ul, bounds_ll, bounds_ul = getBilevel(fnum)
      x, y, best, nevals, _ = hbo(F, f, D_ul, D_ll, bounds_ul, bounds_ll; showResults = false)
      @printf("SMD%d \t F = %e \t f = %e  \t nfes = %d \n", fnum, best.F, best.f, nevals)
    end

    return 
end
