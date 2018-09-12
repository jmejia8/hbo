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

function feasible_map(fnum)
    p = q = 3
    r = 2
    s = 0
    n = 100

    f, F, lower_D, upper_D, lower_bounds, upper_bounds = getBilevel(fnum)

    Nsamples = 1000
    a = upper_bounds[1,:]
    b = upper_bounds[2,:]
    
    x = a' .+ (b -a)' .* rand(Nsamples, upper_D)


    a = lower_bounds[1,:]; b = lower_bounds[2,:]
    y = zeros(Nsamples, lower_D);

    if fnum == 1
        y[:, q+1:end] = atan.(x[:,p+1:end])
    elseif fnum == 2
        y[:, q+1:end] = exp.(x[:,p+1:end])
    elseif fnum == 3
        y[:, q+1:end] = atan.(x[:,p+1:end].^2)
    elseif fnum == 4
        y[:, q+1:end] = exp.(abs.(x[:,p+1:end])) - 1.0
    elseif fnum == 5
        y[:,1:q] = 1.0
        y[:, q+1:end] = sqrt.(abs.(x[:,p+1:end]))
    elseif fnum == 6
        q = 1; s = 2;
        y[:, q+s+1:end] = x[:,p+1:end]
    elseif fnum == 7
        y[:, q+1:end] = exp.(x[:,p+1:end])
    elseif fnum == 8
        y[:, q+1:end] = (x[:,p+1:end]).^(1/3)
    end

    return x, y
end

function test()

    for fnum = 1:8
      f, F, D_ll, D_ul, bounds_ll, bounds_ul = getBilevel(fnum)
      x, y, best, nevals, _,_ = hbo(F, f, D_ul, D_ll, bounds_ul, bounds_ll; showResults = false)
      @printf("SMD%d \t F = %e \t f = %e  \t nfes = %d \n", fnum, best.F, best.f, nevals)
      # println("x: ", x)
      # println("y: ", y)

      # if fnum == 7
      #     println(">>>>> F = ", smd7(x, y))
      #     println(">>>>> f = ", smd7_ll(x, y))
      # end
      # break
    end

    return 
end
