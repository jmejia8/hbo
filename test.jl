using BilevelBenchmark
import Printf.@printf
import Random.randperm
using CSVanalyzer
import DelimitedFiles.writedlm

include("hbo.jl")



# configures problem
function getBilevel(fnum::Int)

    if fnum == -1
       F_(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
       f_(x, y) = sum((x.^2 + y.^2 - 1.0).^2)
       D_ul = 1
       D_ll = 1

       bounds_ul = Array{Float64, 2}([-1.0  1.0])'
       
       bounds_ll = bounds_ul

       return f_, F_, D_ll, D_ul, bounds_ll, bounds_ul
    elseif fnum == -2
        g(fxy) = cos(fxy/3) + sin(fxy)

        f1_(x, y) = sum(y.^2)/2 + dot(y, x)
        F1_(x, y) = g(f1_(x, y))

        D_ul = 1
        D_ll = 1

        bounds_ul = Array{Float64, 2}([-5  5])'
       
        bounds_ll = bounds_ul

        return f1_, F1_, D_ll, D_ul, bounds_ll, bounds_ul
    end
    
    D_ul = D_ll = 5

    bounds_ul, bounds_ll = bilevel_ranges(D_ul, D_ll, fnum)

    # leader
    F(x::Array{Float64}, y::Array{Float64}) =  bilevel_leader(x, y, fnum)

    # follower
    f(x::Array{Float64}, y::Array{Float64}) =bilevel_follower(x, y, fnum)
    
    return f, F, D_ll, D_ul, bounds_ll, bounds_ul
end

function feasible_map(fnum;nsamples=500)
    p = q = 3
    r = 2
    s = 0
    n = 100

    f, F, lower_D, upper_D, lower_bounds, upper_bounds = getBilevel(fnum)

    Nsamples = nsamples
    a = upper_bounds[1,:]
    b = upper_bounds[2,:]
   
    # k = 5
    # g = linspace(0, 1, k)
    # R = zeros(upper_D^k, upper_D)
    
    # l = 1
    # for i = 1:upper_D^k
    #     for j = 1:upper_D
    #         R[i, j] = g[l]
    #     end
    #     l = 1 + (l % k)

    # end


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
    elseif fnum == -2
        y = -x
    elseif fnum == -1
        θ = linspace(0,2π, Nsamples)
        x[:,1] = cos.(θ)
        y[:,1] = sin.(θ)
    end

    return x, y
end

function best_known_values(fnum)
    if fnum==-2
        return -8.0, -1.87868481483642
    end

    return 1.9209739856726114e-6,  0.8028017347965496
end

function all_map(fnum)
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
    y = a' .+ (b -a)' .* rand(Nsamples, lower_D)

    return x, y
end

function test()

    for fnum = 1:8
        f, F, D_ll, D_ul, bounds_ll, bounds_ul = getBilevel(fnum)
        x, y, best, nevals, _,_ = hbo(F, f, bounds_ul, bounds_ll; showResults = false)
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


function makeStats(nruns = 31)
    configure()
    nfuns = 8

    errors_UL = zeros(nfuns, nruns)
    errors_LL = zeros(nfuns, nruns)

    evals_UL = zeros(nfuns, nruns)
    evals_LL = zeros(nfuns, nruns)

    for fnum = 1:nfuns

        for i = 1:nruns
            f, F, D_ll, D_ul, bounds_ll, bounds_ul = getBilevel(fnum)
            
            x, y, best, nevals, _,_ = hbo(F, f, bounds_ul, bounds_ll; showResults = false)
            
            @printf("SMD%d \t F = %e \t f = %e  \t nfes = %d \n", fnum, best.F, best.f, nevals)

            errors_UL[fnum, i] = best.F; errors_LL[fnum, i] = best.f
            evals_UL[fnum, i] = evals_LL[fnum, i] = nevals

        end

        println("-----------------------------------------------")

    end

    writedlm("output/accuracy_UL.csv", errors_UL, ',')
    writedlm("output/accuracy_LL.csv", errors_LL, ',')

    writedlm("output/evals_UL.csv", evals_UL, ',')
    writedlm("output/evals_LL.csv", evals_LL, ',')

    ##
    println("Upper level")
    statsToLatex("output/accuracy_UL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")

    println("Lower level")
    statsToLatex("output/accuracy_LL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")


    println("Evaluations")
    statsToLatex("output/evals_UL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")

end

# makeStats()
test()
