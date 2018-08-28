using BilevelBenchmark

include("hbo.jl")

function genBounds(uBounds, lBounds, p, q, r, s)

    return [repmat(uBounds[:,1], 1, p) repmat(uBounds[:,2], 1, r)],
           [repmat(lBounds[:,1], 1, q) repmat(lBounds[:,2], 1, r + s)]
end

# configures problem
function getBilevel(fnum::Int)

    if fnum == 1 || fnum == 3
        ub = [-5 10; -5 10.0]
        lb = [-5 10; -π/2 π/2]
    elseif fnum == 2 || fnum == 7
        ub = [-5 10; -5 1.0]
        lb = [ -5.0  10; 0.0  e]
    elseif fnum == 4
        ub = [-5.0 10; -1  1]
        lb = [-5.0 10;  0  e]
    elseif fnum == 5  || fnum == 6 || fnum == 8
        ub = [-5.0 10; -5.0  10.0]
        lb = [-5.0 10; -5.0  10.0]
    end

    if fnum == 6
        p, q, r, s = 3, 1, 2, 2
    else
        p, q, r, s = 3, 3, 2, 0
    end
    
    upper_D, lower_D = p + r, q + r + s

    upper_bounds, lower_bounds = genBounds( ub', lb', p, q, r, s ) 

    # leader
    F(x::Array{Float64}, y::Array{Float64}) =  SMD_leader(x, y, fnum, p, q, r, s)

    # follower
    f(x::Array{Float64}, y::Array{Float64}) =SMD_follower(x, y, fnum, p, q, r, s)
    
    return f, F, lower_D, upper_D, lower_bounds, upper_bounds
end

function test()
    fnum = 1

    f, F, lower_D, upper_D, lower_bounds, upper_bounds = getBilevel(fnum)
    hbo(F, f, upper_D, lower_D, upper_bounds, lower_bounds)
end

test()