function is_better_mass(S1, S2, searchType=:minimize)
    # S1 is better than S2
    α = 10.0

    m1 = S1.F + α*S1.f
    m2 = S2.F + α*S2.f

    if searchType == :minimize
        return  m1 < m2
    end

    return m1 > m2
end

function Selection(Old::xf_indiv, New::xf_indiv, searchType::Symbol=:minimize; leq::Bool=false)
    return is_better_mass(New, Old, searchType)
end

# Deb rules (selection)
function Selection(Old::xfgh_indiv, New::xfgh_indiv, searchType::Symbol=:minimize; leq::Bool=false)

    old_vio = violationsSum(Old.g, Old.h) + violationsSum(Old.G, Old.H)
    new_vio = violationsSum(New.g, New.h) + violationsSum(New.G, New.H)

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    return is_better_mass(New, Old, searchType)
end

# Deb rules (selection)
function Selection(Old::xfg_indiv, New::xfg_indiv, searchType::Symbol=:minimize; leq::Bool=false)
    old_vio = violationsSum(Old.g, []) + violationsSum(Old.g, [])
    new_vio = violationsSum(New.g, []) + violationsSum(Old.G, [])

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    return is_better_mass(New, Old, searchType)
end

function getBest(Population, searchType::Symbol = :minimize)
    best = Population[1]

    for i = 2:length(Population)
        if Selection(best, Population[i])
            best = Population[i]
        end
    end

    return best
end

function getWorst(Population, searchType::Symbol = :minimize)
    worst = Population[1]

    for i = 2:length(Population)
        if Selection(Population[i], worst)
            worst = Population[i]
        end
    end

    return worst
end

function getWorstInd(Population, searchType::Symbol = :minimize)
    worst = 1

    for i = 2:length(Population)
        if Selection(Population[i], Population[worst], searchType)
            worst = i
        end
    end

    return worst
end

function is_better(x, y)
    # x better than y
    return Selection(y, x)
end

function getBestInd(Population, searchType::Symbol = :minimize)
    j = 1

    for i = 2:length(Population)
        if Selection(Population[j], Population[i])
            j = i
        end
    end

    return j
end

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Float64, fResult::Float64)
    return xf_indiv(x, y, FResult, fResult)
end

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Tuple{Float64,Array{Float64,1}}, fResult::Tuple{Float64,Array{Float64,1}})
    F, G = FResult
    f, g = fResult
    return xfg_indiv(x, y, F, f, G, g)
end

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}}, fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    F, G, H = FResult
    f, g, h = fResult
    return xfgh_indiv(x, y, F, f, G, g, H, h)
end

function inferType(fVal::Tuple{Float64})
    return xf_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1}})
    return xfg_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    return xfgh_indiv
end