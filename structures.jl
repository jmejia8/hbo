mutable struct xf_indiv # Single Objective
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
end

mutable struct xfg_indiv # Single Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
end

mutable struct xfgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end

mutable struct xFgh_indiv # multi Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Vector{Float64}
    f::Vector{Float64}
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end