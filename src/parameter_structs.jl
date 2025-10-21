struct DistributionParams
    alpha1::Number
    beta1::Number
    alpha2::Number
    beta2::Number
    p1::Number
    save_dir::String
end

Base.@kwdef struct BaseParams
    fecundity::Int64 = 4
    delta::Float64 = 0.1
    population_size::Int64 = 100
    use_educated_guess::Bool = true
    analyze_centers::Bool = true
    tol::Float64 = 1e-8
    seed::Number = 0
    max_runs::Int = 1000
    p_cutoff::Float64 = 0.05
end

Base.@kwdef struct HierarchicalParams
    q0::Float64 = 0.5
    num_generations::Int64 = 75
    num_runs::Int64 = 100
end

Base.@kwdef struct Run
    dist::DistributionParams
    base::BaseParams = BaseParams()
    hier::Union{Nothing,HierarchicalParams} = nothing  # nothing => non-hierarchical
end