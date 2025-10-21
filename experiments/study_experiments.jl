import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using HierarchicalVariation

function collect_bimodal_runs(seed::Int||Nothing = nothing)
    alpha1s = [0.6, 0.9]
    p1s = 0.5:1/40:1
    base_save_dir = "study/hierarchical/bimodal"
    if seed
        base_save_dir = joinpath(base_save_dir, "seed_$(seed)")
    end

    runs = DistributionParams[]
    
    for alpha1 in alpha1s
        for p1 in p1s
            beta1 = 1 - alpha1
            alpha2 = beta1
            beta2 = alpha1
            save_dir = joinpath(base_save_dir, "alpha1_$(alpha1)__p1_$(p1)")
            push!(runs, DistributionParams(alpha1, beta1, alpha2, beta2, p1, save_dir))
        end
    end
    return runs
end

function collect_central_mode__varying_width_runs(seed::Int||Nothing = nothing)
    
    pairs_ = [(4,2), (7,3), (20,10), (23,7)] 
    p1s = 0.5:1/40:1
    base_save_dir = "study/hierarchical/central_mode__varying_width"
    if seed
        base_save_dir = joinpath(base_save_dir, "seed_$(seed)")
    end

    runs = DistributionParams[]

    for pair in pairs_
        for p1 in p1s
            alpha1 = pair[1]
            beta1 = pair[2]
            
            alpha2 = beta1
            beta2 = alpha1

            save_dir = joinpath(base_save_dir, "alpha1_$(alpha1)__p1_$(p1)")

            push!(runs, DistributionParams(alpha1, beta1, alpha2, beta2, p1, save_dir))
        end
    end
    return runs
end

function collect_extreme_mode__varying_mass_everywhere_else_runs(seed::Int||Nothing = nothing)
    beta1s = [0.05, 0.95]
    p1s = 0.5:1/40:1
    base_save_dir = "study/hierarchical/extreme_mode__varying_mass_everywhere_else"
    if seed
        base_save_dir = joinpath(base_save_dir, "seed_$(seed)")
    end
    
    runs = DistributionParams[seed::Int||Nothing = nothing]
    
    for beta1 in beta1s
        for p1 in p1s
            alpha1 = 1
            alpha2 = beta1
            beta2 = alpha1
            save_dir = joinpath(base_save_dir, "beta1_$(beta1)__p1_$(p1)")
            push!(runs, DistributionParams(alpha1, beta1, alpha2, beta2, p1, save_dir))
        end
    end
    return runs
end

function collect_extreme_mode__varying_slope_runs(seed::Int||Nothing = nothing)
    alpha1_p1_combs = Dict(1.2 => [], 1.8 => 0.7:1/40:1, 2 => 0.7:1/40:1, 3 => 0.6:1/40:1, 4 => 0.5:1/40:1, 5 => 0.5:1/40:1, 7.5 => 0.5:1/40:1, 10 => 0.5:1/40:1)

    base_save_dir = "study/hierarchical/extreme_mode__varying_slope"
    if seed
        base_save_dir = joinpath(base_save_dir, "seed_$(seed)")
    end
    
    runs = DistributionParams[]

    for (alpha1, p1s) in pairs(alpha1_p1_combs)
        for p1 in p1s
            beta1 = 1
            
            alpha2 = beta1
            beta2 = alpha1

            save_dir = joinpath(base_save_dir, "alpha1_$(alpha1)__p1_$(p1)")

            push!(runs, DistributionParams(alpha1, beta1, alpha2, beta2, p1, save_dir))
        end
    end
    return runs
end

function collect_all_runs_for_seed(seed::Int||Nothing = nothing)
    dist_params_runs = vcat(
        collect_bimodal_runs(),
        collect_central_mode__varying_width_runs(),
        collect_extreme_mode__varying_mass_everywhere_else_runs(),
        collect_extreme_mode__varying_slope_runs()
    )
    base_params = BaseParams(fecundity=4, delta=0.1, population_size=100, use_educated_guess=true, analyze_centers=true, tol=1e-8, seed=seed, max_runs=1000, p_cutoff=0.05)
    hierarchical_params = HierarchicalParams(q0=0.5, num_generations=75, num_runs=100)
    runs = [Run(base=base_params, dist=el, hier=hierarchical_params) for el in dist_params_runs]
    return runs
end

function collect_all_runs()
    runs = vcat([collect_all_runs_for_seed(seed) for seed in 0:6]...)
    return runs
end



