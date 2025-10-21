import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using HierarchicalVariation
using Base.Threads

include(joinpath(@__DIR__, "..", "experiments", "study_experiments.jl"))

function has_been_previously_run(run::Run)
    return isdir(run.dist.save_dir) && length(readdir(run.dist.save_dir)) == 2
end

function run_experiments(runs::Vector{Run})
    @threads for i in eachindex(runs)
        if !has_been_previously_run(runs[i])
            minimize_extinction_probability(runs[i])
        end
    end
end

study_experiments = collect_all_runs()
run_experiments(study_experiments)