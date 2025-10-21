module HierarchicalVariation
    include("parameter_structs.jl")
    include("extinction_probability_multipatch_functions.jl")

    export minimize_extinction_probability, DistributionParams, BaseParams, HierarchicalParams, Run
end