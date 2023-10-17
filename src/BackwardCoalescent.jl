module BackwardCoalescent

using Distributions
using StatsBase
using TreeTools

include("objects.jl")
export Coalescent
export EFCoalescent, KingmanCoalescent, SEFCoalescent, YuleCoalescent

include("evolve.jl")
export coalescence_times, genealogy

end # module
