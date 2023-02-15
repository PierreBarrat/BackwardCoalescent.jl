module BackwardCoalescent

using Distributions
using TreeTools

include("objects.jl")
export Coalescent
export EFCoalescent, KingmanCoalescent, SEFCoalescent

include("evolve.jl")
export coalescence_times, genealogy

end # module
