using BackwardCoalescent
using TreeTools

n = 50

# For the EF coalescent
β = 0.2
ρ = 0.025
EFC = EFCoalescent(n, β, ρ)
@info "EF coalescent: effective population size Ne=$(BackwardCoalescent.Ne(EFC))"

t_ef = BackwardCoalescent.genealogy(EFC)
write_newick("example_trees/tree_ef.nwk", t_ef)

# Kingman coalescent
N = 1_000
KC = KingmanCoalescent(n, N)

t_k = BackwardCoalescent.genealogy(KC)
write_newick("example_trees/tree_kingman.newk", t_k)
