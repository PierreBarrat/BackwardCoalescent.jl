"""
	genealogy(C::Coalescent; coalescence_times = false)

Return a tree sampled from the coalescent `C`. If `coalescence_times`, also return the
times `Tn` during which `n` lineages are present, in the form of a dictionary `n => Tn`.

## Examples
For the EF coalescent with parameters `β` and `ρ` and `n=100` lineages:
```
tree = genealogy(EFCoalescent(100, β, ρ))
```
"""
function genealogy(C::Coalescent; coalescence_times = false)
	# Leaves
	nodes = [TreeNode(label="$i") for i in 1:C.n]
	reset_id()
	r, T = genealogy!(nodes, deepcopy(C))
	if coalescence_times
		return node2tree(r), T
	else
		return node2tree(r)
	end
end

function genealogy!(nodes::Vector{<:TreeNode}, C::Coalescent, T = Dict())
	ids, τ = choose_event(C) # nodes to be coalesced and time to event
	@debug "Coalescence of nodes $([n.label for n in nodes[ids]]) with time $τ"
	T[C.n] = τ
	new_nodes = merge!(nodes, ids, τ)
	k = length(ids) # size of the coalescence
	C.n += -k + 1 # k lineages merge into one
	@assert C.n == length(nodes) "Inconsistent number of remaining lineages: \
		$(C.n) vs. $(length(nodes))"

	if C.n == 1
		return nodes[1], T
	else
		return genealogy!(nodes, C, T)
	end
end

"""
	merge!(nodes, ids, τ)

Merge `nodes[ids]` into one common ancestor, with branch length `τ`, and remove them from
  `nodes`. The new ancestor is pushed at the end of `nodes`.
"""
function merge!(nodes, ids, τ)
	# Updating branch length
	for n in nodes
		n.tau = ismissing(n.tau) ? τ : n.tau + τ
	end
	# Merging nodes
	r = TreeNode(label="internal_$(get_id())")
	for n in nodes[ids]
		graftnode!(r, n)
	end
	deleteat!(nodes, ids)
	push!(nodes, r)
	return nothing
end

function choose_event(C::EFCoalescent)
	@assert C.n > 1 "Cannot choose coalescence event for one lineage." C
	τ = 0.
	ids = Int[]
	k = 0
	while k < 2
		ids = Int[]
		τ += rand(Exponential(1/C.ρ))
		for i in 1:C.n
			if rand() < C.β
				push!(ids, i)
			end
		end
		k = length(ids)
		ids
	end
	return ids, τ
end

function choose_event(C::SEFCoalescent)
	@assert C.n > 1 "Cannot choose coalescence event for one lineage." C
	τ = 0.
	ids = Int[]
	k = 0
	while k < 2
		ids = Int[]
		τ += rand(Exponential(1/C.ρ))
		βval = rand(C.β)
		@debug "`SEFCoalescent`: using β=$(βval)"
		for i in 1:C.n
			if rand() < βval
				push!(ids, i)
			end
		end
		k = length(ids)
		ids
	end
	return ids, τ
end

function choose_event(C::KingmanCoalescent)
	@assert C.n > 1 "Cannot choose coalescence event for one lineage." C
	τ = rand(Exponential(2*C.N/C.n/(C.n-1)))
	ids = rand(1:C.n, 2)
	while ids[2] == ids[1]
		ids[2] = rand(1:C.n)
	end
	return sort(ids), τ
end

let i = 0
	global get_id() = (i+=1)
	global reset_id() = (i=0)
end
