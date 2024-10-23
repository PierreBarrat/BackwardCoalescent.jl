abstract type Coalescent end


@kwdef mutable struct KingmanCoalescent <: Coalescent
    n::Int = 2
    N::Float64 = 1
end
function choose_event(C::KingmanCoalescent)
    @assert C.n > 1 "Cannot choose coalescence event for one lineage." C
    τ = rand(Exponential(2*C.N/C.n/(C.n-1)))
    return 2, τ
end

"""
    mutable struct YuleCoalescent

```julia
n::Int = 2
b::Float64 = 1 # birth rate
```
Rate of coalescence `(n-1)b`.
The expected heigh of the tree is (I think) `~log(n)/b`
"""
@kwdef mutable struct YuleCoalescent <: Coalescent
    n::Int = 2
    b::Float64 = 1 # birth rate
end
function choose_event(Y::YuleCoalescent)
    @assert Y.n > 1 "Cannot choose coalescence event for one lineage." Y
    τ = rand(Exponential(1/(Y.n-1)/Y.b))
    return 2, τ
end

"""
	EFCoalescent

```
mutable struct EFCoalescent
	n::Int
	β::Float64
	ρ::Float64
end
```
"""
mutable struct EFCoalescent <: Coalescent
	n::Int
	β::Float64
	ρ::Float64
end
Ne(C::EFCoalescent) = 1/C.ρ/C.β/C.β
function choose_event(C::EFCoalescent)
    @assert C.n > 1 "Cannot choose coalescence event for one lineage." C
    τ = 0.
    merger = 0
    coin = Binomial(C.n, C.β)
    while merger < 2
        τ += rand(Exponential(1/C.ρ))
        merger = rand(coin)
    end
    return merger, τ
end

"""
	mutable struct SEFCoalescent

Stochastic EF coalescent. At each step, β is chosen from a distribution.

## Fields
```
	n::Int
	β::Distribution
	ρ::Float64
```
"""
mutable struct SEFCoalescent{T<:UnivariateDistribution} <: Coalescent
	n::Int
	β::T
	ρ::Float64
end
function choose_event(C::SEFCoalescent)
    @assert C.n > 1 "Cannot choose coalescence event for one lineage." C
    τ = 0.
    merger = 0
    while merger < 2
        τ += rand(Exponential(1/C.ρ))
        βval = rand(C.β)
        coin = Binomial(C.n, βval)
        merger = rand(coin)
        @debug "`SEFCoalescent`: using β=$(βval)"
    end
    return merger, τ
end
