abstract type Coalescent end

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

mutable struct KingmanCoalescent <: Coalescent
	n::Int
	N::Float64
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
mutable struct SEFCoalescent <: Coalescent
	n::Int
	β::Distribution
	ρ::Float64
end
