using JuMP, Gurobi #, GLPK

# A natural formulation of the problem of partitioning the nodes of a graph
function model1(model::Model, m::Int64, n::Int64, K::Int64, P::Int64, V::Vector{Int64}, E::Vector{Tuple{Int64, Int64}}, c::Matrix{Float64})
	# Descision variables
	@variable(model,y[1:m,1:n,1:K], Bin)# Equal to one if and only if edge e ∈ E joins two nodes in subset k
    @variable(model,z[i in V,1:K], Bin)	# Equal to one if and only if i ∈ V belongs to subset k

	# Objective function
    @objective(model, Max, sum(sum(c[e...]*y[e...,k] for e in E) for k=1:K))# Maximize the sum of the costs of the edges inside the subsets

	# Constraints
	@constraint(model, c1[i in V], sum(z[i,k] for k=1:K)==1)		# Ensure each node is in exactly one subset
	@constraint(model, c2[(i,j) in E,k=1:K], y[i,j,k] <= z[i,k])	# link y and z variables
	@constraint(model, c3[(i,j) in E,k=1:K], y[i,j,k] <= z[j,k])	# link y and z variables
	@constraint(model, c4[k=1:K], 1 <= sum(z[i,k] for i in V) <= P)	# Make sure subsets are not empty and contain at most P nodes

	return model
end
