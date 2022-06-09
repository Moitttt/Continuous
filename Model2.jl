using TSPLIB
using JuMP, Gurobi, GLPK

# A new formulation based on the notion of representatives of the problem of partitioning the nodes of a graph
function model2(model::Model,solverSelected::DataType, m::Int64, n::Int64, K::Int64, P::Int64, V::Vector{Int64}, E::Vector{Tuple{Int64, Int64}}, c::Matrix{Float64})
	# Descision variables
	@variable(model,y[i in 1:m,j in 1:n,1:min(i,j)], Bin)# Equal to one if and only if edge e ∈ E joins two nodes in subset k
    @variable(model,z[i in V,1:i], Bin)	# Equal to one if and only if i ∈ V belongs to subset k

	# Objective function
    @objective(model, Max, sum(sum(c[i,j]*y[i,j,k] for k=1:min(i,j)) for (i,j) in E))# Maximize the sum of the costs of the edges inside the subsets

	# Constraints
	@constraint(model, c1[i in V], sum(z[i,k] for k=1:i)==1)			# Ensure each node is in exactly one subset
	@constraint(model, c2[(i,j) in E,k=1:min(i,j)], y[i,j,k] <= z[i,k])	# link y and z variables
	@constraint(model, c3[(i,j) in E,k=1:min(i,j)], y[i,j,k] <= z[j,k])	# link y and z variables
	@constraint(model, c4[k in V], sum(z[i,k] for i in V if i>k) <= (P-1)*z[k,k])	# Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets with at most P nodes.
	@constraint(model, c5, sum(z[k,k] for k in V) == K)					# Make sure subsets are not empty and contain at most P nodesmake sure we have exactly K subsets in the partition.

	return model
end
