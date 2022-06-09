using TSPLIB
using JuMP, Gurobi#, GLPK

include("Model1.jl")
#include("Model2.jl")
include("LagrangianHeuristic.jl")


function solve_models()
	# Reading the data
	tsp = readTSPLIB(:ulysses16)
	c::Matrix{Float64}= tsp.weights

	# Weights matrix
	m = size(c,1) # m rows
	n = size(c,2) # n column

	V::Vector{Int64} = [i for i=1:n]
	E::Vector{Tuple{Int64, Int64}} = [(i,j) for j=1:n for i=1:m if i<j]

	solverSelected::DataType = Gurobi.Optimizer
	time_limit::Int64 = 600
	model::Model = Model(solverSelected)

	for K::Int64=3:5
		P::Int64 = ceil(n/K)
		println("K = $K,P = $P")

		# Create the model
		#model = model1(model, m, n, K, P, V, E, c)
		#model = model2(model, m, n, K, P, V, E, c)
        # x_best, y_best, UB, LB
		Z_UB = LR_test( model, m, n, K, P, V, E, c)
		#println("LR = $Z_UB")

		# Suppression of the Gurobi verbose mode and use of the callback function
		if solverSelected == Gurobi.Optimizer
			set_optimizer_attribute(model, "LogToConsole", 0)
		end

		set_time_limit_sec(model, time_limit) # Resolution time limit at time_limit seconds

	    # Resolution
	    t = @elapsed optimize!(model)

	    # Display of results
	    status = termination_status(model)

	    if status == MOI.OPTIMAL
	        println("Problem solved to optimality")
	        println("obj = ",objective_value(model)) # Display of the optimum value
	        #println("z = ",value.(model[:z])) # Displaying the values of the descision variables
			#println("y = ",value.(model[:z]))
			println("t = ",t) # CPUt(s) display
	    elseif status == MOI.INFEASIBLE
	        println("Impossible problem")
	    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
	        println("Unbounded problem")
		else
			println(status)
	    end
		empty!(model)
	end
end
solve_models()
#
