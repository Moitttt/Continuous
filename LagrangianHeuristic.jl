using TSPLIB
using JuMP, Gurobi#, GLPK

# A Lagrangian Relaxation of model 2
function LR_model2(model::Model, m::Int64, n::Int64, K::Int64, P::Int64, V::Vector{Int64}, E::Vector{Tuple{Int64, Int64}}, c::Matrix{Float64}, λ::Vector{Float64})
	# Descision variables
	@variable(model,y[i in 1:m,j in 1:n,1:min(i,j)], Bin)# Equal to one if and only if edge e ∈ E joins two nodes in subset k
    @variable(model,z[i in V,1:i], Bin)	# Equal to one if and only if i ∈ V belongs to subset k

	# Objective function
    @objective(model, Max, sum(sum(c[i,j]*y[i,j,k] for k=1:min(i,j)) for (i,j) in E)  -  sum( λ[i] * (1-sum(z[i,k] for k=1:i)) for i in V) - sum(λ)*(K-sum(z[k,k] for k in V)))

	# Constraints
	@constraint(model, c2[(i,j) in E,k=1:min(i,j)], y[i,j,k] <= z[i,k])	# link y and z variables
	@constraint(model, c3[(i,j) in E,k=1:min(i,j)], y[i,j,k] <= z[j,k])	# link y and z variables
	@constraint(model, c4[k in V], sum(z[i,k] for i in V if i>k) <= (P-1)*z[k,k])	# Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets with at most P nodes.

	return model
end

function LR_test(model::Model, m::Int64, n::Int64, K::Int64, P::Int64, V::Vector{Int64}, E::Vector{Tuple{Int64, Int64}}, c::Matrix{Float64})
    λ = zeros(size(V))
    #λ₂ = zeros(1)
    θ = 1.0
    Zlb = Inf
    Zub = -Inf

    for g=1:100
        #println(g)
        println("λ = $λ")
        #println("Zlb = $Zlb")

        model = LR_model2(model, m, n, K, P, V, E, c, λ)
        set_optimizer_attribute(model, "LogToConsole", 0)
        optimize!(model)

        Zlbtemp = objective_value(model)
        ylb = value.(model[:y])
		z = value.(model[:z])

		"""
		for i=1:m
            for j=1:n
                for k=1:min(i,j)
                    if ylb[i, j, k] == 1
                        println(ylb[i, j, k])
                    end
                end
            end
        end
		"""

        Zubtemp, yub, z = upper_bound(z, ylb, m, n, K, P, V, E, c)

        if Zub < Zubtemp
            Zub = Zubtemp
            besty = yub
            bestz = z
        end

        if Zlb > Zlbtemp
            Zlb = Zlbtemp
        end

        sommez = zeros(m)
        for k=1:m
            for i=k:n
                sommez[k]+=z[i,k]
            end
        end
		residual = 1 .- sommez
		t = θ * (Zub - Zlb) / sum(residual.^2)
		λ = λ + vec(t * residual)

		empty!(model)
    end
    println("res: $Zub")
    return Zub
end

function is_feasible(z, y, m, n, K, P, V, E, c)

    for i in V
        if !(sum(z[i,k] for k=1:i)==1)
            #println("hello1")
            return false
        end
    end


	for (i,j) in E
        for k=1:min(i,j)
            if !(y[i,j,k] <= z[i,k] &&  y[i,j,k] <= z[j,k])
                #println("i , j : $i $j $k")
				#println("yyyyyyyyyy")
                return false
            end
        end
    end
    for k in V
        if !(sum(z[i,k] for i in V if i>k) <= (P-1)*z[k,k])
            """println("k: ",k)
            println("sum =", (sum(z[i,k] for i in V if i>k)))
            println("coooo = ",(P-1)*z[k,k])
			println("zzzzzzzzzzzzzzz")"""
            return false
        end
    end

    if !(sum(z[k,k] for k in V) == K)
        #println("hello2")
        return false
    end
    return true
end

function is_feasible2(z, y, m, n, K, P, V, E, c)
    notfeasible = []
    notfeasible2=[]
    for i in V
        if !(sum(z[i,k] for k=1:i)==1)
            temp = []
            for k=1:i
                if z[i,k]==1
                    push!(temp,(i,k))
                end
            end
            push!(notfeasible, temp)
        end
    end

	if !(sum(z[k,k] for k in V) == K)

		temp1 = []
		temp0 = []
		for k in V
			if z[k,k] == 1
				push!(temp1, k)
			else
				push!(temp0, k)
			end
		end
		notfeasible2=[temp1,temp0]
	end

    return notfeasible,notfeasible2
end


function upper_bound(z, y, m, n, K, P, V, E, c)
  # Computing y, given z
    notfeasible,notfeasible2 = is_feasible2(z, y, m, n, K, P, V, E, c)
    """println("notfeasible : $notfeasible")
    println("notfeasible2 : $notfeasible2")"""

	# Constraints (1)
	len1=length(notfeasible)
	for it=1:len1
		len2=length(notfeasible[it])-1
        for it2=1:len2
			rnd = rand(notfeasible[it])
			z[rnd]=0
			deleteat!(notfeasible[it],findall(x->x==rnd,notfeasible[it]))
		end
	end

    notfeasible,notfeasible2 = is_feasible2(z, y, m, n, K, P, V, E, c)

	#println(notfeasible2)
    #println(length(notfeasible2))
	#println(K)

	# Constraints (5)
	nbK= length(notfeasible2[1])
	if nbK < K
        len1 = K - nbK
        for it=1:len1
        	i = rand(notfeasible2[2])
			#println("i = $i")
        	temp = []
            println(z)
            for k=1:i
                if z[i,k]==1
                    push!(temp,(i,k))
                end
            end
        	rnd = rand(temp)
			#println("rnd = $rnd, temp = $temp")

			#println(z)
        	z[rnd]=0
        	z[i,i]=1
			#println(z)

        	deleteat!(notfeasible2[2],findall(x->x==i,notfeasible2[2]))
	        push!(notfeasible2[1],i)
        end

	elseif nbK > K
		len2 = nbK - K
		for it=1:len2
            #println(notfeasible2[1])

			i = rand(notfeasible2[1])
            while i ==1
                i = rand(notfeasible2[1])
            end
            #println(i)
			z[i,i]=0
			rnd=rand(1:i)
			while rnd == i
				rnd=rand(1:i)
			end
			z[i,rnd]=1

			deleteat!(notfeasible2[1],findall(x->x==i,notfeasible2[1]))
			push!(notfeasible2[2],i)
        end
	end

	# Constraints (2) and (3)
	for i=1:m
		for j=1:n
			for k=1:min(i,j)
				if !(y[i,j,k] <= z[i,k] && y[i,j,k] <= z[j,k])
					y[i,j,k] =0
				end
			end
		end
	end

    """if is_feasible(z, y, m, n, K, P, V, E, c)
        println("VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV")
    end"""

  # Computing Z
    Z = 0.0
    for i=1:m
        for j=1:n
            for k=1:min(i,j)
                Z = Z + c[i,j]*y[i,j,k]
            end
        end
    end
    #println("Z : $Z")

  return Z, y, z
end
