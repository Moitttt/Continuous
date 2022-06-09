function t()
    c = Matrix{Float64}(undef, 2, 2)
    c[1,1]=23
    c[1,2]=53
    c[2,1]=98
    c[2,2]=1000
    print(c[:,1])
    maximum(c)
end
