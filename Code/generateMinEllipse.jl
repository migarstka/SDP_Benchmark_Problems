
using JuMP
using Mosek
using PyPlot

# We will use three ellipses
m = 3
# Two "simple" ones
As = [[2.0  0.0;0.0  1.0],[1.0  0.0;0.0  3.0]]
# and a random one
randA = rand(2,2)
push!(As, (randA' * randA) * (rand()*2+1))

# We change the weights to see different solutions, if they exist
W = [1.0 0.0;
     0.0 1.0]

model = Model(solver=MosekSolver())
@variable(model, X[1:2,1:2], SDP)
@objective(model, Min, trace(W*X))
for i = 1:m
    @SDconstraint(model, X >= As[i])
end
JuMP.solve(model)

X_val = getValue(X)
println(X_val)

# Setup the figure
fig = PyPlot.figure(1,facecolor="white",figsize=(12,5))


# Draw provided ellipses
for i = 1:m
    xs = Float64[]
    ys = Float64[]
    for angle in linspace(0, 2*pi, 100)
        u = [cos(angle),sin(angle)]
        x = As[i] * u
        push!(xs, x[1])
        push!(ys, x[2])
    end
    PyPlot.plot(xs, ys, "b", linewidth=2.0)
end

# Draw bounding ellipse
xs = Float64[]
ys = Float64[]
for angle in linspace(0, 2*pi, 100)
    u = [cos(angle),sin(angle)]
    x = X_val * u
    push!(xs, x[1])
    push!(ys, x[2])
end
PyPlot.plot(xs, ys, "r", linewidth=2.0)
PyPlot.axis("equal")
