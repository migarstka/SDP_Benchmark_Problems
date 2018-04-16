# SDP Benchmark Problems
In the following a set of benchmark problems is presented. These convex optimization problems all involve positive semidefinite constraints and occur in various fields, 
such as portfolio optimization, graph theory, or robust control. These problems were collected to benchmark conic solvers 
against each other and to study the performance gains of various extensions to conic solvers. The problems are given in the following format:
```
min 1/2x'Px+q'x+r
s.t. Ax+s=b,
     s in K
```
with decision variables `x ϵ R^n`, `s ϵ R^m` and data matrices `P=P'>=0`, `q ϵ R^n`, `A ϵ R^(m×n)`, and `b ϵ R^m`. The convex cone K is a composition of the zero cone, 
the non-negative orthant, a set of second order cones, and a set of positive semidefinite cones. The following list gives an overview over the considered problems (problem types with nonzero `P` in italic):
1. _Random SDP with Quadratic Objective_
2. _Nearest Correlation Matrix_
3. Smallest Sphere around multiple Ellipsoids
4. LMI-based Robust Control Problems
5. Semidefinite Relaxation of MIQO Problems
6. Lovasz Theta Function in Graph Theory

The problems are available as Julia JLD files and Matlab MAT files (see respective directories). Each file describes one problem and contains at least the following variables (some contain extra problem information):

Variable | Description |Type
--- | --- | --- |
m,n | problem dimension | Int |
A,b | constraint data |  SparseMatrixCSC, Vector |
P,q,r | objective function data |  SparseMatrixCSC, Vector, Float |
objTrue | solution to problem from MOSEK with standard tolerance | Float |
problemType | descriptive name of the problem | String |
problemName | short problem tag with number | String |
Kf | Number of variables in zero-cone | Int |
Kl | Number of variables in nonnegative orthant | Int |
Kq | Dimensions of second-order cone variables | Array{Int} |
Ks | Dimensions of semidefinite cone variables | Array{Int} |

The `Code`folder contains the Julia scripts used to generate the `.jld`-files. Furthermore, the script `convertToMAT.jl`can be used to convert `.jld`-files into Matlab `.mat`-files.

## Installation / Usage
- In Julia the .jld files can be loaded by using the `JLD` package. An example for a random SDP with quadratic objective is given here:
```julia
using JLD

filePath = #insert filepath here, e.g. "/Users/User1/DataFiles/Julia/SDPQuad/SDPQuad01.jld"
data = JLD.load(filePath)
P = data["P"]
q = data["q"]
r = data["r"]

A = data["A"]
b = data["b"]

m = data["m"]
n = data["n"]

Kf = data["Kf"]
Kl = data["Kl"]
Kq = data["Kq"]
Ks = data["Ks"]

objTrue = data["objTrue"]
```

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
For questions or suggestions for other psd-problems (especially with quadratic cost), send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!	
