# SDP Benchmark Problems
In the following a set of benchmark problems is presented. These convex optimization problems all involve positive semidefinite constraints and occur in various fields, 
such as portfolio optimization, graph theory, robust control, or polynomial optimization. These problems were collected to benchmark conic solvers 
against each other and to study the performance gains of various extensions to our solver. The problems are in the following format:
```
min 1/2x'Px+q'x+r
s.t. Ax+s=b,
     s in K
```
with decision variables x in Rn, s in Rm and data matrices P=P'>=0, q in Rn, A in Rm×n, and b in Rm. The convex cone K is a composition of the zero cone, 
the non-negative orthant, a set of second order cones, and a set of positive semidefinite cones. The following list gives an overview over the considered problems:
1. Random SDP with Quadratic Objective
2. Nearest Correlation Matrix
3. Smallest Sphere around multiple Ellipsoids
4. LMI-based Robust Control Problems
5. Semidefinite Relaxation of MIQO Problems
6. Lovasz Theta Function in Graph Theory

The problems are available as Julia JLD files and Matlab MAT files. Each file contains at least the following variables (some contain extra problem information):

Variable | Description |
--- | --- |
m,n | problem dimension | 
A,b | constraint data | 
P,q,r | objective function data | 
objTrue | solution to problem from MOSEK with standard tolerance |
problemType | descriptive name of the problem | 
problemName | short problem tag with number |
Kf, Kl, Kq, Ks | dimensions of the cone K |

## Installation / Usage

## Tasks / Future Work

## Licence
This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact
Send an email :email: to [Michael Garstka](mailto:michael.garstka@eng.ox.ac.uk) :rocket:!	
