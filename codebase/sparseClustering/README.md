Sparse convex optimization for Clustering
=======================

This github directory contains our experiments about self-initiated research
project: a method to achieve acceleration by decomposing optimization to two subproblem.

DEVELOPERS
---------------

    Jimmy Lin - jimmylin@utexas.edu
    Ian Yen - a061105@gmail.com


Inventory
--------------

    Makefile                for compilation
    README.md               records basic information

    util.h                  auxiliary facility function
    mat.h                   contains auxiliary functions for mat data structure
    sparseClustering.h      header file for algorithmic script
    sparseClustering.cpp    codebase as our algorithmic script

TODOs
---------------
1. Fix bugs about result inconsistency according to the input data shuffling.
   Possible solution involves in tracing the computational process.

2. The case of indistinguishable cluster assignment (0.50 to one cluster and
   \0.5 to the other one)

3. Compare and improve current performance. 

4. Apply similar approach to other model.

DEVELOPMENT LOGS
---------------

1. [June 30, 2014 by Jimmy] Develop exact line search to accelerate the convergence of
   frank-wolfe algorithm. This modification significantly enhance the
   performance of frank-wolfe algorithm. We compare the iteration number of
   converging first subproblem objective function to 22 with 35(**exact line search**) / 82(**inexact line search**). 
   
2. [June 30, 2014 by Jimmy] We trace the finally converged result and found that the
   element of resulted **W** is not stricly either 0 or 1. That means
   belonging of instances is fuzzier than our expectation. (**why???**) Furthermore, we
   trace the number of cluster centroids varies with regard to the **Overall
   iteration number** as follows:  

   | iter | nCentroid |
   | -----|:---------:|
   | 100  | about 60  |
   | 200  | about 40  |
   | 300  | about 30  |

   (**why???**) However, the synthetic dataset is generated by means of four
   acutal clustered centroid with some noises. 

3. [June 30, 2014 by Jimmy] We refine the mechanism of **managing variable
   dumping works** with macros definition at the very beginning of code script.

4. [July 1, 2014 by Jimmy] Fix up a bug in frank-wolfe algorithm. It resolves
   the problem of glic's dumping double free and corruption problem. It is
   technically math problem. Previously, we did not consider the case of 
   **|| w - z || ^ 2 = 0 **, which causes the crash in frank-wolfe algorithm.
   This condition occurs as convergence indicator of frank-wolfe algorithm,
   saying that there is no further growth for **w** as w^(k+1) = w^(k) + gamma
   \* (s - w^(k)). 

5. [July 1, 2014 by Jimmy] the found optima would vary according to the
   placement ordering of the instances in dataset. I swap the order of a few
   data instances and find that **the converged objective varies!**.  

6. [July 1, 2014 by Jimmy] fix up a minor problem about L2norm. previously,
   have mix up N and D as dimensionality of instances. 
