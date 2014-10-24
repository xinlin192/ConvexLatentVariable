Sparse convex optimization for Clustering
=======================

This github directory contains our experiments about self-initiated research
project: a method to achieve acceleration by decomposing optimization to two subproblem.

Inventory
--------------

    Makefile                for compilation
    README.md               records basic information

    util.h                  auxiliary facility function
    mat.h                   contains auxiliary functions for mat data structure
    sparseClustering.h      header file for algorithmic script
    sparseClustering.cpp    codebase as our algorithmic script

Executive Commands
-----------------

    lambda = 10, 3 clusters
    ./cvx_clustering ../../dataset/iris/iris.scale.1 30 3000 10 20
