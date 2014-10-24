Sparse convex optimization for Clustering
======================================

This github directory contains our experiments about self-initiated research
project: a method to achieve acceleration by decomposing optimization to two subproblem.

Inventory of relevant scripts
--------------------------------------

    Makefile                for compilation
    README.md               records basic information

    util.h                  auxiliary facility function
    mat.h                   contains auxiliary functions for mat data structure

    DP_convex.h             header file for algorithmic script
    DP_convex.cpp           codebase as our algorithmic script

Executive Commands
-------------------------------------
Command to be executed for DP\_convex program:

    DP_convex [dataFile] [fw_max_iter] [ADMM_max_iter] [lambda]

e.g. 

    ./DP_convex ../../datasets/iris/iris.scale.1 30 3000 2


Recommended Parameters 
-------------------------------------
Now we give a few selected lambdas that yield integer solution for DP-convex.

Iris: 

    lambda=4~10, K=3
    lambda=12~50, K=2
    lambda=60~, K=1

Glass:

    lambda=4~5, K=6
    lambda=5.5~12, K=4
    lambda=18~, K=1

Wine: 

    lambda=9~10, K=4
    lambda=12~30, K=3
    lambda=35~90, K=2
    lambda=100~, K=1
