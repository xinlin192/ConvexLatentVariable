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

Executive Commands
---------------
    
    ./PAM ../../dataset/wine.scale 13 3

DEVELOPMENT LOGS
---------------

1. [Tue Aug 19 21:54:54 2014 by Jimmy] take four hours to develop PAM
   algorithm for clustering.

1. [Tue Aug 19 22:59:45 2014 by Jimmy] take one hour to debug.
