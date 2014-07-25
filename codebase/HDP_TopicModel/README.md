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
    
    HDP.cpp                 HDP codes
    HDP.h                   HDP header file

    exSparseMat.h           extensible sparse matrix API

TODOs
---------------
1. test availability of esmat API functions
2. modify frank\_wolfe\_solver by utilizing esmat
3. modify frank\_wolfe\_solver to respect the new loss function (dummy only)
4. modify group\_lasso\_solver by using esmat

DEVELOPMENT LOGS
---------------

1. [Mon Jul 21 22:19:58 2014 by Jimmy] repository setup

2. [Mon Jul 21 22:20:22 2014 by Jimmy] modifies HDP function and frank\_wolfe

3. [Wed Jul 23 23:00:00 2014 by Jimmy and Ian] identify clearly definition of loss and three regularization terms

4. [Wed Jul 23 23:00:00 2014 by Jimmy and Ian] design data structure esmat for
   w, z and y

5. [Fri Jul 25 02:06:40 2014 by Jimmy] implement initial version of esmat and
   do corresponding modification on HDP.cpp based on sparseClustering.cpp
