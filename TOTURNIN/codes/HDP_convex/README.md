Convex hdp medoids 
=======================

Code Inventory 
--------------

    Makefile                for compilation
    README.md               records basic information

    HDP_convex.h      header file for algorithmic script
    HDP_convex.cpp    codebase as our algorithmic script

Executive Commands
-----------------
The command to get HDP\_convex program run:

    cvx_hdp_medoids [word_dataFile] [lambda_global] [lambda_local] [FW_MAX_ITER] [ADMM_MAX_ITER]

e.g.
    ./HDP_convex ../../datasets/wholesale.scale.libsvm 1.0 1.0 30 100000

