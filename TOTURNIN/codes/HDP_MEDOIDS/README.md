DP\_MEANS code
=======================

Inventory
--------------

    Makefile                for compilation
    README.md               basic information record

    util.h                  auxiliary facility function
    mat.h                   contains auxiliary functions for mat data structure
    HDP_MEANS.h             header file for algorithmic script
    HDP_MEANS.cpp           codebase as our algorithmic script

Executive Commands
---------------
The execute command

    HDP_MEDOIDS [dataFile] [nRuns] [lambda_global] [lambda_local]

e.g.

    ./HDP_MEDOIDS ../../datasets/wholesale.scale.libsvm 1 1.0 1.0

