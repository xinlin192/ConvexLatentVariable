Dataset colon-cancer
=========================


Algorithm: PAM
-------------------------

    ./PAM ../../dataset/colon-cancer/colon-cancer 10 10000 2
    min_objective: 1567.96

Algorithm: DP\_MEDOIDS
--------------------------

    ./DP_MEDOIDS ../../dataset/colon-cancer/colon-cancer 1 1000 80
    min_objective: 1567.96

Algorithm: DP\_Cvx\_MEDOIDS
------------

    ./cvx_clustering ../../dataset/colon-cancer/colon-cancer 100 3000 100
    loss: 1564.75
