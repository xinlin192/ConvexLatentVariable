Dataset 
=========================

Small Scale
-----------
1. colon-cancer N=60 D=2000(high) C=2
2. iris N=150 D=4(low) C=3
3. glass N=214 D=8 C=6
4. vehicle N=846 D=18 C=4


Medium scale
----------
1. DNA N=2000 D=180 C=3 
2. segment N=2310 D=19(low) C=7
3. sector N=6412 D=55197(high) C=105(high)

Large scale
----------
1. Epsilon N=400,000 D=2000(high) C=2
2. covtype N=581,012 D=54(low) C=7
3. shuttle N=43,500 D=9 C=7(low)


       Dataset & CVX_DP(MEDOIDS) & CVX_DP(MEANS) & DP_MEDOIDS & DP_MEAN 
    ---------------------------------------------------------------------
    Iris(N=150,D=4,L=9) & 41.8408(K=3) & 41.1715(K=3) & 93.5825(K=1) & 91.2763(K=1)
    Iris(N=150,D=4,L=3) & 23.8408(K=3) & 23.1715(K=3) & 31.1684(K=2) & 30.2874(K=2)
    Glass(N=214,D=8,L=9) & 88.3696(K=4) & 83.4829(K=4) & 96.7075(K=2) & 86.3317(K=2)
    Wine(N=178,D=13,L=15) & 164.614(K=3) & 144.466(K=4) & 255.018(K=1) & 206.199(K=1)
    Breast-cancer(N=683,D=10,L=9) & 145.275(K=4) & 127.896(K=4) & 458.584(K=6) & 380.842(K=6)



