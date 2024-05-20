This code contains the algorithm described in the paper: 
Ye Liu, Xuelei Lin, Yejia Chen, and Reynold Cheng, Multi-order Graph Clustering with Adaptive Node-level Weight Learning, submitted. 

The main file is Clusterwithisolated.m that multi-order graph clustering model (MOGC) to integrate multiple motifs and edge information 
in order to resolve fragmentation issue and achieve more accurate partition results simultaneously; namely, it solves
min_{Lambda,U} trace(U^TL_fU)+alpha*|Lambda|_F^2, s.t., U^TD_fU=I, Lambda 1_m=1_n, Lambda>=0

To try the code, you can run: 
1. main.m  that runs the experiments on the football data set from section 7.4 described in the paper above.  