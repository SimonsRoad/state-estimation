clc;clear;
%LOADING THE DATA
load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=4;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid




%SOLTUION ALGORITHM STARTS HERE
%getting Jacobian sparsity pattern (done only once)
[topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2017(topo, meas);

%obtaining the admittance matrix in the compressed row form (done only once)
Y_bus = f_Y_bus_compressed_v2017(topo);



%constructing vector z, which is a vector of available measurements (hint: use structure ind_meas)
%NOTE: different types of measurements should be in the same order as in the handout
%STUDENT CODE 1
%%combine derivations and measurements into one matrix Com
Com=[meas.BusV,meas.std_BusV;
     meas.BusInjP,meas.std_BusInjP;
     meas.BusInjQ,meas.std_BusInjQ;
     meas.LineP_FT,meas.std_LineP_FT;
     meas.LineP_TF,meas.std_LineP_TF;
     meas.LineQ_FT,meas.std_LineQ_FT;
     meas.LineQ_TF,meas.std_LineQ_TF;
    ];
Com = Com(all(~isnan(Com),2),:); % for nan - rows
z = Com(:,1);



%constructing the matrix of weights W and its square root Wsqrt using standard deviations provided in structure 'meas'
%NOTE: matrices W and Wsqrt must be constructed in sparse format!!!
%STUDENT CODE 2

W=diag(Com(:,2).^(-2));
W=sparse(W);
Wsqrt=diag(Com(:,2).^(-1));

%choosing the initial point (use flat start)
%STUDENT CODE 3
a=size(topo.busNumbers,1);
V=ones(a,1);
theta=zeros(a,1);


%TASK 1: SOLVING SE PROBLEM WITH GN ALGORITHM
%STUDENT CODE 4
%NOTE: use the following function (the description of its inputs and
%outputs can be found inside the function)
[ V, theta, eps_all, time, convergence] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
            ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );



%TASK 2: COMPUTING CONDITION NUMBERS AND DENSITY FACTORS OF MATRICES
%STUDENT CODE 5




%TASK 3: COMPARING PERFORMANCE OF FULL/DECOUPLED GN
%STUDENT CODE 6
%NOTE: use the following function (the description of its inputs and
%outputs can be found inside the function)
%[ V, theta, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
%            ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );






