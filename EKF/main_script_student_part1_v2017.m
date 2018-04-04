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
alpha=0.775;%smoothing parameter
beta=0.1;%smoothing parameter


%initial start point
Pe=0.1*eye(topo.nBus*2-1,topo.nBus*2-1);
Q=zeros(topo.nBus*2-1,topo.nBus*2-1);
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c,1);
Vf=ones(c,1);
thetaf=zeros(c-1,1);
a=zeros(c*2-1,1);
b=zeros(c*2-1,1);

for i=1:10
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



%constructing the matrix of measurements covariance matrix R and state vector covariance matrix Pe using standard deviations provided in structure 'meas'
%NOTE: matrices R must be constructed in sparse format!!!
%STUDENT CODE 2
R=diag(Com(:,2).^(2));
R=sparse(R);



%TASK 1: forcasted state vectors and measurements
%Vf 14*1, thetaf 13*1
[ Vf, thetaf, a, b, zf, Tf, H] = forcast (a, b, alpha, beta, Vf, thetaf, Ve, thetae, Pe, Q,...
    topo, Y_bus, ind_meas, N_meas, H_decoupled, H_sparse );


%TASK 2: state correction

Pe = inv(H.'/R*H);%state vector covariance matrix at time k+1 
K = Pe*H.'/R;
v = z-zf;

xf = [thetaf;Vf];
xe = xf+K*v;%xe is 27*1 here
thetae=[0;xe(1:topo.nBus-1,1)];%add reference theta inside
Ve=xe(topo.nBus:topo.nBus*2-1,1);

end


%TASK 3: COMPARING PERFORMANCE OF FULL/DECOUPLED GN
%STUDENT CODE 6
%NOTE: use the following function (the description of its inputs and
%outputs can be found inside the function)
%[ V, theta, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
%            ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );






