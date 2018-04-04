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
n=topo.nBus*2-1;

%initial start point
%for prediction
Pe=0.1*eye(topo.nBus*2-1,topo.nBus*2-1);
Q=ones(topo.nBus*2-1,topo.nBus*2-1);
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c-1,1);
Vf=ones(c,1);
thetaf=zeros(c-1,1);
a=zeros(c*2-1,1);
b=zeros(c*2-1,1);
Xf=ones(2*topo.nBus-1,2*n+1);
%for sigma points
alpha_s=0.1;%selected from 10^(-4) to 1
beta_s=2;

if n<=3
    k=3-n;
else
    k=0;
end
lamda=alpha_s^2*(n+k)-n;
Wm=1/(2*(n+lamda))*ones(2*n+1,1);
Wc=Wm;
Wm(1)=lamda/(n+lamda);
Wc(1)=lamda/(n+lamda)+1-alpha_s^2+beta_s;
Wc=diag(Wc);
for o=1:10

    %SOLTUION ALGORITHM STARTS HERE
    
    %build sigma points equation (17) on PSL paper
    
    x0=[thetae;Ve];
    X0=repmat(x0,1,2*n+1);
    c=n+lamda;
    zerox(1:n,1)=0;
    Xe=X0+sqrt(c)*[zerox, sqrt(Pe), -sqrt(Pe)];
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
    [ Xf, a, b] = forcast (a, b, alpha, beta, Xf, Xe, Q,...
     topo);
    
    xf=Xf*Wm;%transform forcasted state vectors from sigma form Xf to xf
    xfp=repmat(xf,1,2*n+1);%for Xf minis 2*n+1 columns of xf
    Pf=(Xf-xfp)*Wc*(Xf-xfp).'+Q;%covariance of the forcasted state vectors
    

    %TASK 2: state correction
    X1=repmat(xf,1,2*n+1);
    Xep=X1+sqrt(c)*[zerox, sqrt(Pf), -sqrt(Pf)];%new sigma points in step 3
    Xep_theta=Xep(1:topo.nBus-1,1);
    Xep_v=Xep(topo.nBus:n,1);    
    Y0=f_measFunc_h_v2017( Xep_v, [0;Xep_theta], Y_bus, topo, ind_meas, N_meas);
    Y=Y0;
    for m=2:2*n+1
        Xep_theta=Xep(1:topo.nBus-1,m);
        Xep_v=Xep(topo.nBus:n,m);
        Yk=f_measFunc_h_v2017( Xep_v, [0;Xep_theta], Y_bus, topo, ind_meas, N_meas);
        Y=[Y,Yk];        
    end
    y=Y*Wm;%miu in PSL paper
    yp=repmat(y,1,2*n+1);%for Xf minis 2*n+1 columns of xf
    S=(Y-yp)*Wc*(Y-yp).'+R;
    C=(Xep-xfp)*Wc*(Y-yp).';
    K=C*inv(S);
    xe=xf+K*(z-y);
    Pe=Pf-K*S*K.';
    thetae=xe(1:topo.nBus-1,1);%add reference theta inside
    Ve=xe(topo.nBus:topo.nBus*2-1,1);
end

thetae=[0;xe(1:topo.nBus-1,1)];%add reference theta inside
Ve=xe(topo.nBus:topo.nBus*2-1,1);







