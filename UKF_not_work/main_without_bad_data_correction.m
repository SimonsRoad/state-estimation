clear all
%LOADING THE DATA
%load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;
%load Dataset_A.mat;
load Measurement.mat

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
alpha=0.775;%smoothing parameter
beta=0.5;%smoothing parameter
n=topo.nBus*2-1;

%initial start point
%for prediction
Pe=0.1*eye(topo.nBus*2-1,topo.nBus*2-1);
Q=0.1*eye(topo.nBus*2-1,topo.nBus*2-1);
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c-1,1);
a=zeros(c*2-1,1);
b=zeros(c*2-1,1);
Xf=ones(2*topo.nBus-1,2*n+1);
%for sigma points
alpha_s=0.1;%for UT, selected from 10^(-4) to 1
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

d=size(meas.BusV,2);%getting the total number of measurements
meas_new=meas;
for i=1:d
    meas.BusInjP=meas_new.BusInjP(:,i);
    meas.BusInjQ=meas_new.BusInjQ(:,i);
    meas.BusV=meas_new.BusV(:,i);
    meas.LineP_FT=meas_new.LineP_FT(:,i);
    meas.LineP_TF=meas_new.LineP_TF(:,i);
    meas.LineQ_FT=meas_new.LineQ_FT(:,i);
    meas.LineQ_TF=meas_new.LineQ_TF(:,i);
    meas.std_BusInjP=meas_new.std_BusInjP(:,i);
    meas.std_BusInjQ=meas_new.std_BusInjQ(:,i);
    meas.std_BusV=meas_new.std_BusV(:,i);
    meas.std_LineP_FT=meas_new.std_LineP_FT(:,i);
    meas.std_LineP_TF=meas_new.std_LineP_TF(:,i);
    meas.std_LineQ_FT=meas_new.std_LineQ_FT(:,i);
    meas.std_LineQ_TF=meas_new.std_LineQ_TF(:,i); 
    %getting measurement
    [topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2017(topo, meas);
    %obtaining the admittance matrix in the compressed row form
    Y_bus = f_Y_bus_compressed_v2017(topo);

    %constructing vector z, which is a vector of available measurements
    %NOTE: different types of measurements should be in the same order as in the handout
    %STUDENT CODE 1
    z=[meas.BusV(ind_meas.ind_meas_V); meas.BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.BusInjQ(ind_meas.ind_meas_Qinj); meas.LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.LineP_TF(ind_meas.ind_meas_Pji); meas.LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.LineQ_TF(ind_meas.ind_meas_Qji)];
    N_total=length(z);

    %constructing the matrix of measurements covariance matrix R and state vector covariance matrix Pe using standard deviations provided in structure 'meas'
    %NOTE: matrices R must be constructed in sparse format!!!
    %STUDENT CODE 2
    R_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
    R=sparse(1:N_total,1:N_total,((R_diag).^2),N_total,N_total);


    %SOLTUION ALGORITHM STARTS HERE
    
    %build sigma points equation (17) on PSL paper
    
    x0=[thetae;Ve];
    X0=repmat(x0,1,2*n+1);
    c=n+lamda;
    zerox(1:n,1)=0;
    L1=chol(Pe,'lower');
    Xe=X0+sqrt(c)*[zerox, L1, -L1];
    %TASK 1: forcasted state vectors and measurements
    %Vf 14*1, thetaf 13*1
    [ Xf, a, b] = forcast (a, b, alpha, beta, Xf, Xe, Q,...
     topo);
    
    xf=Xf*Wm;%transform forcasted state vectors from sigma form Xf to xf
    xfp=repmat(xf,1,2*n+1);%for Xf minis 2*n+1 columns of xf
    Pf=(Xf-xfp)*Wc*(Xf-xfp)'+Q;%covariance of the forcasted state vectors
    

    %TASK 2: state correction
    X1=repmat(xf,1,2*n+1);%sigma points for predicted state mean vector
    L2=chol(Pf,'lower');
    Xep=X1+sqrt(c)*[zerox, L2, -L2];%new sigma points in step 3
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
    S=(Y-yp)*Wc*(Y-yp)'+R;
    C=(Xep-xfp)*Wc*(Y-yp)';
    K=C/S;
    xe=xf+K*(z-y);
    Pe=Pf-K*S*K';
    thetae=xe(1:topo.nBus-1,1);%add reference theta inside
    Ve=xe(topo.nBus:topo.nBus*2-1,1);  
end









