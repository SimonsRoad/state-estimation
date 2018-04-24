function [ Vf, thetaf, a, b, zf, Tf, H] = forcast (a0, b0, alpha, beta, Vf, thetaf, Ve, thetae, Pe, Q,...
    topo, Y_bus, ind_meas, N_meas, H_decoupled, H_sparse )

%this function is for forcasting the state vectors at time k+1
%exmaple 14 buses
%output Vf(14*1) and thetaf(13*1) are forcasted values at time k+1
%a(27*1) and b(27*1) are values at time k
%output Pf is the error covariance of forcasted state vector x at time k+1
%zf and Tf are forcasted measurement and corresponding covariance matrix at time k
%a0 and b0 are values at time k-1
%Vf and thetaf are forcasted values at time k, dimension 27*1
%Ve and thetae are estimated value at time k, dimension 28*1
%imput Pe is the error covariance of estimated state vector x at time k
%imput Q is the error covariance of modeling uncertainty at time k


xf=[thetaf;Vf];%dimension 27*1
xe=[thetae(2:topo.nBus);Ve];%dimension 27*1, not include reference phase
a=alpha*xe+(1-alpha)*xf;%dimension 27*1
b=beta*(a-a0)+(1-beta)*b0;%dimension 27*1
F=alpha*(1+beta)*eye(topo.nBus*2-1);%%dimension 27*27, reference theta not included
G=(1+beta)*(1-alpha)*xf-beta*a0+(1-beta)*b0;%dimension 27*1
x=F*xe+G;%dimension 27*1
thetaf=x(1:topo.nBus-1);%%dimension 13*1, not include the reference theta
Vf=x(topo.nBus:topo.nBus*2-1);%dimension 14*1
Mf=F*Pe*F.'+Q;%dimension 27*27, M(k+1) forcast

%calculate H and h by using forcast value of V and theta at k+1
[ H ] = f_measJac_H_v2017( Vf, [0;thetaf], Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
[ h ] = f_measFunc_h_v2017( Vf, [0;thetaf], Y_bus, topo, ind_meas, N_meas);
zf=h;
Tf=H*Mf*H.';









