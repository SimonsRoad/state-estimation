function [ Xf, a, b] = forcast (a0, b0, alpha, beta, Xf, Xe, Q,...
    topo)

%this function is for forcasting the state vectors at time k+1
%exmaple 14 buses
%a(27*55) and b(27*55) are values at time k
%a0 and b0 are values at time k-1
%imput Q is the error covariance of modeling uncertainty at time k
%Xf dimension 27*(27*2+1)
%Xe dimension 27*55, not include reference phase

a=alpha*Xe+(1-alpha)*Xf;%dimension 27*55
b=beta*(a-a0)+(1-beta)*b0;%dimension 27*55
F=alpha*(1+beta)*eye(topo.nBus*2-1);%%dimension 27*27, reference theta not included
G=(1+beta)*(1-alpha)*Xf-beta*a0+(1-beta)*b0;%dimension 27*55
Xf=F*Xe+G;%dimension 27*55
%Mf=F*Pe*F.'+Q;%dimension 27*27, M(k+1) forcast











