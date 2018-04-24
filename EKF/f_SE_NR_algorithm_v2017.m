function [ Ve, thetae, xe] = f_SE_NR_algorithm_v2017 (topo, z, zf, H, R, xf)

%SOLTUION ALGORITHM STARTS HERE

v=z-zf;
Pe = inv(H.'/R*H);%state vector covariance matrix at time k+1 
K = Pe*H.'/R;

xe = xf+K*v;%xe is 27*1 here
thetae=[0;xe(1:topo.nBus-1,1)];%add reference theta inside
Ve=xe(topo.nBus:topo.nBus*2-1,1);

