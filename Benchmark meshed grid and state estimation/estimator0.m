% =======================================
%   Initial state estimation
% =======================================

clear all
tic
global h cons Flows Buses
global Jh Jc 

% Network data input
% ===========================
A=load('branches.txt');
Bus=load('buses.txt');
Shunt=load('shunts.txt');

[m,n]=size(A);
ci=A(1:m,1); 
cj=A(1:m,2);
circ=A(1:m,3);
rgr=A(1:m,4);
xgr=A(1:m,5);
bgro=A(1:m,6);
a=A(1:m,8);

zgr=rgr+xgr*i;  % Branch impedance



for k=1:m;
    ygr(k,1)=1/zgr(k);    % Branch admitance
    ggr(k,1)=real(ygr(k));    
    bgr(k,1)=imag(ygr(k));
end

nbus=length(Bus(:,1)); % Number of network nodes
B=[ci,cj,circ,ggr,bgr,bgro,a];

for i=1:length(Bus(:,1))
     if Bus(i,2)==3 || Bus(i,2)==2
        slack=i;
     end
 end

% Branch power flows - measures
% =====================
Cp=load('Conventional Active Flow.txt');
Cq=load('Conventional Reactive Flow.txt');

[ntp,mtp]=size(Cp);
[ntq,mtq]=size(Cq);
ptm=Cp(:,4:end-1);
qtm=Cq(:,4:end-1);
sigmatp=Cp(:,end);
sigmatq=Cq(:,end);

state1=zeros(2*length(Bus(:,1)),mtp-4);%for saving all the results


% Injected powers - measures
% ==========================
Dp=load('Conventional Pinj.txt');
Dq=load('Conventional Qinj.txt');
[nip,mip]=size(Dp);
[niq,miq]=size(Dq);
pim=Dp(:,2:end-1);
qim=Dq(:,2:end-1);
sigmaip=Dp(:,end);
sigmaiq=Dq(:,end);
%  
% Node voltages - measures
% ========================
E=load('Conventional Voltage.txt');
nu=length(E(:,1));
ciu=E(:,1);
um=E(:,3:end-1); 
sigmau=E(:,end);
% 
% R matrix
% ================

sigma=[sigmaip;sigmaiq;sigmatp;sigmatq;sigmau].^2;
R=diag(sigma);

% 
% State variables initial values
% =====================================
F=load('state0.txt');
ce=F(1:nbus,1);
V=F(1:nbus,2);
teta=F(1:nbus,3);
States=[teta;V];

for time=1:length(ptm(1,:))
tic
%States = [teta;V]; % flat starts
time
% Vector Z calculation
% ======================
Z=[pim(:,time);qim(:,time);ptm(:,time);qtm(:,time);um(:,time)];



% Iterative procedure for state estimation
% ========================================
stop=0;
iter=0;
itermax=50;

while stop==0 && iter<itermax
%
iter=iter+1;
    if iter==itermax
       disp('PROGRAM CAN NOT CONVERGE AFTER 50 ITERATIONS');
    end
    
% Vector h calculation
%========================

evalfunc0(States,B,Bus,Cp,Cq,Dp,Dq,E,Shunt)

% Jacobian Calculation
%========================
jacobianES0(States,B,Bus,Cp,Cq,Dp,Dq,E,Shunt)

 i1=1;
 H=zeros(length(Jh(:,1)),length(States)-1);
 %C=zeros(length(Jc(:,1)),length(States)-1);
  
 for i=1:length(States)
     if i~=slack 
       H(:,i1)=Jh(:,i); 
%       C(:,i1)=Jc(:,i);
       i1=i1+1;
     end
 end
 
% State variables change calculation
% =========================================

delta=inv(H'*inv(R)*H)*H'*inv(R)*(Z-h);

deltax=delta(1:length(States)-1);
i2=1;
for i=1:length(States)
    if i~=slack
        States(i)=States(i)+deltax(i2);
        i2=i2+1;
    end
end

% Convergence test
% ==================
if max(abs(deltax))<1E-6
        stop=1;
end
end % end while
toc
x=[States(1:nbus);States(nbus+1:2*nbus)];
state1(:,time)=States;

%ReportSE(x,Bus,A,Shunt)

J=(Z-h)'*inv(R)*(Z-h) % this is the objective
iter

end

toc