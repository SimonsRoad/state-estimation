% =======================================
%   Initial state estimation
% =======================================

clear i Jh Jc h  cons Z
tic
global h cons 
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
pt=Cp(:,4);
qt=Cq(:,4);
sigmatp=Cp(:,5);
sigmatq=Cq(:,5);

% Injected powers - measures
% ==========================
Dp=load('Conventional Pinj.txt');
Dq=load('Conventional Qinj.txt');
[nip,mip]=size(Dp);
[niq,miq]=size(Dq);
pi=Dp(:,2);
qi=Dq(:,2);
sigmaip=Dp(:,3);
sigmaiq=Dq(:,3);
 
% Node voltages - measures
% ========================
E=load('Conventional Voltage.txt');
nu=length(E(:,1));
ciu=E(:,1);
u=E(:,3); 
sigmau=E(:,4);

% R matrix
% ================

sigma=[sigmaip;sigmaiq;sigmatp;sigmatq;sigmau];

R=eye(length(sigma));   
for i=1:length(sigma)
    R(i,i)=sigma(i)^2;
end

% Vector Z calculation
% ======================
Z=[pi;qi;pt;qt;u];
% 
% State variables initial values
% =====================================
F=load('state0.txt');
ce=F(1:nbus,1);
V=F(1:nbus,2);
teta=F(1:nbus,3);
States=[teta;V];


% Iterative procedure for state estimation
% ========================================
stop=0;
iter=0;
itermax=50;

while stop==0 && iter<itermax
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
 C=zeros(length(Jc(:,1)),length(States)-1);
  
 for i=1:length(States)
     if i~=slack 
       H(:,i1)=Jh(:,i); 
       C(:,i1)=Jc(:,i);
       i1=i1+1;
     end
 end
 
% State variables change calculation
% =========================================

alpha=max(diag((R)));

delta=inv([alpha*H'*inv(R)*H -C';-C zeros(length(cons))])*[alpha*H'*inv(R)*(Z-h);cons];  

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

end


%%
% Ztrue=[Piog;Qiog;Ptog;Qtog;uog];
% evalfunc0(States,B,Bus,Cp,Cq,Dp,Dq,E,Shunt)
% sum1=0;sum2=0;sum3=0;
% vecaux=(Ztrue-h).^2./sigma.^2;
% for u=1:length(Z)
%      sum1=sum1+(h(u)-Ztrue(u)).^2;
%      sum2=sum2+(Z(u)-Ztrue(u)).^2;
%      sum3=sum3+vecaux(u);
% end
% 
% Index1=sum1/sum2;
% Index2=sum3;

x=[States(1:nbus);States(nbus+1:2*nbus)];

ReportSE(x,Bus,A,Shunt)

J=(Z-h)'*inv(R)*(Z-h);
iter;
toc