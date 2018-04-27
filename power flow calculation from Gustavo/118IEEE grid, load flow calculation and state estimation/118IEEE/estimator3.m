% =======================================
%   Initial state estimation
% =======================================

clear i  Jh Jc h  cons Z  
tic
global h cons Flows Buses 
global  Jh Jc

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

nbus=length(Bus(:,1));% Number of network nodes
B=[ci,cj,circ,ggr,bgr,bgro,a];

for i=1:length(Bus(:,1))
     if Bus(i,2)==3 || Bus(i,2)==2
        slack=Bus(i,1);
        slackid=i;
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
Pi=Dp(:,2);
Qi=Dq(:,2);
sigmaip=Dp(:,3);
sigmaiq=Dq(:,3);
 
% Bus voltages - measures
% ========================
E=load('Conventional Voltage.txt');
nu=length(E(:,1));
ciu=E(:,1);
u=E(:,3); 
sigmau=E(:,4); 
 

% V phasor- PMU
%========================
F=load('V phasor PMU Measured.txt');
nupmu=length(F(:,1));
Vp=F(:,3);
sigmaup=F(:,5);

i1=1;slackyes=0;
for i=1:length(F(:,1))
    if F(i,2)~=3 && F(i,2)~=2
       cp(i1,1)=F(i,4);
       sigmac(i1,1)=F(i,6);
       cog(i1,1)=vpmuog_(i,4);
       i1=i1+1;
    else
        slackyes=1;
    end
end


% I Phasor PMU
%=========================
K=load('I phasor PMU Measured.txt');
nI=length(K(:,1));
nIc=length(K(:,1));
im=K(:,4);
ic=K(:,5);
sigmaim=K(:,6);
sigmaic=K(:,7);

% Pseudo Measurements
%=========================
S=load('Pseudo Voltage PMU Measured.txt');
ps=S(:,4);
sigps=S(:,5:6);

Ra=zeros(2*nI);
i1=1; 
i2=1;
for i=1:nI
    Ra(i1,i1)=sigps(i2,1); 
    Ra(i1,i1+1)=sigps(i2,2);
    Ra(i1+1,i1)=sigps(i2+1,1); 
    Ra(i1+1,i1+1)=sigps(i2+1,2);
    i1=i1+2;
    i2=i2+2;
end

% R matrix
% ================
sigma=[sigmaip;sigmaiq;sigmatp;sigmatq;sigmau;sigmaup;sigmac];

Rb=eye(length(sigma));   
for i=1:length(sigma)
    Rb(i,i)=sigma(i)^2;
end


R=[Rb zeros(length(Rb),length(Ra));zeros(length(Ra),length(Rb)) Ra];

% Vector Z calculation
% ======================
Z=[Pi;Qi;pt;qt;u;Vp;cp;ps];
 
% 
% State variables initial values
% =====================================
Xo=load('state0.txt');
ce=Xo(1:nbus,1);
V=Xo(1:nbus,2);
teta=Xo(1:nbus,3);

States=[teta;V];
% 

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

evalfunc3(States,B,Bus,Cp,Cq,Dp,Dq,E,F,K,S,Shunt)

% Jacobian Calculation
%========================
jacobianES3(States,B,Bus,Cp,Cq,Dp,Dq,E,F,K,S,Shunt)
% 
if slackyes==1
   H=zeros(length(Jh(:,1)),length(States)-1);
   C=zeros(length(Jc(:,1)),length(States)-1);
else
   H=zeros(length(Jh(:,1)),length(States));
   C=zeros(length(Jc(:,1)),length(States));
end
% % 
if slackyes==1
   i1=1;
   for i=1:length(Jh(1,:))
       if i~=slackid
          H(:,i1)=Jh(:,i); 
          C(:,i1)=Jc(:,i);
          i1=i1+1;
       end
   end
else
   H=Jh;
   C=Jc;
end

%
% State variables change calculation
% =========================================

alpha=max(diag((R)));
delta=inv([alpha*H'*inv(R)*H -C';-C zeros(length(cons))])*[alpha*H'*inv(R)*(Z-h);cons];  

d1=1;
if slackyes==1
   deltax=delta(1:length(States)-1);
else
   deltax=delta(1:length(States));
end

% 
i2=1;
if slackyes==1
   for i=1:length(States)
       if i~=slackid
          States(i)=States(i)+deltax(i2);
          i2=i2+1;
       end
   end
else
    States=States+deltax;
end
% 
% Convergence test
% ==================
if max(abs(deltax))<1E-6
        stop=1;
end

end
%%

Ztrue=[Piog;Qiog;Ptog;Qtog;uog;vog;cog;Sog];
evalfunc3(States,B,Bus,Cp,Cq,Dp,Dq,E,F,K,S,Shunt)
sum1=0;sum2=0;sum3=0;
vecaux=(Ztrue-h).^2./diag(R);
for u=1:length(Z)
     sum1=sum1+(h(u)-Ztrue(u)).^2;
     sum2=sum2+(Z(u)-Ztrue(u)).^2;
     sum3=sum3+vecaux(u);
end

Index1=sum1/sum2;
Index2=sum3;

x=[States(1:nbus);States(nbus+1:2*nbus)];
J=(Z-h)'*inv(R)*(Z-h);
iter;
% 
% % Final Estimation Report
% % ========================
% % Estimation=zeros(length(Bus(:,1)),4);
% % Estimation(:,1:2)=Bus(:,1:2);
% % Estimation(:,3)=States(nbus+1:2*nbus);
% % Estimation(:,4)=States(1:nbus)
% 
% % ReportSE (States,Bus,A,Shunt)
% % dat3=fopen('Flows Final Estimation.txt','w');
% % for k=1:length(A(:,1))*2
% %    fprintf(dat3,'%3i %3i %3i %8.8f %8.8f \n',Flows(k,1:5));
% % end   
% % fclose(dat3);
% % 
% % dat4=fopen('Buses Final Estimation.txt','w');
% % for k=1:length(Bus(:,1))
% %    fprintf(dat4,'%3i %3i %8.6f %8.6f %8.8f %8.8f \n',Buses(k,1:6));
% % end   
% % fclose(dat4);
% % 
% % Estimation=[Bus(:,1:2),States(nbus+1:2*nbus),States(1:nbus)]
% % 
% % dat5=fopen('Bus States.txt','w');
% % for k=1:nbus
% %    fprintf(dat5,'%3i %3i %8.10f %8.10f \n',Estimation(k,1:4));
% % end   
% % fclose(dat5);
% 
 toc