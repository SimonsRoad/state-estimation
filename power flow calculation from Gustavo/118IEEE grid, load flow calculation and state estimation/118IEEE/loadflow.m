function loadflow(Bus,Branch,Gentor,Load,Motor,Shunt)
% This file perform a load flow analysis
% The output of this function is LF
%
global Ym LF GentorLF;
global Pmist Qmist Pmistmax Qmistmax J Ang Vol MotorLF;
%
Ymatrix(Bus,Branch,Shunt)% have to use again in equilibrium program
%
if length(Motor)~=0
    motorlf(Bus,Motor)
    Motor=MotorLF;
end
mismatchlf(Bus,Ym,Gentor,Load,Motor)
%
r=1;k=1;
% 
% NR formulation  
iter=0;
while r==1
jacobianlf(Bus,Ym)
Yg=[Ang;Vol];
Yg=Yg+inv(J)*[Pmist;Qmist];
i1=1;i2=1; % auxiliar
for i=1:length(Bus(:,1)) % assign new values to Bus matrix
    if Bus(i,2)==0 || Bus(i,2)==4
    Bus(i,3)=Yg(i1+length(Ang));
    i1=i1+1;
    end
    if Bus(i,2)~=2 && Bus(i,2)~=3
    Bus(i,4)=Yg(i2);
    i2=i2+1;
    end
end
%
if length(Motor)~=0
    motorlf(Bus,Motor)
    Motor=MotorLF;
end
mismatchlf(Bus,Ym,Gentor,Load,Motor)
if abs(Pmistmax)<=0.00001 && abs(Qmistmax)<=0.00001
    r=0;
    disp('LOAD FLOW CONVERGED');
    k=0;
end
if iter==100
r=0;
disp('NOT CONVERGED AFTER 100 iterations, CHECK LOADFLOW DATA');
end
iter=iter+1;
end % while
%
%Once main calculation has been performed
%
% SLACK BUS CALCULATION
%
infbus=0;
for i=1:length(Bus(:,1))% to identify slack bus
    if Bus(i,2)==3 || Bus(i,2)==2
    m=Bus(i,1);
    mp=i;
    if Bus(i,2)==2
    infbus=1; % Indicates there is an infinite bus
    end
    end
end
% Slack bus power calculation
Pslack=0;Qslack=0;
for f=1:length(Bus(:,1)) 
      Pslack=Pslack+Bus(mp,3)*Bus(f,3)*(real(Ym(mp,f))*cos(Bus(mp,4)-Bus(f,4))+imag(Ym(mp,f))*sin(Bus(mp,4)-Bus(f,4)));
      Qslack=Qslack+Bus(mp,3)*Bus(f,3)*(real(Ym(mp,f))*sin(Bus(mp,4)-Bus(f,4))-imag(Ym(mp,f))*cos(Bus(mp,4)-Bus(f,4)));
end
%
for i=1:length(Load(:,1))
    if Load(i,1)==Bus(mp,1)
    Pslack=Pslack+Load(i,4); 
    Qslack=Qslack+Load(i,5); 
    end
end
%
if length(Motor)~=0
    for i=1:length(Motor(:,1))
       if Motor(i,1)==Bus(mp,1)
       Pslack=Pslack+Motor(i,5);
       Qslack=Qslack+Motor(i,6);
       end
    end
end
Slack=[Bus(mp,1),Pslack,Qslack,infbus];
%
% Q GENERATED FOR ALL PV BUSES
%
d=0;
for i=1:length(Bus(:,1))% to quantify PV buses
    if Bus(i,2)==1
    d=d+1;
    end
end
%
 Qg=zeros(d,2);% Gives reactive powers and buses where they are connected
 ia=1;%aux
for i=1:length(Bus(:,1))% Reactive power calculation of PV buses 
    if Bus(i,2)==1
     Qg(ia,1)=Bus(i,1);
     for z=1:length(Bus(:,1)) 
     Qg(ia,2)=Qg(ia,2)+Bus(i,3)*Bus(z,3)*(real(Ym(i,z))*sin(Bus(i,4)-Bus(z,4))-imag(Ym(i,z))*cos(Bus(i,4)-Bus(z,4)));
     end
     for u=1:length(Load(:,1))
         if Load(u,1)==Qg(ia,1)
          Qg(ia,2)=Qg(ia,2)+Load(u,5);
         end
     end
     %
     if length(Motor)~=0
        for u=1:length(Motor(:,1))
            if Motor(u,1)==Qg(ia,1)
            Qg(ia,2)=Qg(ia,2)+Motor(u,6);
            end
        end
    end
     ia=ia+1;
    end
end 
%
% POWER FLOW IN ALL BRANCHES
%
for f=1:length(Branch(:,1))
    i=Branch(f,1);k1=Branch(f,2);
    for h=1:length(Bus(:,1))
        if Bus(h,1)==i;
           Vi=Bus(h,3);
           ci=Bus(h,4);
        end
        if Bus(h,1)==k1;
           Vk=Bus(h,3);
           ck=Bus(h,4);
        end
    end
    Yik=1/((Branch(f,4))+j*(Branch(f,5)));    
    Gik=real(Yik);
    Bik=imag(Yik);
    circ=Branch(f,3);
    Bii=Branch(f,6)/2;
    aik=Branch(f,8);
    Ptrs1=+Vi^2*1/aik^2*Gik-Vi*Vk*1/aik*(Gik*cos(ci-ck)+Bik*sin(ci-ck)); 
    Qtrs1=-Vi^2*1/aik^2*(Bii+Bik)-Vi*Vk*1/aik*(Gik*sin(ci-ck)-Bik*cos(ci-ck)); 
    Ptrs2=+Vk^2*Gik-Vk*Vi*1/aik*(Gik*cos(ck-ci)+Bik*sin(ck-ci));  
    Qtrs2=-Vk^2*(Bii+Bik)-Vk*Vi*1/aik*(Gik*sin(ck-ci)-Bik*cos(ck-ci)); 
    %
    flows1(f,:)=[i,k1,circ,Ptrs1,Qtrs1];
    flows2(f,:)=[k1,i,circ,Ptrs2,Qtrs2];
    %
    I1=conj((Ptrs1+j*Qtrs1)/(Vi*(cos(ci)+j*sin(ci))));
    I2=conj((Ptrs2+j*Qtrs2)/(Vk*(cos(ck)+j*sin(ck))));
    Icur1rect(f,:)=[i,k1,circ,real(I1),imag(I1)];
    Icur2rect(f,:)=[k1,i,circ,real(I2),imag(I2)];
    Icur1(f,:)=[i,k1,circ,sqrt(real(I1)^2+imag(I1)^2),atan2(imag(I1),real(I1))];
    Icur2(f,:)=[k1,i,circ,sqrt(real(I2)^2+imag(I2)^2),atan2(imag(I2),real(I2))];
   
end
    Flows=[flows1;flows2] % 
    Currents=[Icur1;Icur2];

    Rectangular=[Icur1rect;Icur2rect]
    
%
%Inclusion of gen powers in Generators
%
for i=1:length(Qg(:,1)) % 
    for g=1:length(Gentor(:,1))
       if Qg(i,1)==Gentor(g,1)
       Gentor(g,5)=Qg(i,2);
       end
    end
end
%
if infbus==0
    for g=1:length(Gentor(:,1))
        if Gentor(g,1)==Bus(mp,1)
        Gentor(g,4)=Pslack;
        Gentor(g,5)=Qslack;
        end
    end
end
%
LF=Bus;
GentorLF=Gentor;
MotorLF=Motor;
LF(:,6)=LF(:,3);
%
if k==0
 reportLF(LF,Flows,Currents,GentorLF,Load,Motor,Shunt,Slack); % Provides output of function "loadflow"
end