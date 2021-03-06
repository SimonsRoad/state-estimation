function evalfuncHSE(States,B,Bus,Cp,Cq,Dp,Dq,E,F,K,Shunt)
global h cons

clear c V;

nbranch=length(B(:,1));
nbus=length(Bus(:,1));
nflowsp=length(Cp(:,1));
nflowsq=length(Cq(:,1));
npinj=length(Dp(:,1));
nqinj=length(Dq(:,1));
nvol=length(E(:,1));
nvolp=length(F(:,1));
nang=length(F(:,1));
nI=length(K(:,1))/2;



c=States(1:nbus);
V=States(nbus+1:2*nbus);

Gik=zeros(length(Bus(:,1)));
Bik=zeros(length(Bus(:,1)));

for i=1:nbus
    for k=1:nbranch
        if B(k,1)==Bus(i,1)           
           Gik(i,i)=Gik(i,i)+B(k,4)/B(k,7)^2;
           Bik(i,i)=Bik(i,i)+B(k,5)/B(k,7)^2+B(k,6)/2;
        end
         if B(k,2)==Bus(i,1)
           Gik(i,i)=Gik(i,i)+B(k,4);
           Bik(i,i)=Bik(i,i)+B(k,5)+B(k,6)/2;
        end
    end
end

for r=1:nbranch
    for t=1:nbus
        if Bus(t,1)==B(r,1)
           i=t;
        end
        if Bus(t,1)==B(r,2)
           k=t; 
        end
    end
    aik=B(r,7);
    Gik(i,k)=Gik(i,k)-B(r,4)/aik;
    Bik(i,k)=Bik(i,k)-B(r,5)/aik;
    Gik(k,i)=Gik(k,i)-B(r,4)/aik;
    Bik(k,i)=Bik(k,i)-B(r,5)/aik;
end

if isempty(Shunt)==0
    for i=1:nbus
        for r=1:length(Shunt(:,1))
            if Shunt(r,1)==Bus(i,1)
               Gik(i,i)=Gik(i,i)+Shunt(r,3);
               Bik(i,i)=Bik(i,i)+Shunt(r,4);
            end
        end
    end
end

Pinj=zeros(npinj,1);Qinj=zeros(nqinj,1);

% Injected powers - measures
% =====================
for r=1:npinj
    for i=1:nbus
        if Bus(i,1)==Dp(r,1)
           for k=1:nbus
               Pinj(r)=Pinj(r)+V(i)*V(k)*(Gik(i,k)*cos(c(i)-c(k))+Bik(i,k)*sin(c(i)-c(k)));
           end
        end
    end
end

for r=1:nqinj
    for i=1:nbus
        if Bus(i,1)==Dq(r,1)
           for k=1:nbus
            Qinj(r)=Qinj(r)+V(i)*V(k)*(Gik(i,k)*sin(c(i)-c(k))-Bik(i,k)*cos(c(i)-c(k))); 
           end
        end
    end
end


% Branch power flows - measures
% =====================

Ptrs=zeros(nflowsp,1);Qtrs=zeros(nflowsq,1);

for f=1:nflowsp
    for r=1:nbus
        if Bus(r,1)==Cp(f,1)
           i=r;
        end
        if Bus(r,1)==Cp(f,2)
           k=r;
        end
    end
   for s=1:nbranch
       if  Cp(f,1)==B(s,1) && Cp(f,2)==B(s,2)
           if B(s,3)==Cp(f,3)
           gik=B(s,4);
           bik=B(s,5);
           Bii=B(s,6)/2;
           aik=B(s,7);
           aki=1;
           end
       end
       if Cp(f,1)==B(s,2) && Cp(f,2)==B(s,1)
           if B(s,3)==Cp(f,3)
           gik=B(s,4);
           bik=B(s,5);
           Bii=B(s,6)/2;
           aik=1;
           aki=B(s,7);
           end
       end
    end
    Ptrs(f)=+V(i)^2*1/aik^2*gik-V(i)*V(k)*1/aik*1/aki*(gik*cos(c(i)-c(k))+bik*sin(c(i)-c(k)));
end

for f=1:nflowsq
    for r=1:nbus
        if Bus(r,1)==Cq(f,1)
           i=r;
        end
        if Bus(r,1)==Cq(f,2)
           k=r;
        end
    end
 for s=1:nbranch
       if  Cq(f,1)==B(s,1) && Cq(f,2)==B(s,2)
           if B(s,3)==Cq(f,3)
           gik=B(s,4);
           bik=B(s,5);
           Bii=B(s,6)/2;
           aik=B(s,7);
           aki=1;
           end
       end
       if Cq(f,1)==B(s,2) && Cq(f,2)==B(s,1)
           if B(s,3)==Cq(f,3)
           gik=B(s,4);
           bik=B(s,5);
           Bii=B(s,6)/2;
           aik=1;
           aki=B(s,7);
           end
       end
 end
    Qtrs(f)=-V(i)^2*1/aik^2*(Bii+bik)-V(i)*V(k)*1/aik*1/aki*(+gik*sin(c(i)-c(k))-bik*cos(c(i)-c(k)));
end


% Node voltages - measures
% =====================
fV=zeros(nvol,1);
for i=1:nvol
    for r=1:nbus    
       if E(i,1)==Bus(r,1)
       fV(i)=V(r); 
       end
    end
end


% V phasor - PMU 
%=====================
pV=zeros(nvolp,1);
cd=zeros(nang,1);
yes=0;

for i=1:length(F(:,1))
    if F(i,2)==2 || F(i,2)==3
        yes=1;
    end
end

if yes==1
    cd=zeros(nang-1,1);
end

i1=1;
for i=1:nvolp
    for r=1:nbus
        if F(i,1)==Bus(r,1)
           pV(i)=V(r);
           if F(i,2)~=2 && F(i,2)~=3
               cd(i1)=c(r);
               i1=i1+1;
            end
        end
    end
end


% I phasor - PMU
%===========================

imc=zeros(2*nI,1);


i1=1;
for f=1:nI
     for r=1:nbus
        if K(i1,1)==Bus(r,1)
           i=r;
        end
        if K(i1,2)==Bus(r,1)
           k=r; 
        end
    end
    
    for s=1:nbranch
        if K(i1,1)==B(s,1) && K(i1,2)==B(s,2)
           if K(i1,3)==B(s,3)
              aik=B(s,7);
              gik=1/aik*B(s,4);
              bik=1/aik*B(s,5);
              Bii=B(s,6)/2+(1/aik-1)*1/aik*B(s,5);
              Gii=(1/aik-1)*1/aik*B(s,4);
           end
        end
        if K(i1,2)==B(s,1) && K(i1,1)==B(s,2)
           if K(i1,3)==B(s,3)
              aik=B(s,7);
              gik=1/aik*B(s,4);
              bik=1/aik*B(s,5);
              Bii=B(s,6)/2+(1-1/aik)*B(s,5);
              Gii=(1-1/aik)*B(s,4);
           end
        end   
    end
   k1=gik*(V(i)*cos(c(i))-V(k)*cos(c(k)))-bik*(V(i)*sin(c(i))-V(k)*sin(c(k)))-Bii*V(i)*sin(c(i))+Gii*V(i)*cos(c(i));
   k2=gik*(V(i)*sin(c(i))-V(k)*sin(c(k)))+bik*(V(i)*cos(c(i))-V(k)*cos(c(k)))+Bii*V(i)*cos(c(i))+Gii*V(i)*sin(c(i));
   imc(i1)=k1;
   imc(i1+1)=k2;
   i1=i1+2;
end

h=[Pinj;Qinj;Ptrs;Qtrs;fV;pV;cd;imc];


% Null Sinj Buses
% ========================

Busnull=[];
i1=1;
for i=1:nbus
    if Bus(i,2)==4
       Busnull(i1,1)=i; 
       i1=i1+1;
    end
end

if isempty(Busnull)~=1
    
nnull=length(Busnull);
Pnull=zeros(nnull,1);Qnull=zeros(nnull,1);

for r=1:nnull
     for k=1:nbus
         i=Busnull(r);
         Pnull(r)=Pnull(r)+V(i)*V(k)*(Gik(i,k)*cos(c(i)-c(k))+Bik(i,k)*sin(c(i)-c(k)));
     end        
end

for r=1:nnull
    for k=1:nbus
        i=Busnull(r);
        Qnull(r)=Qnull(r)+V(i)*V(k)*(Gik(i,k)*sin(c(i)-c(k))-Bik(i,k)*cos(c(i)-c(k))); 
    end    
end

cons=[Pnull;Qnull];

end
