function ReportSE (x,Bus,A,Shunt)
global Flows Buses

nbus=length(Bus(:,1));
V=x(nbus+1:2*nbus);
c=x(1:nbus);

% Injected Power Calculation
%===========================

Pi=zeros(nbus,1);
Qi=zeros(nbus,1);

[m,n]=size(A);
ci=A(1:m,1);
cj=A(1:m,2);
circ=A(1:m,3);
rgr=A(1:m,4);
xgr=A(1:m,5);
bgro=A(1:m,6);
a=A(1:m,8);

zgr=rgr+xgr*j;  % Branch impedance

for k=1:m;
    ygr(k,1)=1/zgr(k);    % Branch admitance
    ggr(k,1)=real(ygr(k));    
    bgr(k,1)=imag(ygr(k));
end

B=[ci,cj,circ,ggr,bgr,bgro,a];

Gik=zeros(length(Bus(:,1)));
Bik=zeros(length(Bus(:,1)));

for i=1:length(Bus(:,1))
    for k=1:length(B(:,1))
        if B(k,1)==Bus(i,1)           
           Gik(i,i)=Gik(i,i)+B(k,4)/B(k,7)^2;
           Bik(i,i)=Bik(i,i)+(B(k,5)+B(k,6)/2)/B(k,7)^2;
        end
         if B(k,2)==Bus(i,1)
           Gik(i,i)=Gik(i,i)+B(k,4);
           Bik(i,i)=Bik(i,i)+B(k,5)+B(k,6)/2;
        end
    end
end

for r=1:length(B(:,1))
    for m=1:length(Bus(:,1))
        if Bus(m,1)==B(r,1)
           i=m;
        end
        if Bus(m,1)==B(r,2)
           k=m; 
        end
    end
    aik=B(r,7);
    Gik(i,k)=Gik(i,k)-B(r,4)/aik;
    Bik(i,k)=Bik(i,k)-B(r,5)/aik;
    Gik(k,i)=Gik(k,i)-B(r,4)/aik;
    Bik(k,i)=Bik(k,i)-B(r,5)/aik;
end

if isempty(Shunt)==0
   for r=1:length(Shunt(:,1))
       for m=1:length(Bus(:,1))
           if Bus(m,1)==Shunt(r,1)
              i=m;
           end
       end
       Gik(i,i)=Gik(i,i)+Shunt(r,3);
       Bik(i,i)=Bik(i,i)+Shunt(r,4);
   end
end

for i=1:nbus
    for k=1:nbus
        Pi(i)=Pi(i)+V(i)*V(k)*(Gik(i,k)*cos(c(i)-c(k))+Bik(i,k)*sin(c(i)-c(k)));
        Qi(i)=Qi(i)+V(i)*V(k)*(Gik(i,k)*sin(c(i)-c(k))-Bik(i,k)*cos(c(i)-c(k)));
    end  
end

Buses=[Bus(:,1:2),x(1:nbus),x(nbus+1:2*nbus),Pi,Qi]


% Transferred Power Calculation
%===========================

for f=1:length(A(:,1))
    i=A(f,1);k1=A(f,2);
    for h=1:length(Bus(:,1))
        if Bus(h,1)==i;
           Vi=V(h);
           ci=c(h);
        end
        if Bus(h,1)==k1;
           Vk=V(h);
           ck=c(h);
        end
    end
    Yik=1/((A(f,4))+j*(A(f,5)));    
    Gik=real(Yik);
    Bik=imag(Yik);
    circ=A(f,3);
    Bii=A(f,6)/2;
    aik=A(f,8);
    Ptrs1=+Vi^2*1/aik^2*Gik-Vi*Vk*1/aik*(Gik*cos(ci-ck)+Bik*sin(ci-ck)); 
    Qtrs1=-Vi^2*1/aik^2*(Bii+Bik)-Vi*Vk*1/aik*(Gik*sin(ci-ck)-Bik*cos(ci-ck)); 
    Ptrs2=+Vk^2*Gik-Vk*Vi*1/aik*(Gik*cos(ck-ci)+Bik*sin(ck-ci));  
    Qtrs2=-Vk^2*(Bii+Bik)-Vk*Vi*1/aik*(Gik*sin(ck-ci)-Bik*cos(ck-ci)); 
    flows1(f,:)=[i,k1,circ,Ptrs1,Qtrs1];
    flows2(f,:)=[k1,i,circ,Ptrs2,Qtrs2];
end
    Flows=[flows1;flows2]