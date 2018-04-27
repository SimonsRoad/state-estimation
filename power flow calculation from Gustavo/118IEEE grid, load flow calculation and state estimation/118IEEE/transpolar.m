function transpolar(ipmu,nim,nic)
global IPMUn

IPMUn=zeros(2*length(ipmu),6);

for i=1:length(ipmu(:,1))
    Ireal(i,1:3)=ipmu(i,1:3);
    Iimag(i,1:3)=ipmu(i,1:3);
    Ireal(i,4)=ipmu(i,4)*cos(ipmu(i,5));
    Iimag(i,4)=ipmu(i,4)*sin(ipmu(i,5));
    
    P=[nim(i) 0;0 nic(i)]^2;
    J=[cos(ipmu(i,5)) -ipmu(i,4)*sin(ipmu(i,5));sin(ipmu(i,5)) ipmu(i,4)*cos(ipmu(i,5))];
    Pnew=J*P*J';
    Ireal(i,5:6)=Pnew(1,:);
    Iimag(i,5:6)=Pnew(2,:);
end

 i1=1;
 for i=1:length(Ireal(:,1))
     IPMUn(i1,:)=Ireal(i,:);
     IPMUn(i1+1,:)=Iimag(i,:);
     i1=i1+2;
 end