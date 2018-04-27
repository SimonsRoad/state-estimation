clear all;
initiate;
Xtrue=load('V Phasor PMU.txt');
Strue=[Xtrue(:,4);Xtrue(:,3)];

trials=100;
vs=zeros(trials,4);
Jmin=zeros(trials,4);
IND1=zeros(trials,4);
IND2=zeros(trials,4);

for mc=1:trials
    measurements; 
    clear States 
    estimator; % Constrained
    for mcj=1:length(x)
         vs(mc,1)=vs(mc,1)+(x(mcj)-Strue(mcj))^2;
    end
    Error=x-Strue;
    Err1A(:,mc)=abs(Error(1:nbus));
    Err1V(:,mc)=abs(Error(nbus+1:2*nbus));
    Jmin(mc,1)=J;
    IND1(mc,1)=Index1;
    IND2(mc,1)=Index2;
    ITER(mc,1)=iter;
    
    clear States
    HybridSE; % Rectangular
    for mcj=1:length(x)
         vs(mc,2)=vs(mc,2)+(x(mcj)-Strue(mcj))^2;
    end
    Error=x-Strue;
    Err2A(:,mc)=abs(Error(1:nbus));
    Err2V(:,mc)=abs(Error(nbus+1:2*nbus));
    Jmin(mc,2)=J;
    IND1(mc,2)=Index1;
    IND2(mc,2)=Index2;
    ITER(mc,2)=iter;
    
    clear States
    estimator3; % Pseudo Only linear
    for mcj=1:length(x)
         vs(mc,3)=vs(mc,3)+(x(mcj)-Strue(mcj))^2;
    end
    Error=x-Strue;
    Err3A(:,mc)=abs(Error(1:nbus));
    Err3V(:,mc)=abs(Error(nbus+1:2*nbus));
    Jmin(mc,3)=J;
    IND1(mc,3)=Index1;
    IND2(mc,3)=Index2;
    ITER(mc,3)=iter;
    
    clear States
    estimator0; % Classical State Estimator
    for mcj=1:length(x)
         vs(mc,4)=vs(mc,4)+(x(mcj)-Strue(mcj))^2;
    end
    Error=x-Strue;
    Err4A(:,mc)=abs(Error(1:nbus));
    Err4V(:,mc)=abs(Error(nbus+1:2*nbus));
    Jmin(mc,4)=J;
    IND1(mc,4)=Index1;
    IND2(mc,4)=Index2;
    ITER(mc,4)=iter;
end

VS=zeros(1,4);
JMIN=zeros(1,4);
ind1=zeros(1,4);
ind2=zeros(1,4);

for k=1:4
    for i=1:trials
        VS(1,k)=VS(1,k)+vs(i,k);
        JMIN(1,k)=JMIN(1,k)+Jmin(i,k);
        ind1(1,k)=ind1(1,k)+IND1(i,k);
        ind2(1,k)=ind2(1,k)+IND2(i,k);
    end
end

Suv1=zeros(118,1);
Suv2=zeros(118,1);
Suv3=zeros(118,1);
Suv4=zeros(118,1);
Sua1=zeros(118,1);
Sua2=zeros(118,1);
Sua3=zeros(118,1);
Sua4=zeros(118,1);

for i=1:trials
    Sua1=Sua1+Err1A(:,i);
    Sua2=Sua2+Err2A(:,i);
    Sua3=Sua3+Err3A(:,i);
    Sua4=Sua4+Err4A(:,i);
    Suv3=Suv3+Err3V(:,i);
    Suv2=Suv2+Err2V(:,i);
    Suv1=Suv1+Err1V(:,i);
    Suv4=Suv4+Err4V(:,i);
end

ErAngle=[Sua2,Sua3,Sua1];
ErMag=[Suv2,Suv3,Suv1];

VarSt=VS/trials
J=JMIN/trials
Idx1=ind1/trials
Idx2=ind2/trials

mean(ITER)

