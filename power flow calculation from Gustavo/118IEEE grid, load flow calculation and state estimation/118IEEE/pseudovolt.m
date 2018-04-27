function pseudovolt(VPMU,IPMU,B,Bus)
global zpseudo

nbranch=length(B(:,1));
nbus=length(Bus(:,1));
nvolp=length(VPMU(:,1));
nI=length(IPMU(:,1));


for i=1:length(Bus(:,1))
     if Bus(i,2)==3 || Bus(i,2)==2
        slack=i;
     end
end

% Derivatives Calculation
%============================

% J1 for magnitude and J2 for angles
% 1-angle vol 2-mag vol 3-ang current 4-mag current

sigV=zeros(nI,2);
sigc=zeros(nI,2);
idefV=zeros(nI,2);
idefc=zeros(nI,2);

sV=zeros(nI,1);
sc=zeros(nI,1);
ivol=zeros(nI,1);
iang=zeros(nI,1);

i1=1;
 for f=1:length(IPMU(:,1))
      for r=1:nbus
          if Bus(r,1)==IPMU(f,1)
             i=r;
          end
          if Bus(r,1)==IPMU(f,2)
             k=r;
          end
      end 
      for s=1:nbranch
           if IPMU(f,1)==B(s,1) && IPMU(f,2)==B(s,2)
              if B(s,3)==IPMU(f,3)
                 aik=B(s,7);
                 gik=1/aik*B(s,4);
                 bik=1/aik*B(s,5);
                 Bii=B(s,6)/2+(1/aik-1)*1/aik*B(s,5);
                 Gii=(1/aik-1)*1/aik*B(s,4);
              end
           end
           if IPMU(f,1)==B(s,2) && IPMU(f,2)==B(s,1)
              if B(s,3)==IPMU(f,3)
                 aik=B(s,7);
                 gik=1/aik*B(s,4);
                 bik=1/aik*B(s,5);
                 Bii=B(s,6)/2+(1-1/aik)*B(s,5);
                 Gii=(1-1/aik)*B(s,4);
              end
           end
      end
       
    Iik=IPMU(f,4);
    cik=IPMU(f,5);
    for w=1:nvolp
        if VPMU(w,1)==IPMU(f,1)
           id=w;
        end
    end
    Vi=VPMU(id,3);
    ci=VPMU(id,4);
    a=gik*(Gii+gik)+bik*(Bii+bik);
    b=gik*(Bii+bik)-bik*(Gii+gik);
    cp=gik^2+bik^2;
    
    VjR=(a*Vi*cos(ci)-b*Vi*sin(ci)-gik*Iik*cos(cik)-bik*Iik*sin(cik))/cp;
    VjI=(b*Vi*cos(ci)+a*Vi*sin(ci)+bik*Iik*cos(cik)-gik*Iik*sin(cik))/cp; 
    
    J11=(VjI*(-b/cp*Vi*sin(ci)+a/cp*Vi*cos(ci))+VjR*(-a/cp*Vi*sin(ci)-b/cp*Vi*cos(ci)))/sqrt(VjR^2+VjI^2);
    J12=(VjI*(b/cp*cos(ci)+a/cp*sin(ci))+VjR*(a/cp*cos(ci)-b/cp*sin(ci)))/sqrt(VjR^2+VjI^2);
    J13=(VjI*(-bik*Iik/cp*sin(cik)-gik*Iik/cp*cos(cik))+VjR*(gik*Iik/cp*sin(cik)-bik*Iik/cp*cos(cik)))/sqrt(VjR^2+VjI^2);
    J14=(VjI*(bik/cp*cos(cik)-gik/cp*sin(cik))+VjR*(-gik/cp*cos(cik)-bik/cp*sin(cik)))/sqrt(VjR^2+VjI^2);
    J21=(VjR*(-b/cp*Vi*sin(ci)+a/cp*Vi*cos(ci))-VjI*(-a/cp*Vi*sin(ci)-b/cp*Vi*cos(ci)))/(VjR^2+VjI^2);
    J22=(VjR*(b/cp*cos(ci)+a/cp*sin(ci))-VjI*(a/cp*cos(ci)-b/cp*sin(ci)))/(VjR^2+VjI^2);
    J23=(VjR*(-bik*Iik/cp*sin(cik)-gik*Iik/cp*cos(cik))-VjI*(gik*Iik/cp*sin(cik)-bik*Iik/cp*cos(cik)))/(VjR^2+VjI^2);  
    J24=(VjR*(bik/cp*cos(cik)-gik/cp*sin(cik))-VjI*(-gik/cp*cos(cik)-bik/cp*sin(cik)))/(VjR^2+VjI^2);
   
    P=[VPMU(id,5) 0 0 0;0 VPMU(id,6) 0 0;0 0 IPMU(f,6) 0; 0 0 0 IPMU(f,7)]^2;
    J1=[J12,J11,J14,J13];
    J2=[J22,J21,J24,J23];
    J=[J1;J2];
    
    Pnew=J*P*J';
     
    sV(f)=sqrt(VjR^2+VjI^2); 
    sigV(f,:)=Pnew(1,:);
    idefV(f,:)=[k,f];
    ivol(f,1)=1; % define if Pseudo measurement correspond to magnitude
    sc(f,1)=atan2(VjI,VjR);
    sigc(f,:)=Pnew(2,:);
    idefc(f,:)=[k,f];
    iang(f,1)=0;
    
   
 end

 zpseudo1=[idefV  ivol sV sigV];
 zpseudo2=[idefc iang sc sigc];

 i1=1;
 for i=1:length(zpseudo1(:,1))
     zpseudo(i1,:)=zpseudo1(i,:);
     zpseudo(i1+1,:)=zpseudo2(i,:);
     i1=i1+2;
 end
 
 zpseudo;

