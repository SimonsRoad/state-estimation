function UTtransvolt(VPMU,IPMU,B,Bus)
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
    
    % Parameters
    n=4;
    k1=(3-n);
    c1=(n+k1);
    
    % Weigths
    Wom=k1/(n+k1);
    Woc=k1/(n+k1);
    Wi=1/(2*(n+k1));
    Wii=zeros(2*n,1);
    for j=1:2*n
        Wii(j,1)=Wi;
    end
    wm=[Wom;Wii];
    wc=[Woc;Wii];
    Wm=zeros(length(wm(:,1)),2*n+1);
    for j=1:2*n+1 
        Wm(:,j)=wm;
    end
    W=(eye(size(Wm))-Wm)*diag(wc)*(eye(size(Wm))-Wm)';
    
    % Create sigma points
     x=[Vi;ci;Iik;cik];
     P=[VPMU(id,5) 0 0 0;0 VPMU(id,6) 0 0;0 0 IPMU(f,6) 0; 0 0 0 IPMU(f,7)]^2;
     X=repmat(x,1,2*n+1) + sqrt(c1)*[zeros(n,1) sqrt(P) -sqrt(P)];
     
     Y=zeros(2,length(X(1,:)));
     for j=1:length(X(1,:))
          VjR=(a*X(1,j)*cos(X(2,j))-b*X(1,j)*sin(X(2,j))-gik*X(3,j)*cos(X(4,j))-bik*X(3,j)*sin(X(4,j)))/cp;
          VjI=(b*X(1,j)*cos(X(2,j))+a*X(1,j)*sin(X(2,j))+bik*X(3,j)*cos(X(4,j))-gik*X(3,j)*sin(X(4,j)))/cp; 
          Y(1,j)=sqrt(VjR^2+VjI^2);
          Y(2,j)=atan2(VjI,VjR);
     end
     
     uk=Y*wm;
     Sk=Y*W*Y';
     Pnew=Sk;  % covariance matrix
    
     sV(f)=uk(1); 
     sigV(f,:)=Pnew(1,:);
     idefV(f,:)=[k,f];
     ivol(f,1)=1; % define if Pseudo measurement correspond to magnitude
     
     sc(f,1)=uk(2);
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
 
