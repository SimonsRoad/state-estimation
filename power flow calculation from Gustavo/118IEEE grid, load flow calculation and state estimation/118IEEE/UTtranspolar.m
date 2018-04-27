function UTtranspolar(ipmu,nim,nic)
global IPMUUT

IPMUUT=zeros(2*length(ipmu),6);

for i=1:length(ipmu(:,1))
    r=ipmu(i,4);
    c=ipmu(i,5);
    sdr=nim(i);
    sdc=nic(i);
    x=[r;c];
    P=[sdr 0;0 sdc]^2;
    % Parameters
    n=2;
    a=1;
    k=(3-n)*0;
    c1=(n+k);
    % Weigths
    Wom=k/(n+k);
    Woc=k/(n+k);
    Wi=1/(2*(n+k));
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
    X=repmat(x,1,2*n+1) + sqrt(c1)*[zeros(n,1) sqrt(P) -sqrt(P)];
    Y=zeros(size(X));
    for j=1:length(X(1,:))
        Y(1,j)=X(1,j)*cos(X(2,j));
        Y(2,j)=X(1,j)*sin(X(2,j));
    end

    uk=Y*wm;
    Sk=Y*W*Y';
    Pnew=(Sk);
    Ireal(i,1:3)=ipmu(i,1:3);
    Iimag(i,1:3)=ipmu(i,1:3);
    Ireal(i,4)=uk(1);
    Iimag(i,4)=uk(2);
    Ireal(i,5:6)=Pnew(1,:);
    Iimag(i,5:6)=Pnew(2,:);
end

 i1=1;
 for i=1:length(Ireal(:,1))
     IPMUUT(i1,:)=Ireal(i,:);
     IPMUUT(i1+1,:)=Iimag(i,:);
     i1=i1+2;
 end

