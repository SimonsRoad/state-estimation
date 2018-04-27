% This is needed for load flow and post contingency
function Ymatrix(Bus,Branch,Shunt)
global Ym 
%
for i=1:length(Branch(:,1))
     Ybranch(i)=1/(Branch(i,4)+j*Branch(i,5));
end
%
% SIZE OF Ym
Ym=zeros(length(Bus(:,1)));
%
%Y MATRIX CONSTRUCTION
%
for i=1:length(Bus(:,1))% Y(i,i) calculation
    for l=1:length(Branch(:,1))
        r=Branch(l,8);
        if Branch(l,1)==Bus(i,1)
           Ym(i,i)=Ym(i,i)+(1/r-1)/r*Ybranch(l)+1/r*Ybranch(l)+j*Branch(l,6)/2;
       end  
       if Branch(l,2)==Bus(i,1)
          Ym(i,i)=Ym(i,i)+(1-1/r)*Ybranch(l)+1/r*Ybranch(l)+j*Branch(l,6)/2;
       end
    end
    if length(Shunt)~=0
        for s=1:length(Shunt(:,1))
            if Shunt(s,1)==Bus(i,1)
            Ym(i,i)=Ym(i,i)+Shunt(s,3)+j*Shunt(s,4); %Compensation inclusion in Ymatrix
            end
        end
    end
end
%
for i=1:length(Bus(:,1))%Y(i,j) calculation
     for u=1:length(Bus(:,1))
         if Bus(u,1)~=Bus(i,1)
            for m=1:length(Branch(:,1))
                r=Branch(m,8);
            if Branch(m,1)==Bus(i,1) && Branch(m,2)==Bus(u,1)
               Ym(i,u)=Ym(i,u)-Ybranch(m)/r; 
            end 
            if Branch(m,1)==Bus(u,1) && Branch(m,2)==Bus(i,1)
               Ym(i,u)=Ym(i,u)-Ybranch(m)/r; 
            end
            end 
         end
     end
end
