function jacobianlf(Bus,Ym)
global J Ang Vol
clear ang vol;
%
% vectors of angles
i1=1;
for i=1:length(Bus(:,1))
    if Bus(i,2)~=2 && Bus(i,2)~=3
    ang(i1,1)=Bus(i,4);
    angnum(i1,1)=i;
    i1=i1+1;
    end
end
%
% vectors of voltages
ix=1;
for i=1:length(Bus(:,1))
    if Bus(i,2)==0 || Bus(i,2)==4
    vol(ix,1)=Bus(i,3);
    volnum(ix,1)=i;
    ix=ix+1;
    end
end
%
Vol=vol;
Ang=ang;
%
Q=zeros(length(Bus(:,1)),1);P=zeros(length(Bus(:,1)),1);
J1c=[];J1v=[];J2c=[];J2v=[];
% 
i1c=1;i1v=1;i2c=1;i2v=1;% aux
for i=1:length(Bus(:,1)) 
    for j=1:length(Bus(:,1))
        Q(i)=Q(i)+Bus(i,3)*Bus(j,3)*(real(Ym(i,j))*sin(Bus(i,4)-Bus(j,4))-imag(Ym(i,j))*cos(Bus(i,4)-Bus(j,4)));
        P(i)=P(i)+Bus(i,3)*Bus(j,3)*(real(Ym(i,j))*cos(Bus(i,4)-Bus(j,4))+imag(Ym(i,j))*sin(Bus(i,4)-Bus(j,4)));
    end
    i1=1;i2=1; % auxiliar
    if Bus(i,2)~=3 &&  Bus(i,2)~=2
       for j=1:length(ang(:,1))
       u1=angnum(j,1);
       J1c(i1c,i1)=Bus(i,3)*Bus(u1,3)*(real(Ym(i,u1))*sin(Bus(i,4)-Bus(u1,4))-imag(Ym(i,u1))*cos(Bus(i,4)-Bus(u1,4))); % 4.22 Wood
         if u1==i
         J1c(i1c,i1)=-Q(i)-imag(Ym(i,i))*Bus(i,3)^2;  
         end
         i1=i1+1;
       end 
       i1c=i1c+1; 
    end
%     
    if Bus(i,2)==0 || Bus(i,2)==4
       for j=1:length(ang(:,1))
       u1=angnum(j,1);
       J2c(i2c,i2)=-Bus(i,3)*Bus(u1,3)*(real(Ym(i,u1))*cos(Bus(i,4)-Bus(u1,4))+imag(Ym(i,u1))*sin(Bus(i,4)-Bus(u1,4)));
         if u1==i
         J2c(i2c,i2)=P(i)-real(Ym(i,i))*Bus(i,3)^2; 
         end
         i2=i2+1;
       end
       i2c=i2c+1;
    end
    %
    i3=1;i4=1; % auxiliar
    if Bus(i,2)~=3 && Bus(i,2)~=2
        for j=1:length(vol(:,1))
        u1=volnum(j,1);
        J1v(i1v,i3)=Bus(i,3)*(real(Ym(i,u1))*cos(Bus(i,4)-Bus(u1,4))+imag(Ym(i,u1))*sin(Bus(i,4)-Bus(u1,4)));
          if u1==i
          J1v(i1v,i3)=P(i)/Bus(i,3)+real(Ym(i,i))*Bus(i,3);
          end
          i3=i3+1;
        end     
    i1v=i1v+1;
    end
%     
    if Bus(i,2)==0 || Bus(i,2)==4
        for j=1:length(vol(:,1))
       u1=volnum(j,1); 
       J2v(i2v,i4)=Bus(i,3)*(real(Ym(i,u1))*sin(Bus(i,4)-Bus(u1,4))-imag(Ym(i,u1))*cos(Bus(i,4)-Bus(u1,4)));
        if u1==i
        J2v(i2v,i4)=Q(i)/Bus(i,3)-imag(Ym(i,i))*Bus(i,3);
        end
        i4=i4+1;
        end
    i2v=i2v+1;
    end
 end % initial for
J=[J1c J1v;J2c J2v];
