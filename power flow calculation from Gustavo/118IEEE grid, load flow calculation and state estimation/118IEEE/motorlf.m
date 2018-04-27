function motorlf(Bus,Motor)
global Pmot Qmot x2 MotorLF


 M=length(Motor(:,1));
 MotorLF=Motor;
 %
 for i=1:M
    inmot (Bus,Motor(i,:))
    MotorLF(i,5)=Pmot;
    MotorLF(i,6)=Qmot;
    MotorLF(i,7)=x2;
 end
 



