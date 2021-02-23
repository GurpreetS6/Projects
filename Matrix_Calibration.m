 %% Camera Calibration
% November 2011
% This program determines all calibartion coefficients Lij of A & B 
% Input data  xs, ys and non-coplanar Xw, Yw, Zw
%It yields all camera paameter
clear all
close all
clc

[points_A_Grid]=load('Done.txt');
%[points_AWhiteDots]=load('StaticAluWhiteout_xsysXYZw_Leftt1a_1DEF.txt');
%naw=length(points_AWhiteDots);
lp=0;% Points to skip from end i.e. Bolts

[points_A]=[points_A_Grid];%points_AWhiteDots(1:naw-lp,:)];
%xc0A=194.5; yc0A=264.5;%Point of Origin of Calibartion target  on image A (left)
%Image Right
% [points_BGrid]=load('StaticALGridout_xsysXYZw_Right1a_1DEF.txt');
% %[points_BWhiteDots]=load('StaticAluWhiteout_xsysXYZw_Right1a_1DEF.txt');
% [LHomograpgy2D]=load('NewStaticCompCoeff_2D_0negYDEF.txt');% IMport from Homography.m
% LH2D=LHomograpgy2D(:,1:3)
% %nbw=length(points_BWhiteDots)
% [points_B]=[points_BGrid];%points_BWhiteDots(1:nbw-lp,:)];

%save('Alu_HomoAB_3D_120203_0_1negY.txt','L','-ascii')

% AAA=points_A;
% save('Left_View.txt','AAA','-ascii')
% BBB=points_B;
% save('Right_View.txt','BBB','-ascii')
[mA,nA]=size(points_A);
NA=length(points_A); %Number of Data Points%
XwAref=0.0;
YwAref=0.0;
XwA=points_A(:,3);
YwA=points_A(:,4);
ZwA=points_A(:,5);%-points_A(10,5);%Relative Zw ie.put Zw=0mm at the valleys and the 5mm at crests!
xsAref=0;% Transforming back to original image due to Imcrop
ysAref=0;
xsA=points_A(:,1)+xsAref;
ysA=points_A(:,2)+ysAref;
figure(122)
plot(xsA, -ysA,'*b')
       NT=11;% Number of coefficients to be determined
       MA=zeros(2*NA,NT);
       bA=zeros(2*NA,1);
       
       j=0;
       JL=0;
       JP=1;
       Lmax=10;
       AL=zeros(NT,Lmax);
       LA(9)=0;
       LA(10)=0;
       LA(11)=0;
       kA=0;
       CXA=0;%283;
       CYA=0;%420;
       
       for L=1:1:Lmax
           kA=kA+1
           j=0;
       for i=1:2*JP:2*NA-JL
           j=j+1;  
           R=LA(9)*XwA(j)+LA(10)*YwA(j)+LA(11)*ZwA(j)+1;
           MA(i,1)=XwA(j)/R;
           MA(i,2)=YwA(j)/R;
           MA(i,3)=ZwA(j)/R;
           MA(i,4)=1./R;
           MA(i+1,5)=XwA(j)/R;
           MA(i+1,6)=YwA(j)/R;
           MA(i+1,7)=ZwA(j)/R;
            MA(i+1,8)=1./R;
           MA(i,9)=-xsA(j)*XwA(j)/R;
           MA(i,10)=-xsA(j)*YwA(j)/R;
           MA(i,11)=-xsA(j)*ZwA(j)/R; 
                      
           MA(i+1,9)=-ysA(j)*XwA(j)/R;
           MA(i+1,10)=-ysA(j)*YwA(j)/R;
           MA(i+1,11)=-ysA(j)*ZwA(j)/R;
           
           ks(j)=xsA(j)-CXA;% ?? Which values
           et(j)=ysA(j)-CYA;
           r2(j)=(ks(j)*ks(j)+et(j)*et(j));
%              MA(i,12)=ks(j)*r2(j);
%              MA(i,13)=ks(j)*r2(j)*r2(j);
%             MA(i,14)=ks(j)*r2(j)*r2(j)*r2(j);     
%             MA(i,15)=(r2(j)*r2(j)+2*ks(j)*ks(j));
%             MA(i,16)=ks(j)*et(j);
% %             
%             MA(i+1,12)=et(j)*r2(j);
%            MA(i+1,13)=et(j)*r2(j)*r2(j);
%             MA(i+1,14)=et(j)*r2(j)*r2(j)*r2(j);
%             MA(i+1,16)=(r2(j)+2*et(j)*et(j));
%             MA(i+1,15)=ks(j)*et(j);
            
            
           bA(i)=xsA(j)/R;
           bA(i+1)=ysA(j)/R;
       end
       
       
HA=MA'*MA; %{H}={M}^T*{M}%
KA=MA'*bA; %{K}={M}^T*{b}%

%DD=(diag(H))
%Gauss-Seidel Method to Find {L} in the System of Equations {H}{L}={K} %
%LA=zeros(11,1);

%init=ones(length(HA),1); %Initial Gauss-Seidel Guess%
%LA=pinv(HA)*KA;%
LA=HA\KA;
 % LA=Gauss_Seidel(H,K,tol,itermax,init);
%Corresponding the variables with the pinhole thoery:
% [U,S,V]=svd(MA)
%MAverf=U*S*V';
% MAinv=V*pinv(S)*U'
% LA=MAinv*bA;
% [Ua,Sa,Va]=svd(HA);
% % %MAverf=U*S*V';
% HAinv=Va*pinv(Sa)*Ua';
% LA=HAinv*KA;
%LA=pinv(HAinv*HA)*LA;
%%LB=pinv(HBinv*HB)*LB;
LA11=LA(1);
LA12=LA(2);
LA13=LA(3);
LA14=LA(4);
LA21=LA(5);
LA22=LA(6);
LA23=LA(7);
LA24=LA(8);
LA31=LA(9);
LA32=LA(10);
LA33=LA(11);
%  LAd1=LA(12);
% LAd2=LA(13);
% LAd3=LA(14);
% LAd4=LA(15);
% LAd5=LA(16);
R3Ascale=(LA31^2+LA32^2+LA33^2)

CXA=[LA(1) LA(2) LA(3)]*[LA(9);LA(10);LA(11)]/R3Ascale;% Orthogonality between LA1i and LA3i
CYA=[LA(5) LA(6) LA(7)]*[LA(9);LA(10);LA(11)]/R3Ascale;% Orthogonality between LA2i and LA3i Also see Kwon's
CXA
CYA
CxA(kA)=CXA;
CyA(kA)=CYA;

R3sA(kA)=R3Ascale;

alfa=[LA31*XwA+LA32*YwA+ LA33*ZwA+1];
AL(:,kA)=LA;%(:,1);
       end
FXA2=(LA11-CXA*LA31)^2+(LA12-CXA*LA32)^2+(LA13-CXA*LA33)^2;
FXA=sqrt(FXA2);
FYA2=(LA21-CYA*LA31)^2+(LA22-CYA*LA32)^2+(LA23-CYA*LA33)^2;
FYA=sqrt(FYA2);
determinantA=det([LA(1) LA(2) LA(3);LA(5) LA(6) LA(7);LA(9) LA(10) LA(11)])
 RA11=(LA11-CXA*LA31)/FXA
RA12=(LA12-CXA*LA32)/FXA
RA13=(LA13-CXA*LA33)/FXA
RA21=(LA21-CYA*LA31)/FYA
RA22=(LA22-CYA*LA32)/FYA
RA23=(LA23-CYA*LA33)/FYA
TZA=1
TXA=(LA14-CXA*TZA)/FXA
TYA=(LA24-CYA*TZA)/FYA
RA31=LA31;
RA32=LA32;
RA33=LA33;
RA11^2+RA12^2+RA13^2
RA21^2+RA22^2+RA23^2
RA31^2+RA32^2+RA33^2
RAn1=[RA11;RA12;RA13]
RAn2=[RA21;RA22;RA23]
RAn3=[RA31;RA32;RA33]
RAc=cross(RAn1',RAn2')%Checking orthogonality or parallelism between Rc & Rn3
%RA3T=cross(R1',R2')
RAn=[RAn1';RAn2';RAn3']
[UArn,SARn,VARn]=svd(RAn)
RRAn=UArn*SARn*VARn'
I=[1,0,0;0,1,0;0,0,1];%Diagonal only
RAfn=UArn*I*VARn' % There are some differences between the 
RAfnn=UArn*VARn'% see page 48 of Sutton
IARn=RAfn'*RAfn
Out1A=[CXA,CYA,FXA,FYA,TXA,TYA,TZA]
Out2A=[RAfnn]
       GA=[LA11-CXA*LA31, LA12-CXA*LA32, LA13-CXA*LA33
           LA21-CYA*LA31, LA22-CYA*LA32, LA23-CYA*LA33
           LA31, LA32, LA33]
       gA=[-LA14;-LA24;-1     ]
       XoA=inv(GA)*gA  %Principal point According to Kwon
%        
       % Principal point
       
%    InputA=[XwA;YwA;ZwA];
%    alfa=alf*InputA;
% %IMAGE B
% [mB,nB]=size(points_B);
% NB=length(points_B); %Number of Data Points%
% XwBref=0;
% YwBref=0;
% xsBref=0;%1280/2;
% ysBref=0;
% XwB=points_B(:,3)+XwBref;
% YwB=points_B(:,4)+YwBref;
% ZwB=points_B(:,5);%-points_B(10,5);
% xsB=points_B(:,1)+xsBref;
% ysB=points_B(:,2)+ysBref;
figure (1)
plot3(XwA,YwA,ZwA,'*b')
grid on
hold on 
% plot3(XwB,YwB,ZwB,'Ob')
%        %ZwB=base_points_B(:,3);
%     MB=zeros(2*NB,NT);
%        bB=zeros(2*NB,1);
%        
%        j=0;
%        JL=0;
%        JP=1;
%        Lmax=10;
%        BL=zeros(NT,Lmax);
%        LB(9)=0;
%        LB(10)=0;
%        LB(11)=0;
%        kB=0;  
%        CXB=0;%283;
%        CYB=0;%420;
       
%        for L=1:1:Lmax
%            kB=kB+1
%            j=0;
%        for i=1:2*JP:2*NB-JL
%            j=j+1;     
%            Rb=LB(9)*XwB(j)+LB(10)*YwB(j)+LB(11)*ZwB(j)+1;
%            MB(i,1)=XwB(j)/Rb;
%            MB(i,2)=YwB(j)/Rb;
%            MB(i,3)=ZwB(j)/Rb;
%            MB(i,4)=1./Rb;
%            MB(i+1,5)=XwB(j)/Rb;
%            MB(i+1,6)=YwB(j)/Rb;
%            MB(i+1,7)=ZwB(j)/Rb;
%            MB(i+1,8)=1./Rb;%
%            MB(i,9)=-xsB(j)*XwB(j)/Rb;
%            MB(i+1,9)=-ysB(j)*XwB(j)/Rb;
%            MB(i,10)=-xsB(j)*YwB(j)/Rb;
%            MB(i+1,10)=-ysB(j)*YwB(j)/Rb;
%            MB(i,11)=-xsB(j)*ZwB(j)/Rb;
%            MB(i+1,11)=-ysB(j)*ZwB(j);
%            
%            ksb(j)=xsB(j)-CXB;
%            etb(j)=ysB(j)-CYB;
%             rb2(j)=(ksb(j)*ksb(j)+etb(j)*etb(j));
% %             MB(i,12)=ksb(j)*rb2(j);
% %             MB(i,13)=ksb(j)*rb2(j)*rb2(j);
% %             MB(i,14)=ksb(j)*rb2(j)*rb2(j)*rb2(j);
% %             MB(i,15)=(rb2(j)*rb2(j)+2*ksb(j)*ksb(j));
% %             MB(i,16)=ksb(j)*etb(j);
% % %             
% %              MB(i+1,12)=etb(j)*rb2(j);
% %              MB(i+1,13)=etb(j)*rb2(j)*rb2(j);
% %             MB(i+1,14)=etb(j)*rb2(j)*rb2(j)*rb2(j);
% %             MB(i+1,16)=(rb2(j)+2*etb(j)*etb(j));
% %             MB(i+1,15)=ksb(j)*etb(j);
%            
%     
%            bB(i)=xsB(j)/Rb;
%            bB(i+1)=ysB(j)/Rb;
%        end
%        
%       
%        HB=MB'*MB; %{H}={M}^T*{M}%
%        KB=MB'*bB; %{K}={M}^T*{b}%

%DD=(diag(H))
%Gauss-Seidel Method to Find {L} in the System of Equations {H}{L}={K} %
%LB=zeros(11,1);
% init=ones(length(HB),1); %Initial Gauss-Seidel Guess%
% LB=pinv(HB)*KB;%
 % LB=Gauss_Seidel(H,K,tol,itermax,init);
 %LB=HB\KB;
% [Ub,Sb,Vb]=svd(HB);
% %MAverf=U*S*V';
% HBinv=Vb*pinv(Sb)*Ub';
% LB=HBinv*KB;
%LB=pinv(HBinv*HB)*LB;
% LB=LB*10^10;
% LA=LA*10^10;
%Corresponding the variables with the pinhole thoery:
% LB11=LB(1);
% LB12=LB(2); 
% LB13=LB(3);
% LB14=LB(4);
% LB21=LB(5);
% LB22=LB(6);
% LB23=LB(7);
% LB24=LB(8);
% LB31=LB(9);
% LB32=LB(10);
% LB33=LB(11);
% %  LBd1=LB(12);
% % LBd2=LB(13);
% % LBd3=LB(14);
% % LBd4=LB(15);
% % LBd5=LB(16);
% 
% R3Bscale=(LB31^2+LB32^2+LB33^2);% For orthonormality rescaling. It is =1!! 
% 
% CXB=[LB(1) LB(2) LB(3)]*[LB(9);LB(10);LB(11)]/R3Bscale;% Orthogonality between LB1i and LB3i
% CYB=[LB(5) LB(6) LB(7)]*[LB(9);LB(10);LB(11)]/R3Bscale;% Orthogonality between LB2i and LB3i
% CxB(kB)=CXB;
% CyB(kB)=CYB;
% R3sB(kB)=R3Bscale;
% determinantB=det([LB(1) LB(2) LB(3);LB(5) LB(6) LB(7);LB(9) LB(10) LB(11)])
% 
% beta=[LB31*XwB+LB33*YwB+ LB33*ZwB+1];
% BL(:,kB)=LB;
% 
%        end
%     % LAa=[LA(1:4)';LA(5:8)'; LA(9:11)' 1]
% %LBb=[LB(1:4)';LB(5:8)'; LB(9:11)' 1]      
%        
            
    
LAa=[LA(1:4)';LA(5:8)'; LA(9:11)' 1]
%LBb=[LB(1:4)';LB(5:8)'; LB(9:11)' 1]
LAH=[LA11,LA12,LA14;LA21,LA22,LA24;LA31,LA32,1]
%Lh2D=LBH*inv(LAH)
%  Lh2D=Lh2D/Lh2D(3,3);
% LL1=LBb*LAa'
% LL2=inv(LAa*LAa');%3D Homography see Zisserman p. 
% Lhomog3D=LL1*LL2;
% LL1=LAa'
% LL2=inv(LAa'*LAa)
% % on(1:1:NA)=1;
%Lhomog3D=LBb*inv(LAa'*LAa)*LAa';
% Lhomog3D=LBb*LAa'*inv(LAa*LAa');
% Lhomog3D=Lhomog3D/Lhomog3D(3,3)
on(1:NA)=1;
% xa=[xsA,ysA,on']';
% for k=1:1:NA-1
%     X1=xa(:,k);
%    Wb(k)=Lh2D(3,:)*X1;
%    xsb(k)=Lh2D(1,:)*X1/Wb(k);
%     ysb(k)=Lh2D(2,:)*X1/Wb(k);
% end
% OUT2=[XwB(1:NB-1),YwB(1:NB-1),ZwB(1:NB-1),XWb',YWb',ZWb'];
% WN=[XWb' YWb' ZWb' WWb']';
% 
% for k=1:1:NB-1
%     W1=WN(:,k);
%     WWb(k)=LBb(3,:)*W1;
%       xsb(k)=LBb(1,:)*W1/WWb(k);
%     ysb(k)=LBb(2,:)*W1/WWb(k);
% end

%REcover xsa ysa from 3D calibartion matrix
oneA(1:NA)=1;
XA3D=[XwA YwA ZwA oneA']';
for k=1:1:NA
    W1=XA3D(:,k);
    WWa(k)=LAa(3,:)*W1;
      xsa3D(k)=LAa(1,:)*W1/WWa(k);
    ysa3D(k)=LAa(2,:)*W1/WWa(k);
end
out1=[LAa];%, LBb];
figure(22)
plot(xsA, -ysA,'*b')
hold on
grid on
plot(xsa3D',-ysa3D','Or')
meanxsa3D=mean(xsa3D'-xsA);
meanysa3D=mean(ysa3D'-ysA);
stdxsa3D=std(xsa3D'-xsA);
stdysa3D=std(ysa3D'-ysA);
% B
% oneB(1:NB)=1;
% XB3D=[XwB YwB ZwB oneB']';
% for k=1:1:NB
%     W1B=XB3D(:,k);
%     WWb(k)=LBb(3,:)*W1B;
%       xsb3D(k)=LBb(1,:)*W1B/WWb(k);
%     ysb3D(k)=LBb(2,:)*W1B/WWb(k);
% end

% figure(23)
% plot(xsB, -ysB,'*b')
% hold on
% grid on
% plot(xsb, -ysb,'ob')
% plot(xsb3D',-ysb3D','Or')
% meanxsb3D=mean(xsb3D'-xsB);
% meanysb3D=mean(ysb3D'-ysB);
% stdxsb3D=std(xsb3D'-xsB);
% stdysb3D=std(ysb3D'-ysB);
% plot(xsa3D(1:5)',-ysa3D(1:5)','Or')

save('AlCoeff3D_Mine_120515_0negYDEF.txt','out1','-ascii')
save('Astatic_Points_ADEF.txt','points_A','-ascii')
%save('Astatic_Points_BDEF.txt','points_B','-ascii')
DATAa=[xsA, ysA, xsa3D', ysa3D',XwA, YwA, ZwA];
%DATAb=[xsB, ysB, xsb3D', ysb3D',XwB, YwB, ZwB];
DATAERROR=[meanxsa3D, stdxsa3D,meanysa3D,stdysa3D];
           %meanxsb3D,stdxsb3D,meanysb3D, stdysb3D];
%plot(xbb,-ybb,'.b')
%plot(xbB,-ybB,'ob')
% Correct way for Homography 
%  GB=[LB11-CXB*LB31, LB12-CXB*LB32, LB13-CXB*LB33
%            LB21-CYB*LB31, LB22-CYB*LB32, LB23-CYB*LB33
%            LB31, LB32, LB33]
%        gB=[LB14;LB24;-1]
%        XoB=inv(GB)*gB  %Principal pointB??


 LA
%  LB

 N=49;%length(xsA); The rest of the data are Z not 0
 % 2-D Verification for points with Zw=0
MX=zeros(2,2);
   bA=zeros(2,1);
   XA=zeros(N,2); %[Xw,Yw] Matrix
%    initL=[1;1];
   for k=1:N
       MXA11=LA11-xsA(k)*LA31;
       MXA12=LA12-xsA(k)*LA32;
       MXA21=LA21-ysA(k)*LA31;
       MXA22=LA22-ysA(k)*LA32;
       MLAX=[MXA11,MXA12;MXA21,MXA22];
          
       bX1=xsA(k)-LA14;
       bX2=ysA(k)-LA24;
       bLX=[bX1;bX2];
             
       MMA=MLAX'*MLAX;
       bA=MLAX'*bLX;
       XA(k,:)=MMA\bA;
       DXWAR(k)=XA(k,1)-XwA(k);% difference in reconstruction; Kind of error
       DYWAR(k)=XA(k,2)-YwA(k);% difference in reconstruction
       alfa2D(k)=LA31*XA(k,1)+LA32*XA(k,2)+1;
   end
   meDXWAR=mean(DXWAR)
   stdDXWAR=std(DXWAR)
for j=1:1:N
     bMRA=[XA(j,1);XA(j,2);1];
     xyRA=LAH*bMRA;
     xsRA(j)=xyRA(1)/xyRA(3);
     ysRA(j)=xyRA(2)/xyRA(3);
     alfas(j)=xyRA(3);
     dxsRa(j)=xsRA(j)-xsA(j);
     dysRa(j)=ysRA(j)-ysA(j);
end
  figure (6)
  plot(xsRA,ysRA,'*b')
  grid on
  medxsRa=mean(dxsRa)
   stddxsRa=std(dxsRa)
    medysRa=mean(dysRa)
   stdysRa=std(dysRa)
   errorA=[medxsRa, stddxsRa,medysRa, stdysRa]; 
%    DX=DXWAR';
%    DXYr=(XA(:,1).*XA(:,1)+XA(:,2).*XA(:,2))/(XwA(:).*XwA(:)+YwA(:).*YwA(:));
%    dXY=sqrt(DXYr(:,1))-1;
 %hist(DYWAR*100);%,40)
 N=49
%    MXB=zeros(2,2);
%    bB=zeros(2,1);
%    XB=zeros(N,2); %[Xw,Yw] Matrix
%       ;%length(xsB);
%    for k=1:N
%        MXB11=LB11-xsB(k)*LB31;
%        MXB12=LB12-xsB(k)*LB32;
%        MXB21=LB21-ysB(k)*LB31;
%        MXB22=LB22-ysB(k)*LB32;
%        MLBX=[MXB11,MXB12;MXB21,MXB22];
%           
%        bXB1=xsB(k)-LB14;
%        bXB2=ysB(k)-LB24;
%        bLXB=[bXB1;bXB2];
%        
%        %XY(k,:)=Gauss_Seidel(MLX,bLX,tol,itermax,initL);
%        MMB=MLBX'*MLBX;
%        bB=MLBX'*bLXB;
%        XB(k,:)=MMB\bB;
%        DXWBR(k)=XB(k,1)-XwB(k);% difference in reconstruction; Kind of error
%        DYWBR(k)=XB(k,2)-YwB(k);% difference in reconstruction
%        beta2D=[LB31*XB(k,1)+LB33*XB(k,2)+1];
%    end
%    meDXWABR=mean(DXWBR)
%    stdDXWBR=std(DXWBR)
%      meDYWABR=mean(DYWBR)
%    stdDYWBR=std(DYWBR)
%    for j=1:1:N
%      bMRB=[XB(j,1);XB(j,2);1];
%      xyRB=LBH*bMRB;
%      xsRB(j)=xyRB(1)/xyRB(3);
%      ysRB(j)=xyRB(2)/xyRB(3);
%      betas(j)=xyRB(3); 
%      dxsRb(j)=xsRB(j)-xsB(j);
%      dysRb(j)=ysRB(j)-ysB(j); 
%   end
%     figure (7)
%   plot(xsRB,ysRB,'*r')
%   grid on
%   medxsRb=mean(dxsRb)
%    stddxsRb=std(dxsRb)
%     medysRb=mean(dysRb)
%    stdysRb=std(dysRb)
%   errorB=[medxsRb, stddxsRb,medysRb, stdysRb]; 
%   
%   Errorsensor1=[errorA;errorB]
     figure (11)
   plot(XA(:,1),XA(:,2),'or');
  
   grid on
   hold on
   plot(XwA,YwA,'.r')
%    figure (12)
%    plot(XwB,YwB,'.g')
%       hold on
%      %figure (12)
%    plot(XB(:,1),XB(:,2),'og')
%  
%    grid on
   figure (13)
   plot(XA(:,1),XA(:,2),'.r');
   hold on 
   grid on
   %plot(XB(:,1),XB(:,2),'.g')
    PA=[LA11 LA12 LA13 LA14
             LA21 LA22 LA23 LA24
             LA31 LA32 LA33  1];
              [KAA, RAAC, CAA, ppAA, pvAA] = decomposecamera(PA);%Zisserman's book p.155
              KAA=KAA/KAA(3,3);
              LeftCamA= [KAA, RAAC, CAA];
              
              jjjj=pvAA'*pvAA
%       PB=[LB11 LB12 LB13 LB14
%           LB21 LB22 LB23 LB24 
%               LB31 LB32 LB33  1];
%               [KBB, RBBC, CBB, ppBB, pvBB] = decomposecamera(PB);
%               KBB=KBB/KBB(3,3)
%     RightCamB= [KBB, RBBC, CBB];


% Form term Pb*CA which is PB*LB p.162 eq 6.13 Zisserman
lamda=1/50000
CAA1=[CAA;1]
%DL=lamda*PB*CAA1

%L=Lhomog3D/Lhomog3D(3,3);

ona(1:1:NA)=1;
xa1=[xsA, ysA, ona']';
% 
% for k=1:1:NA-1
%     X1=xa1(:,k);
%    Wh1(k)=(Lhomog3D(3,:)*X1);
%    xhb(k)=(Lhomog3D(1,:)*X1/Wh1(k)+DL(1));
%     yhb(k)=(Lhomog3D(2,:)*X1/Wh1(k)+DL(2));
%     Wh1(k)=(Lhomog3D(3,:)*X1+DL(3));
%     
% end

%Zisserman's eq 6.14
    p1 = PA(:,1)
    p2 = PA(:,2)
    p3 = PA(:,3)
    p4 = PA(:,4)  
%     PR1=PB(1,:)
%     PR2=PB(2,:)
%     PR3=PB(3,:)

    Ma = [p1 p2 p3]
    invMa = inv(Ma)
Addedterm=[-invMa*p4;1]
invMa = [inv(Ma);0 0 0]
meu=.100;
for kk=1:1:NA-1
    X1=xa1(:,kk);
   Wh2=invMa(3,:)*X1;%+Addedterm(3);
   xh2b=-invMa(1,:)*X1/Wh2;%+Addedterm(1);
    yh2b=invMa(2,:)*X1/Wh2;%+Addedterm(2);
    Wh2=invMa(3,:)*X1;%+Addedterm(3);
    %cn(kk)=PR3*[xh2b;yh2b;Wh2;0];
    %xn(kk)=PR1*[xh2b;yh2b;Wh2;0]/cn(kk)+Addedterm(1);
    %yn(kk)=PR2*[xh2b;yh2b;Wh2;0]/cn(kk)+Addedterm(2);
    %cn(kk)=PR3*[xh2b;yh2b;Wh2;0]+Addedterm(3);
%  cn(kk)=[xh2b;yh2b;Wh2;0];
%     xn(kk)=[xh2b;yh2b;Wh2;0]/cn(kk)+Addedterm(1);
%     yn(kk)=[xh2b;yh2b;Wh2;0]/cn(kk)+Addedterm(2);
%     cn(kk)=[xh2b;yh2b;Wh2;0]+Addedterm(3);

end


%First=[xsB(1:NB-1), ysB(1:NB-1), xhb',yhb'];
% figure (18)
%  plot(xsB, -ysB, '+b')
%  hold on
%  grid on
% plot(xhb', yhb', 'Or')
% plot((xh2b'), (yh2b'), 'Ob')
% grid on
% %  plot(xn'-588.712, -yn'-1146.2815, 'Og')
% plot(xn', yn', 'Og')
% 
%  grid on
%save('Alu_HomoAB_3D_120203_0_1negY.txt','L','-ascii')
%Mapping all valdxA points onto B view by using Homogarphy: Result is xsBb
%& ysBb  Good for Zw=0 points
 ii1=0;  
kk1=0;
% 
% for k1=1:NA
%     kk1=kk1+1;
%         bXA1=xsA(k1);
%        bXA2=ysA(k1);
%        bCA=[bXA1;bXA2;1];      
%        XAB(k1,:)=L*bCA;
%        xsBb(k1)=XAB(k1,1)/XAB(k1,3);
%        ysBb(k1)=XAB(k1,2)/XAB(k1,3);
%       
%  end
I2 = imread('merged.tiff');  
I1 = imrotate(I2, 90);
[Ir, Ic]=size(I1)
 IARstart=1;
 IARfin=Ir;
 IACstart=1;
 IACfin=Ic/3;
   A = I1(IARstart:IARfin,IACstart:IACfin);
   B= I1(IARstart:IARfin,Ic/2:Ic);
figure (2)
imshow(A)
hold on
% plot (xhb,yhb,'+g')
plot (xsRA,ysRA,'-r')% Reconstruction back on Image A is exact!!
hold on
plot (xsA,ysA,'+r')
%figure (3)
%imshow(B)
%hold on
% plot (xsB,ysB,'+r')
% plot(xsBb',ysBb','oy');% 3D homography It is wrong
% %plot(xbb,ybb,'Ob')
% plot(xhb,yhb,'*g')% 3D homography
% plot(xsb,ysb,'or')% 2D homography
% figure (4)
% imshow(B)
% hold on
% plot (xn,yn,'or')
% plot (xsRB,ysRB,'-r')% Reconstruction back on Image A is exact!!
% 
% %A=ImageA;
[RA,CA]=size(A)
% [RB, CB]=size(B)
RSHA=reshape(A,RA*CA,1);
MXA=zeros(2,2);
bXA=zeros(2,1);
k=0;
% LA31=0;LA32=0;LA33=0;
% LB31=0;LB32=0;LB33=0;
% for k=1:NB;
%     kk=kk+1;
%        MXA11=LA11-xsA(k)*LA31;
%        MXA12=LA12-xsA(k)*LA32;
%        MXA13=LA13-xsA(k)*LA33;
%        MXA21=LA21-ysA(k)*LA31;
%        MXA22=LA22-ysA(k)*LA32;
%        MXA23=LA23-ysA(k)*LA33;
%        
%        MXB11=LB11-xsB(k)*LB31;
%        MXB12=LB12-xsB(k)*LB32;
%        MXB13=LB13-xsB(k)*LB33;
%        MXB21=LB21-ysB(k)*LB31;
%        MXB22=LB22-ysB(k)*LB32;
%        MXB23=LB23-ysB(k)*LB33;
%        
%        MCX=[MXA11,MXA12,MXA13;MXA21,MXA22,MXA23;
%              MXB11,MXB12,MXB13;MXB21,MXB22,MXB23];
%           
%        bXA1=xsA(k)-LA14;
%        bXA2=ysA(k)-LA24;
%        
%        bXB1=xsB(k)-LB14;
%        bXB2=ysB(k)-LB24;
%        
%        bCX=[bXA1;bXA2;bXB1;bXB2];
%        
%        %XY(k,:)=Gauss_Seidel(MLX,bLX,tol,itermax,initL);
%        MC=MCX'*MCX;
%        bAB=MCX'*bCX;
%        XY=inv(MC)*MCX'*bCX;
%        inv(MC)
%        bAB
%        XY2=MC\bAB;
% 
%        DXWC(k)=XY(1);
%        DYWC(k)=XY(2);
%        DZWC(k)=XY(3);
%        
% end
% figure (25)
% plot3(XwA,YwA,ZwA,'*b')
% hold on
% grid on
% plot3(DXWC,DYWC,DZWC,'Or')
% Final=[XwA,YwA,ZwA,DXWC',DYWC',DZWC'];
% save('ALUMOUT_Reconstructed_Mine_120203_0negYDEF.txt','Final','-ascii')
% NormW=[XwA,YwA,ZwA]*[XwA,YwA,ZwA]';
% sum(NormW(:,1));