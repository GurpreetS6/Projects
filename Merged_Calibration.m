%CALIBRATION FILE FOR SPLIT VIEW DIC
%Run first + to find the centerof the plate and then the dot.
%A template can be defined with one element i.e a + in an image called  "X"
%below. This image is convoluted through ffts and iffts with the imageA
%which contains several synbols of +s and dots . Through this convolution
%only the correlated symbols are retained in the imageA ex. only +s.
%This convolution acts like a Dirac's delta which is everywhere zero except
%at the small area of image a with +. Thus it returns the value of the
%fuction at the location of the +.
% Three 3-D target plates is used which requires manual introduction of outof
% plane dispalcements.
clear all 
close all
clc
I1 = imread('merged.tiff');   
[Ir, Ic]=size(I1)
 IARstart=1;
 IARfin=Ir;
 IACstart=1;
 IACfin=Ic/2;
   Ic = I1(IARstart:IARfin,IACstart:IACfin);
   [RA, CA]=size (Ic);
%    figure(2)
%
% Ic = imcrop(I1,[140 150 400 500]) ;
 %%a= imread('right.tif');
% a=imcrop(a1,[82 64 11 9]);
% a=imread('cross2.jpg');
% a=255-a;
% a=255-a1;
imshow(Ic);
[x(1,1),y(1,1)]=ginput(1)
%x(1,1);
hold on
plot(x(1,1),y(1,1),'+b')

[x(2,1),y(2,1)]=ginput(1)
hold on
plot(x(2,1),y(2,1),'+b')

drawnow

xmin = min(x)
xmax = max(x)
ymin = min(y)
ymax = max(y)

lowerline=[xmin ymin; xmax ymin];
upperline=[xmin ymax; xmax ymax];
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

plot(lowerline(:,1),lowerline(:,2),'-r')
plot(upperline(:,1),upperline(:,2),'-r')
plot(leftline(:,1),leftline(:,2),'-r')
plot(rightline(:,1),rightline(:,2),'-r')
%maxA=max(Image Ic(:))
a = Ic(ymin:ymax,xmin:xmax);


%% find the centroids of marks AFTER correlation
%imtool(a)
figure(2), imshow(Ic);

figure(3), imshow(a);

c = normxcorr2(a,Ic);
%figure, surf(c), shading flat
    sc = size(c); C = zeros(size(c));
    for ii=1:sc(1)
        for jj=1:sc(2)
            if c(ii,jj) > 0.55000355 & c(ii,jj)<.99
                C(ii,jj) = 1;
            end
        end
    end
L = bwlabel(C,8) ;%I,[low high])find(0 < x & x < 10*pi)
s  = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid); 
corr_offset = [ centroids(:,1)-size(a,2)/2 centroids(:,2)-size(a,1)/2]; 
%close all; 
figure(3), imshow(Ic); 
hold all; 
plot(corr_offset(:,1),corr_offset(:,2),'+r')
%plot(centroids(1:5,1),centroids(1:5,2),'ob')
xc=corr_offset(:,1);
yc=corr_offset(:,2);
dx=diff(xc);
dy=diff(yc);
dr=sqrt(dx.^2+dy.^2);
% figure(10),hist(dr,90);
% grid on
% We provide by hand the coordinates of the center found by  running first 
% cross-corelation with+
% points=[xc,yc];
% xc0=194; yc0=264.5;
% %% find the point next to origin
% %[idx_ul, xcn1, ycn1]=FindNearPoint(points,xc0, yc0, 10);
% %plot(xul,yul,'+g');
% x=xc-xc0;
% y=yc-yc0;
% dis=abs(x)<8 & abs((yc)-yc0)<8;
% idx=find(dis, 10,'first')
% % use centroid instead of correlation's coordinates to achieve higher 
% % accuracy.
% xcn=xc(idx)
% ycn=yc(idx)
% hold off
% figure(4), imshow(Ic);
% hold on
% plot(xc,yc,'og');
% %plot(xcn1,ycn1,'+r');
% plot(xcn,ycn,'+r');
%DELETING BAD POINTS
Ldata(:,1)=xc;
Ldata(:,2)=yc;
% [x_d y_d]=ginput;
% D(:,1)=x_d
% D(:,2)=y_d
% hold on
% 
% for i=1:length(x_d)
% xL=x_d(i)-Ldata(:,1);
% yL=y_d(i)-Ldata(:,2);
% disLL=abs(xL)<9 & abs(yL)<9;
% idLL=find(disLL, 1,'first');
% p(i)=idLL
% xcL(i)=Ldata(idLL,1);
% ycL(i)=Ldata(idLL,2);
% hold on
% figure (24)
% plot(xcL(i),ycL(i),'ob')
% plot(Ldata(:,1),Ldata(:,2),'+g')
% grid on 
% hold on
% drawnow
% hold on
% % Ldata(idLL,1)=[];
% % Ldata(idLL,2)=[];
% end
% 
% %Kdata(p,:)
% Dx=Ldata(:,1)
% Dy=Ldata(:,2)
% Dx(p) = []
% Dy(p) = []
% Kdata=[Dx,Dy]
% corr_offset=Kdata;
figure (28)
imshow(Ic);
hold all
plot(corr_offset(:,1),corr_offset(:,2),'og')
 out_Centr1=[corr_offset];
%save('OUTA_Left.txt','out_Centr1','-ascii')
xc=corr_offset(:,1)
yc=corr_offset(:,2)
out_xcyc=[xc(1:length(xc-1)),yc(1:length(xc-1))];
%save('StatAlu_out_xs_ys_LeftREF_Crosses.txt','out_xcyc','-ascii')%
save('StaticAlfu_out_xs_ys_WhiteLet1_0.txt','out_xcyc','-ascii')
W = load('OUTA_Left.txt');
D = load('zvalues.txt');
Z = D';
X = W(:,1);
Y = W(:,2);
xc = xc(1:length(X));
yc = yc(1:length(Y));
Done = [xc yc X Y Z];
save('Done.txt','Done','-ascii')
