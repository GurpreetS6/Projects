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
I1 = imread('11x11_3D.tiff');   % reading in 3D image file
[Ir, Ic]=size(I1) % sizing image
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
imshow(Ic);   % selecting and finding first dot on the image
[x(1,1),y(1,1)]=ginput(1)
%x(1,1);
hold on
plot(x(1,1),y(1,1),'+b')

[x(2,1),y(2,1)]=ginput(1) % selecting and finding second dot on the image
hold on
plot(x(2,1),y(2,1),'+b')

drawnow
xmin = min(x) % minimum x value
xmax = max(x) % maximum x value
ymin = min(y) % minimum y value
ymax = max(y) % maximum y value

% assigning lines based on min max values
lowerline=[xmin ymin; xmax ymin]; 
upperline=[xmin ymax; xmax ymax]; 
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

% plotting lines based on min max values
plot(lowerline(:,1),lowerline(:,2),'-r')
plot(upperline(:,1),upperline(:,2),'-r')
plot(leftline(:,1),leftline(:,2),'-r')
plot(rightline(:,1),rightline(:,2),'-r')
%maxA=max(Image Ic(:))
a = Ic(ymin:ymax,xmin:xmax);


%% find the centroids of marks AFTER correlation
%imtool(a)
figure(2), imshow(Ic);  % show updated image
 
figure(3), imshow(a);

c = normxcorr2(a,Ic);
%figure, surf(c), shading flat
    sc = size(c); C = zeros(size(c));
% loop through and find the best value to detect dots on screen
    for ii=1:sc(1)
        for jj=1:sc(2)
            if c(ii,jj) > 0.6488000355 & c(ii,jj)<.99
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
%corr_offset=Kdata;
figure (28) 
imshow(Ic);
hold all
plot(corr_offset(1:50,1),corr_offset(1:50,2),'og') % offset image
 out_Centr1=[corr_offset];
%save('OUTA_Left.txt','out_Centr1','-ascii')
xc=corr_offset(:,1)
yc=corr_offset(:,2)
out_xcyc=[xc(1:length(xc-1)),yc(1:length(xc-1))];
Ng=[11,11,11,11,11,11,11,11,11,11,11]; % 11x11 matrix
dX=20;
dY=20;
for j=1:1:11
k1(j)=sum((Ng(1:1:j)))
end

X(1:k1(1))=-5;
Y(1:k1(1))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];% Positive Y is downwards!
Z(1:k1(1)-7)=0;
Z(k1(1)-6:k1(1))=12;
%2column
% making dots in form of matrix and assigning them values
X(1+k1(1):k1(2))=-4;
Y(1+k1(1):k1(2))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(1):k1(2)-7)=0;
Z(k1(2)-6:k1(2))=12;
% 3rd column
X(1+k1(2):k1(3))=-3;
Y(1+k1(2):k1(3))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(2):k1(3)-7)=0;
Z(k1(3)-6:k1(3))=12;
% 4th Column
X(1+k1(3):k1(4))=-2;
Y(1+k1(3):k1(4))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(3):k1(4)-7)=0;
Z(k1(4)-6:k1(4))=12;
% 5th Column
X(1+k1(4):k1(5))=-1;
Y(1+k1(4):k1(5))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(4):k1(5)-7)=0;
Z(k1(5)-6:k1(5))=12;
% 6th Column
X(1+k1(5):k1(6))=0;
Y(1+k1(5):k1(6))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(5):k1(6))=0;
% 7th Column
X(1+k1(6):k1(7))=1;
Y(1+k1(6):k1(7))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(6):k1(7))=0;
% 8th Column
X(1+k1(7):k1(8))=2;
Y(1+k1(7):k1(8))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(7):k1(8))=0;
% 9th Column
X(1+k1(8):k1(9))=3;
Y(1+k1(8):k1(9))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(8):k1(9))=0;
% 10th Column
X(1+k1(9):k1(10))=4;
Y(1+k1(9):k1(10))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(9):k1(10))=0;
% 11th Column
X(1+k1(10):k1(11))=5;
Y(1+k1(10):k1(11))=-[5,4,3,2,1,0,-1,-2,-3,-4,-5];
Z(1+k1(10):k1(11))=0;
% load values in from 2D answers
W = load('OUTA_Left.txt');
X = W(:,1);
Y = W(:,2);
%outXYZ=[X', Y', Z']
out_xcyc=[xc(1:length(xc-1)),yc(1:length(xc-1)),dX*X, dY*Y, Z'];
% save data into new txt files
save('StatNewPlate_out_xs_ys_LeftREF_Crosses.txt','out_xcyc','-ascii')%
save('zvalues.txt','Z','-ascii')
%save('StaticAlu_out_xs_ys_WhiteLeft1_0.txt','out_xcyc','-ascii')
 
