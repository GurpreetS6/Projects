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
% plane displacements. 
clear all 
close all
clc
I1 = imread('11x11.tiff');   % reading in image
[Ir, Ic]=size(I1) % matrix size
 IARstart=1; % defining starting parameters
 IARfin=Ir;
 IACstart=1;
 IACfin=Ic;
   Ic = I1(IARstart:IARfin,IACstart:IACfin);
   [RA, CA]=size (Ic);
%    figure(2) % figure
%
% Ic = imcrop(I1,[140 150 400 500]) ; % cropping image
 %%a= imread('right.tif');
% a=imcrop(a1,[82 64 11 9]);
% a=imread('cross2.jpg');
% a=255-a;
% a=255-a1; 
imshow(Ic); % display image
[x(1,1),y(1,1)]=ginput(1)
%x(1,1);
hold on
plot(x(1,1),y(1,1),'+b') % plot values on image

[x(2,1),y(2,1)]=ginput(1) % plot values on image depending on image size
hold on
plot(x(2,1),y(2,1),'+b')

drawnow % auto updates figure

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
figure(2), imshow(Ic); % show updated image

figure(3), imshow(a);
% loop through and find the best value to detect dots on screen
c = normxcorr2(a,Ic);
%figure, surf(c), shading flat
    sc = size(c); C = zeros(size(c));
    for ii=1:sc(1)
        for jj=1:sc(2)
            if c(ii,jj) > 0.645088488000355 & c(ii,jj)<.99 % change this
                C(ii,jj) = 1;
            end
        end
    end
L = bwlabel(C,8) ;%I,[low high])find(0 < x & x < 10*pi)
s  = regionprops(L, 'centroid'); % find centroid of dot
centroids = cat(1, s.Centroid); 
corr_offset = [ centroids(:,1)-size(a,2)/2 centroids(:,2)-size(a,1)/2]; 
%close all; 
figure(3), imshow(Ic); 
hold all; 
plot(corr_offset(:,1),corr_offset(:,2),'+r')
%plot(centroids(1:5,1),centroids(1:5,2),'ob')

% correlation offset variables
xc1=corr_offset(1:end,1);
yc1=corr_offset(1:end,2);
lx=length(xc1)
ly=length(yc1)
% dx1=diff(xc1);
% dy1=diff(yc1);
% dr1=sqrt(dx1.^2+dy1.^2);
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
% Ldata(:,1)=xc1;
% Ldata(:,2)=yc1;
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
nd=0;% number of points to delete from bottom i.e. last nd of , yc array.
figure (28) % new image
imshow(Ic);
hold all
plot(corr_offset(1:lx-nd,1),corr_offset(1:ly-nd,2),'og') % plot and determine location of next dots
plot(corr_offset(1:2,1),corr_offset(1:2,2),'oy')
xc=corr_offset(1:lx-nd,1)
yc=corr_offset(1:lx-nd,2)
 out_C1=[xc  yc];
save('OUTA_Left.txt','out_C1','-ascii')
% xc=corr_offset(:,1)
% yc=corr_offset(:,2)
out_xcyc=[xc(1:length(xc-1)),yc(1:length(xc-1))];
Ng=[11,11,11,11,11,11,11,11,11,11,11]; % array depending on the size of dots on screen 11x11
dX=20; % change based on image
dY=20
for j=1:1:11 % iterate through numbers
k1(j)=sum((Ng(1:1:j))) 
end
rxc=reshape(xc,11,11) % reshape matrix to 11x11
rxc=rxc;
%diff(rxc(1,:))
% [Dxx,Dxy]=gradient(rxc)
ryc=reshape(yc,11,11)

[rycn, I]=sort(ryc,'descend'); % descending dots on image
for j2=1:1:11
rxcn(:,j2)=rxc(I(:,j2),j2)
 end
[Dxx,Dxy]=gradient(rxcn)
[Dyx,Dyy]=gradient(rycn)
% for j=1:lx
%     Dx=diff(rxc(j,:))
% end
% X(1:k1(1))=-5;
% Y(1:k1(1))=-[2,1,0,-1,-2];% Positive Y is downwards!
% Z(1:k1(1))=6.4;
% %2column
% X(1+k1(1):k1(2))=-4;
% Y(1+k1(1):k1(2))=-[3,2,1,0,-1,-2,-3];
% Z(1+k1(1):k1(2))=0;
% % 3rd column
% X(1+k1(2):k1(3))=-3;
% Y(1+k1(2):k1(3))=-[4,3,2,1,0,-1,-2,-3, -4];
% Z(1+k1(2):k1(3))=6.4;
% % 4th Column
% X(1+k1(3):k1(4))=-2;
% Y(1+k1(3):k1(4))=-[5,4,3,2,1,0,-1,-2,-3, -4, -5];
% Z(1+k1(3):k1(4))=0;
% % 5th Column
% X(1+k1(4):k1(5))=-1;
% Y(1+k1(4):k1(5))=-[5,4,3,2,1,0,-1,-2,-3, -4, -5];
% Z(1+k1(4):k1(5))=6.4;
% % 6th Column
% X(1+k1(5):k1(6))=0;
% Y(1+k1(5):k1(6))=-[5,4,3,2,1,0,-1,-2,-3, -4, -5];
% Z(1+k1(5):k1(6))=0;
% % 7th Column
% X(1+k1(6):k1(7))=1;
% Y(1+k1(6):k1(7))=-[5,4,3,2,1,0,-1,-2,-3, -4, -5];
% Z(1+k1(6):k1(7))=6.4;
% % 8th Column
% X(1+k1(7):k1(8))=2;
% Y(1+k1(7):k1(8))=-[5,4,3,2,1,0,-1,-2,-3, -4, -5];
% Z(1+k1(7):k1(8))=0;
% % 9th Column
% X(1+k1(8):k1(9))=3;
% Y(1+k1(8):k1(9))=-[4,3,2,1,0,-1,-2,-3, -4];
% Z(1+k1(8):k1(9))=6.4;
% % 10th Column
% X(1+k1(9):k1(10))=4;
% Y(1+k1(9):k1(10))=-[3,2,1,0,-1,-2,-3];
% Z(1+k1(9):k1(10))=0;
% % 11th Column
% X(1+k1(10):k1(11))=5;
% Y(1+k1(10):k1(11))=-[2,1,0,-1,-2];
% Z(1+k1(10):k1(11))=6.4;



%outXYZ=[xc, yc]
%out_xcyc=[xc(1:length(xc-1)),yc(1:length(xc-1)),dX*X', dY*Y', Z'];

%save('StatNewPlate_out_xs_ys_LeftREF_Crosses.txt','out_xcyc','-ascii')%
%save('StaticAlu_out_xs_ys_WhiteLeft1_0.txt','out_xcyc','-ascii')
 
