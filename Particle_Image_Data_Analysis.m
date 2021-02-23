% PIV Data Analysis
clc;
clear all;
close all;

n = 399; % number of text files in folder i.e. number of frames in folder
for i = 1:n; 
    % reading files 
    if i < 10;
        datfile = ['PIVlab_000' num2str(i) '.txt']
        fid = fopen(datfile, 'r');                    
        A = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 3);
        x(i,:) = A{1};
        y(i,:) = A{2};
        u(i,:) = A{3};
        v(i,:) = A{4};
        vor(i,:) = A{5};
        fclose(fid);
    end
    % reading files
    if i >= 10 && i < 100;
        datfile = ['PIVlab_00' num2str(i) '.txt']
        fid = fopen(datfile, 'r');                    
        A = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 3);
        x(i,:) = A{1};
        y(i,:) = A{2};
        u(i,:) = A{3};
        v(i,:) = A{4};
        vor(i,:)= A{5};
        fclose(fid);
    end
    % reading files
    if i >= 100 && i < 1000;
        datfile = ['PIVlab_0' num2str(i) '.txt']
        fid = fopen(datfile, 'r');                    
        A = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 3);
        x(i,:) = A{1};
        y(i,:) = A{2};
        u(i,:) = A{3};
        v(i,:) = A{4};
        vor(i,:) = A{5};
        fclose(fid);
    end
    % reading files
    if i >= 1000;
        datfile = ['PIVlabs_' num2str(i) '.txt']
        fid = fopen(datfile, 'r');                    
        A = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 3);
        x(i,:) = A{1};
        y(i,:) = A{2};
        u(i,:) = A{3};
        v(i,:) = A{4};
        vor(i,:) = A{5};
        fclose(fid);
    end

    t(i) = i*.01; 
    
end

% Making a folder for the results of collected data
mkdir Results

cd Results

save('x.mat','x');
save('y.mat','y');
save('u.mat','u');
save('v.mat','v');
save('vor.mat','vor');
save('t.mat','t');
%% U plot and contour with animation
u = load('u.mat');
u = struct2cell(u);
u = cell2mat(u);
u = u(:,1150);

% formatting
figure();
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.6 0.9]);
subplot(3,1,1)
plot(t,u);
title('Velocity (u) at y/R = 0, with respect to Time');
xlabel('Time (s)');
ylabel('u (m/s)');

% formatting rows/finding unique values and reshaping data
u = load('u.mat');
u = struct2cell(u);
u = cell2mat(u);
x = x(1,:);
y = y(1,:);
x = unique(x);
y = unique(y);
u = u(200,:);
u = reshape(u,[39, 59]);

% plotting and formatting
subplot(3,1,2)
contourf(x,y,u,'LineStyle','none')
grid on;
xlabel('x (m)')
ylabel('y (m)')
title('Contour Plot of u')
h = colorbar();
ylabel(h, 'u (m/s)')

% formatting rows/finding unique values and reshaping data
t = load('t.mat');
t = struct2cell(t);
t = cell2mat(t);
u = load('u.mat');
u = struct2cell(u);
u = cell2mat(u);
u_original = u;

% creating array for specific rows/columns
x = x(1,:);
y = y(1,:);
x_original = unique(x);
y_original = unique(y);

% converting and restructing rows and specific cells
x = load('x.mat');
x = struct2cell(x);
x = cell2mat(x);
y = load('y.mat');
y = struct2cell(y);
y = cell2mat(y);

% plotting and creating the animation file for 60 seconds of the data
subplot(3,1,3)
for i = 1:length(t)
    Frame_Count = i
    u = reshape(u_original(i,:), [39, 59]);
    contourf(x_original,y_original,u,100,'LineStyle','none')
    caxis([-0.03 0.03]);
    grid on;
    xlabel('x (m)')
    ylabel('y (m)')
    title('Contour Animation of u')
    h = colorbar();
    ylabel(h, 'u (m/s)')
    % https://www.mathworks.com/help/matlab/ref/imwrite.html
    % https://www.mathworks.com/help/matlab/ref/getframe.html
    frame = getframe(gcf);
    active_frame =  frame2im(frame);
    [active_frame,colormap] = rgb2ind(active_frame,256);
    % https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
     if i == 1
    imwrite(active_frame,colormap,'u_contour_anim.gif','gif','LoopCount',1,'DelayTime',0.01);
    else
    imwrite(active_frame,colormap,'u_contour_anim.gif','gif','WriteMode','append','DelayTime',0.01);
    end
end

%% V plot and contour with animation

% formatting
v = load('v.mat');
v = struct2cell(v);
v = cell2mat(v);
v = v(:,1150);

% formatting rows/finding unique values and reshaping data
figure();
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.8 0.8]);
subplot(2,2,1)
plot(t,v)
title('Velocity (v) at y/R = 0, with respect to Time');
xlabel('Time (s)');
ylabel('v (m/s)');
subplot(2,2,3)
plot(t,v)
title('Velocity (v) at y/R = 0 zoomed out, with respect to Time');
xlabel('Time (s)');
ylabel('v (m/s)');
ylim([-1 1])

% formatting rows/finding unique values and reshaping data
v = load('v.mat');
v = struct2cell(v);
v = cell2mat(v);
x = x(1,:);
y = y(1,:);
x = unique(x);
y = unique(y);
v = v(200,:);
v = reshape(v,[39, 59]);

% plotting and formatting
subplot(2,2,2)
contourf(x,y,v,'LineStyle','none')
grid on;
xlabel('x (m)')
ylabel('y (m)')
title('Contour Plot of v')
h = colorbar();
ylabel(h, 'v (m/s)')

% formatting rows/finding unique values and reshaping data
t = load('t.mat');
t = struct2cell(t);
t = cell2mat(t);
v = load('v.mat');
v = struct2cell(v);
v = cell2mat(v);
v_original = v;
x = x(1,:);
y = y(1,:);
x_original = unique(x);
y_original = unique(y);

% converting and restructing rows and specific cells
x = load('x.mat');
x = struct2cell(x);
x = cell2mat(x);
y = load('y.mat');
y = struct2cell(y);
y = cell2mat(y);
subplot(2,2,4)

% plotting and creating the animation file for 60 seconds of the data
for i = 1:length(t)
    Frame_Count = i
    v = reshape(v_original(i,:), [39, 59]);
    contourf(x_original,y_original,v,20,'LineStyle','none')
    caxis([-1.5e-3 1.5e-3]);
    grid on;
    xlabel('x (m)')
    ylabel('y (m)')
    title('Contour Animation of v')
    h = colorbar();
    ylabel(h, 'v (m/s)')
    % https://www.mathworks.com/help/matlab/ref/imwrite.html
    % https://www.mathworks.com/help/matlab/ref/getframe.html
    frame = getframe(gcf);
    active_frame =  frame2im(frame);
    [active_frame,colormap] = rgb2ind(active_frame,256);
    % https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
     if i == 1
    imwrite(active_frame,colormap,'v_contour_anim.gif','gif','LoopCount',1,'DelayTime',0.1);
    else
    imwrite(active_frame,colormap,'v_contour_anim.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end
%% Vorticity plot and contour with animation
vor = load('vor.mat');
vor = struct2cell(vor);
vor = cell2mat(vor);
vor = vor(:,1150);
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.6 0.9]);
subplot(3,1,1)
plot(t,vor);
title('Velocity (vorticity) at y/R = 0, with respect to Time');
xlabel('Time (s)');
ylabel('Vorticity (1/s)');

% formatting rows/finding unique values and reshaping data
vor = load('vor.mat');
vor = struct2cell(vor);
vor = cell2mat(vor);
x = x(1,:);
y = y(1,:);
x = unique(x);
y = unique(y);
vor = vor(300,:);
vor = reshape(vor,[39, 59]);

% plotting
subplot(3,1,2)
contourf(x,y,vor,'LineStyle','none')
grid on;

% plotting and formatting
xlabel('x (m)')
ylabel('y (m)')
title('Contour Plot of Vorticity')
h = colorbar();
ylabel(h, 'Vorticity (1/s)')
t = load('t.mat');
t = struct2cell(t);
t = cell2mat(t);

% formatting rows/finding unique values and reshaping data
x = x(1,:);
y = y(1,:);
x_original = unique(x);
y_original = unique(y);
vor = load('vor.mat');
vor = struct2cell(vor);
vor = cell2mat(vor);
vor_original = vor;

% converting and restructing rows and specific cells
x = load('x.mat');
x = struct2cell(x);
x = cell2mat(x);
y = load('y.mat');
y = struct2cell(y);
y = cell2mat(y);
subplot(3,1,3)

% plotting and creating the animation file for 60 seconds of the data
for i = 1:length(t)
    Frame_Count = i
    vor = reshape(vor_original(i,:), [39, 59]);
    contourf(x_original,y_original,vor,100,'LineStyle','none')
    caxis([-50 50]);
    grid on;
    xlabel('x (m)')
    ylabel('y (m)')
    title('Contour Animation of vorticity')
    h = colorbar();
    ylabel(h, 'Vorticity (1/s)')
    % https://www.mathworks.com/help/matlab/ref/imwrite.html
    % https://www.mathworks.com/help/matlab/ref/getframe.html
    frame = getframe(gcf);
    active_frame =  frame2im(frame);
    [active_frame,colormap] = rgb2ind(active_frame,256);
    % https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
     if i == 1
    imwrite(active_frame,colormap,'vor_contour_anim.gif','gif','LoopCount',1,'DelayTime',0.1);
    else
    imwrite(active_frame,colormap,'vor_contour_anim.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end


% Animation data

% t = load('t.mat');
% t = struct2cell(t);
% t = cell2mat(t);
% u = load('u.mat');
% u = struct2cell(u);
% u = cell2mat(u);
% u_original = u;
% v = load('v.mat');
% v = struct2cell(v);
% v = cell2mat(v);
% v_original = v;
% vor = load('vor.mat');
% vor = struct2cell(vor);
% vor = cell2mat(vor);
% vor_original = vor;
% 
% 
% x = x(1,:);
% y = y(1,:);
% x_original = unique(x);
% y_original = unique(y);
% 
% x = load('x.mat');
% x = struct2cell(x);
% x = cell2mat(x);
% y = load('y.mat');
% y = struct2cell(y);
% y = cell2mat(y);
