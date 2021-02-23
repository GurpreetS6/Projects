% Fourier Signals
clc, clear all, close all, tic
syms n analytical_a(n) analytical_b(n) x
% Formulas
T=2*pi; 
w0=2*pi/T; 
a0 = 0
analytical_a(n)= 0
analytical_b(n)= (2/n)*((-1)^(n+1))
final_b(n)= 2*(((sin(n*x))/n)*(-1)^(n+1)); % Final f(x) bn term
an = 0;
n = 1:1:100;
bn = (final_b(n)); %Final b_n value
% double(bn);
x = -5*pi:0.1:5*pi; % Range
TermSum = a0;
for n = 1:1:100 % Summation
 TermSum = (TermSum + 0 + bn(n)); % Adds terms
 if n==1 | n==5 | n==10 | n==25 | n==100 % Terms to plot
 fplot(TermSum); % Plot
 hold on;
 grid on;
 end
end
% formatting
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
ylabel('f(x)')
xlabel('x')
title({'Fourier series of f(x)=x Using Different Amounts of Terms'})
t1 = '1 term';
t2 = '5 terms';
t3 = '10 terms';
t4 = '25 terms';
t5 = '100 terms';
text = {t1, t2, t3, t4,t5};
g =legend (text);
ylim([-4 4])
xlim([-3*pi 3*pi])
figure() % new figure
% Fourier signals for new plot
for n = 100
 TermSum = (TermSum + 0 + bn(n)); % Adds terms
 fplot(TermSum); % Plot
 hold on;
 grid on;
 solved_TS = eval(TermSum); % Convert answer to int for Part 3
end
% formatting
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
ylabel('f(x)')
xlabel('x')
title({'Fourier series of f(x)=x Using 100 terms'})
text = {t5};
g =legend (text);
ylim([-4 4])
xlim([-3*pi 3*pi])
hold off;
fs = 100;
x_length = 315; % Workspace x length
fd = fs*(0:x_length-1)/x_length;
fftTS = fft(solved_TS); % fft
d = abs((2/x_length)*fftTS); % fft formula
figure()
plot(fd,d) % fft plot
% formatting
title('FFT of f(x)')
xlabel('Frequency')
ylabel('Amplitude')
grid on
% adding noise to signal
noise = -0.7 + (0.7+0.7).*rand(size(x)); % Amplitude = 0.7
signal = noise + solved_TS; % Noise added to signal
figure()
% formatting and labelling
subplot(3,1,1)
plot(x,solved_TS)
title('Original Signal')
xlabel('Frequency')
ylabel('Amplitude')
subplot(3,1,2)
plot(x,noise)
title('Noise to be Added to Signal')
xlabel('Frequency')
ylabel('Amplitude')
subplot(3,1,3)
plot(x,signal)
title('Signal with Noise')
xlabel('Frequency')
ylabel('Amplitude')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
cf = 6; % Testing cf values, 6 works well
f = cf/(fs/2);
figure()
[b a] = butter(5, f);
filter = filtfilt(b,a,signal); % filter
subplot(3,1,1)
plot(x,solved_TS)
title('Original Unfiltered Signal with no Noise')
xlabel('Frequency')
ylabel('Amplitude')
subplot(3,1,2)
plot(x,signal)
title('Unfiltered Signal with Noise')
xlabel('Frequency')
ylabel('Amplitude')
subplot(3,1,3)
plot(x,filter)
title('Filtered Signal')
xlabel('Frequency')
ylabel('Amplitude')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 0.5 0.5]);
grid on
% formatting for 8 plots
p=1;
 z=0;
 figure()
for s = 1:2:7
 noise = s+(s+s).*rand(size(x));
 signal = noise + solved_TS;
 fs = 100;
 cf = 6; % Testing cf values, 6 works well
 f = cf/(fs/2);
 [b a] = butter(5, f);
 filter = filtfilt(b,a,signal);
 subplot(4,2,p)
 p=p+2;
 plot(x,signal)
 title('Signal with Noise')
 xlabel('Frequency')
 ylabel('Amplitude')
 z=z+2;
 subplot(4,2,z)
 plot(x,filter)
 title('Filtered Signal')
 xlabel('Frequency')
 ylabel('Amplitude')
 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.75 0.75]);
 grid on
end
% adding annotations for plots
annotation('textbox', [.05, 0.80, 0, 0], 'string', 's=1')
annotation('textbox', [.05, 0.60, 0, 0], 'string', 's=3')
annotation('textbox', [.05, 0.40, 0, 0], 'string', 's=5')