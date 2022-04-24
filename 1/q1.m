close all;
clc;

%% 1.1 - Data Generation:

[v_m,f_s] = audioread('in-the-air.wav'); %load the file
T_s = 1/f_s; %sampling interval
N = length(v_m); %number of sampling points
t = 0:T_s:(N-1)*T_s; %time vector
f = linspace(-f_s/2 , f_s/2 , N); %frequency vector
V_m = fftshift(fft(v_m))/sqrt(N); %v_m fourrier transform

figure;
subplot(2,1,1);
plot(t,v_m);
title('original signal'); ylabel('v_m'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,V_m);
ylabel('V_m'); xlabel('f [Hz]'); grid on;
%sound(v_m,f_s); %listen to the signal


%% 1.2 - Modulator:

f_c = 15000; %carrier frequency
K_AM = 0.02; %modulation index
v_AM = ammod(v_m,f_c,f_s,0,K_AM); %creation of v_AM
V_AM = fft(v_AM); % v_AM forrier tranform

figure;
subplot(2,1,1);
plot(t,v_AM);
title('AM with Carrier'); ylabel('v_A_M'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,V_AM);
ylabel('V_A_M'); xlabel('f [Hz]'); grid on;

%sound(v_AM,f_s); %listem to the signal


%% 1.3 - Channel:

N0 = 0.02; %noise size
z = N0*randn(1,N); %noise generation
z = z';
x_r = z+v_AM; %recieved signal generation
X_r = fftshift(fft(x_r))/sqrt(N); %x_r forrier tranform

figure;
subplot(2,1,1);
plot(t,x_r);
title('recived signal'); ylabel('x_r'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_r);
ylabel('X_r'); xlabel('f [Hz]'); grid on;


%% 1.4 - Demodulator:

x_L = bandpass(x_r,[8500 21000],f_s); %x_r through bandpass for x_L generation
X_L = fftshift(fft(x_L)) / sqrt(N); %x_L forrier tranform

figure;
subplot(2,1,1);
plot(t,x_L);
title('Filtered signal'); ylabel('x_L'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_L);
ylabel('X_L'); xlabel('f [Hz]'); grid on;

x_d = amdemod(x_L,f_c,f_s,0,K_AM); %demodulate x_L
x_d = lowpass(x_d,8500,f_s); %lowpass x_d
X_d = fftshift(fft(x_d)) / sqrt(N); %x_d fourrier transform

figure;
subplot(2,1,1);
plot(t,x_d);
hold on;
plot(t,v_m,'r');
title('Demodulated signel'); xlabel('t [sec]'); grid on;  legend([{'x_d(t)'};{'v_m(t)'}]);
subplot(2,1,2);
plot(f,X_d);
hold on;
plot(f,V_m,'r');
xlabel('f [Hz]'); grid on; legend([{'X_d(f)'};{'V_m(f)'}]);

%sound(x_d,f_s); %listem to the signal

diff = xcorr(x_d, v_m, 0, 'coeff'); %correlation between original signal and decoded signal


%% 1.5 - different noise and modulations:

%% 1.5.1 - different noise:

N1 = 0.1; %noise size
z1 = N1*randn(1,N); %noise generation
z1 = z1';
x_r1 = z1+v_AM; %recieved signal generation
X_r1 = fftshift(fft(x_r1))/sqrt(N); %x_r forrier tranform

figure;
subplot(2,1,1);
plot(t,x_r1);
title('Different noise - recived signal'); ylabel('x_r'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_r1);
ylabel('X_r'); xlabel('f [Hz]'); grid on;

x_L1 = bandpass(x_r1,[8500 21000],f_s); %x_r through bandpass for x_L1 generation
X_L1 = fftshift(fft(x_L1))/sqrt(N); %x_L1 forrier tranform

figure;
subplot(2,1,1);
plot(t,x_L1);
title('Different noise - Filtered signal'); ylabel('x_L'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_L1);
ylabel('X_L'); xlabel('f [Hz]'); grid on;

x_d1 = amdemod(x_L1,f_c,f_s,0,K_AM); %demodulate x_L1
x_d1 = lowpass(x_d1,8500,f_s); %lowpass x_d1
X_d1 = fftshift(fft(x_d1)) / sqrt(N); %x_d1 forrier tranform

figure;
subplot(2,1,1);
plot(t,x_d1);
hold on;
plot(t,v_m,'r');
title('Different noise - Demodulated signel'); xlabel('t [sec]'); grid on;  legend([{'x_d(t)'};{'v_m(t)'}]);
subplot(2,1,2);
plot(f,X_d1);
hold on;
plot(f,V_m,'r');
xlabel('f [Hz]'); grid on; legend([{'X_d(F)'};{'V_m(F)'}]);

%sound(x_d1,f_s); %listem to the signal

diff1 = xcorr(x_d1, v_m, 0, 'coeff'); %correlation between original signal and decoded signal


%% 1.5.2 - FM modulation:

f_d = 10000; %frequency deviation
v_FM = fmmod(v_m,f_c,f_s,f_d); %FM modulation
V_FM = fftshift(fft(v_FM))/sqrt(N); %v_FM forrier tranform

x_r_FM = z+v_FM; %recieved signal generation
X_r_FM = fftshift(fft(x_r_FM))/sqrt(N); %x_r_FM forrier tranform

figure;
subplot(2,1,1);
plot(t,x_r_FM);
title('FM Modulation - recived signal'); ylabel('x_r(t)'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_r_FM);
ylabel('X_r(f)'); xlabel('f [Hz]'); grid on;

x_L_FM = bandpass(x_r_FM,[8500 21000],f_s); %x_r_FM through bandpass for x_L_FM generation
X_L_FM = fftshift(fft(x_L_FM))/sqrt(N); %x_L_FM forrier tranform

figure;
subplot(2,1,1);
plot(t,x_L_FM);
title('FM Modulation - Filtered signal'); ylabel('x_L'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_L_FM);
ylabel('X_L'); xlabel('f [Hz]'); grid on;

x_d_FM = fmdemod(x_L_FM,f_c,f_s,f_d); %demodulate x_L_FM
x_d_FM = lowpass(x_d_FM,8500,f_s); %bandpass x_d_FM
X_d_FM = fftshift(fft(x_d_FM))/sqrt(N); %x_d_FM forrier tranform

figure;
subplot(2,1,1);
plot(t,x_d_FM);
hold on;
plot(t,v_m,'r');
title('FM Modulation - Demodulated signel'); xlabel('t [sec]'); grid on; legend([{'x_d(t)'};{'v_m(t)'}]);
subplot(2,1,2);
plot(f,X_d_FM);
hold on;
plot(f,V_m,'r');
xlabel('f [Hz]'); grid on; legend([{'X_d(F)'};{'V_m(F)'}]);

%sound(x_dFM,f_s); %listem to the signal

diff_FM = xcorr(x_d_FM, v_m, 0, 'coeff');


%% different noise and fm modulation

x_r_FM1 = z1+v_FM; %recieved signal generation
X_r_FM1 = fftshift(fft(x_r_FM1))/sqrt(N); %x_r_FM1 forrier tranform


figure;
subplot(2,1,1);
plot(t,x_r_FM1);
title('Different noise and FM Modulation - recived signal'); ylabel('x_r(t)'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_r_FM1);
ylabel('X_r(f)'); xlabel('f [Hz]'); grid on;

x_L_FM1 = bandpass(x_r_FM1,[8500 21000],f_s); %bandpass x_r_FM1
X_L_FM1 = fftshift(fft(x_L_FM1))/sqrt(N); %x_L_FM1 forrier tranform

figure;
subplot(2,1,1);
plot(t,x_L_FM1);
title('Different noise and FM Modulation - Filtered signal'); ylabel('x_L'); xlabel('t [sec]'); grid on;
subplot(2,1,2);
plot(f,X_L_FM1);
ylabel('X_L'); xlabel('f [Hz]'); grid on;

x_d_FM1 = fmdemod(x_L_FM1,f_c,f_s,f_d); %demodulate x_L_FM1
x_d_FM1 = lowpass(x_d_FM1,8500,f_s); %bandpass x_d_FM1
X_d_FM1 = fftshift(fft(x_d_FM1))/sqrt(N); %x_d_FM1 forrier tranform

figure;
subplot(2,1,1);
plot(t,x_d_FM1);
hold on;
plot(t,v_m,'r');
title('Different noise and FM Modulation - Demodulated signel'); xlabel('t [sec]'); grid on;  legend([{'x_d(t)'};{'v_m(t)'}]);
subplot(2,1,2);
plot(f,X_d_FM1);
hold on;
plot(f,V_m,'r');
xlabel('f [Hz]'); grid on; legend([{'X_d(F)'};{'V_m(F)'}]);

%sound(x_d1FM,f_s); %listem to the signal

diff_FM1 = xcorr(x_d_FM1, v_m, 0, 'coeff'); %correlation between original signal and decoded signal
