% Part 1 ������������ �� ����

close all           % ������� ���� ��� �������� �����������
clear all           % ��������� ��� ���� ��������
clc                 % ��������� �� �������� ������� 

Fs=500;             % ��������� �������������� 500 Hz
Ts=1/Fs;            % �������� ��������������
L=1000;             % ����� ������� (������� ���������)
T=L*Ts;             % �������� ������� 
t=0:Ts:(L-1)*Ts;    % �������� ������� ����������� �� �������
x=sin(2*pi*20*t)+ 0.8*sin(2*pi*70*(t-2)); % ��������� ����

% ��������� �� ���� ��� ����� ��� ������
figure(1)                       % ������� ��������� ��� ������� ���������
plot(t,x)                       % ������� ��������� ��� �������
title('Time domain plot of x')  
% �������������� ����������
xlabel('t (sec)')               % ������� ���� ����� x
ylabel('Amplitude')             % ������� ���� ����� y
pause                           % �������, ������ ��� ������� ��� �� ����������
axis([0 0.2 -2 2])              % �������� ��� ������� ��� 0 ��� 0.2 sec ���
                                % ������� ��� -2 ��� 2
pause                           % �������, ������ ��� ������� ��� �� ����������

% ���������� ��� �������� �������������� Fourier
N = 2^nextpow2(L);      % ����� ��������������� Fourier. 
                        % � nextpow2 ������� ��� ������ ��� ������� ��� 2 ���
                        % ����� ���������� � ��� ��� �� ������ L
                        % �����������, =ceil(log2(L))
Fo=Fs/N;                % ������� ����������
f=(0:N-1)*Fo;           % �������� ���������� 
X=fft(x,N);             % ����������� ����������� ��� ��������� ���������������
                        % Fourier (DFT) ��� � ������

% ��������� �� ���� ��� ����� ����������
% ���� �� ���� ����� ���������� ��������
% �� ���������� ���� ��� ������� ����������
figure(2)                           % ������� ��������� ��� ������� ���������
plot(f(1:N-1),abs(X(1:N-1)))        % ������� ��������� ��� ������� ����������   
title('Frequency domain plot of x') % ������ �������� ����������
xlabel('f  (Hz)')                   % ������� ���� ����� x
ylabel('Amplitude')                 % ������� ���� ����� y
pause                               % ������� ��� �� ����� �� �����
% ������ ��� ������� ��� �� ����������

% ��� �� ������� ��������� ��� ����������� ��������
% ������ �� ��������������� ��� fftshift ���� � ���� ���
% �� ��������� ����� �� ����������� ���� ���� ��� ������
% ����� help fftshift ��� ������������ ������������
figure(3)           % ������� ��������� ��� ������� ���������
f=f-Fs/2;           % �������� ���������� ���� �� �������� ���� �Fs/2
X=fftshift(X);      % �������� ��� ��������� ���������� ��� ������ ��� ��������                
plot(f,abs(X));title('Two sided spectrum of x'); 
xlabel('f  (Hz)'); 
ylabel('Amplitude') 
pause               
% �������, ������ ��� ������� ��� �� ����������
% ���������� ��� ����
power=X.*conj(X)/N/L;       % ����������� ���������� ������
figure(4)                   % ������� ��������� ��� ������� ���������
plot(f,power)               % ����� ��� ��������� ���������
xlabel('Frequency (Hz)')    % ������� ���� ����� x
ylabel('Power')             % ������� ���� ����� y
title('{\bf Periodogram}')  % ������ ������������ �� ����� ��������
pause

disp('Part2') 

n=randn(1,L);                   % ���������� �������
figure(5)                       % ������� ��������� ��� ������� ���������
plot(t,n)                       % ������� ��������� ��� �������
title('Noise')                  % ������������������������
xlabel('t (sec)')               % ������� ���� ����� x
ylabel('Amplitude')             % ������� ���� ����� y
pause                           % �������, ������ ��� ������� ��� �� ����������
axis([0 0.2 -2 2])              % �������� ��� ������� ��� 0 ��� 0.2 sec ���
                                % ������� ��� -2 ��� 2
pause                           % �������, ������ ��� ������� ��� �� ����������

Nx=fft(n,N);
power_n=Nx .* conj(Nx)/N/L;       % ����������� ���������� ������
figure(6)                     % ������� ��������� ��� ������� ���������
plot(f,power_n)               % ����� ��� ��������� ���������
xlabel('Frequency (Hz)')    % ������� ���� ����� x
ylabel('Power')             % ������� ���� ����� y
title('{\bf Periodogram}')  % ������ ������������ �� ����� ��������
pause

s=n+x;
figure(7)                       % ������� ��������� ��� ������� ���������
plot(t,s)                       % ������� ��������� ��� �������
title('Signal with noise')      % ������������������������
xlabel('t (sec)')               % ������� ���� ����� x
ylabel('Amplitude')             % ������� ���� ����� y
axis([0 0.2 -2 2])              % �������� ��� ������� ��� 0 ��� 0.2 sec ���
                                % ������� ��� -2 ��� 2
pause                           % �������, ������ ��� ������� ��� �� ����������

S=fft(s,N);
S=fftshift(S);
figure(8)
plot(f,abs(S));
title('Two sided spectrum of s'); 
xlabel('f  (Hz)'); 
ylabel('Amplitude');
pause               % �������, ������ ��� ������� ��� �� ����������

disp('Part3')

c=sin(2*pi*100*t) .* s;
figure(9)                       % ������� ��������� ��� ������� ���������
plot(t,c)                       % ������� ��������� ��� �������
title('������������ ����')      % ������������������������
xlabel('t (sec)')               % ������� ���� ����� x
ylabel('Amplitude')             % ������� ���� ����� y
axis([0 0.2 -2 2])              % �������� ��� ������� ��� 0 ��� 0.2 sec ���
                                % ������� ��� -2 ��� 2
pause                           % �������, ������ ��� ������� ��� �� ����������

C=fft(c,N);
C=fftshift(C);
figure(8)
plot(f,abs(C));
title('Two sided spectrum of c'); 
xlabel('f  (Hz)'); 
ylabel('Amplitude');
pause               % �������, ������ ��� ������� ��� �� ����������
