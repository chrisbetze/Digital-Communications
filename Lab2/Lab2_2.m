% Bandpass FIR filter - closed-form and Kaiser windowing
clear all; close all;
load sima; %���� s ��� ��������� Fs=8192Hz
figure; pwelch(s,[],[],[],Fs); %����� ���� �� �����������
f1=750; f2=950; %���������� bandpass �������
Ts=1/Fs;
%������� ���������� ���������
H=[zeros(1,f1) ones(1,f2-f1) zeros(1,2*(Fs/2-f2)) ones(1,f2-f1) zeros(1,f1)];
hbp=ifft(H, 'symmetric'); %����������� ��������������� Fourier
%�������� ��� ���������� ���������
middle=length(hbp)/2;
hbp=[hbp(middle+1:end) hbp(1:middle)];
%��������� �������� ����������� ��������� 256+1 ���������
h256=hbp(middle+1-128:middle+129);  
hbp_kaiser=h256.*kaiser(length(h256),5)'; %��������� �������� ��������� Kaiser
wvtool(h256,hbp_kaiser); %������������ ��������� ��� ����� ������ ��� ����������
%����������� �� �� ��� ������
y_rect=conv(s,h256);
figure; pwelch(y_rect,[],[],[],Fs);
y_kai=conv(s,hbp_kaiser);
figure; pwelch(y_kai,[],[],[],Fs);

% Bandpass FIR filter � Equiripple design (Parks Mc Clellan)
f=2*[0 f1*0.95 f1*1.05 f2*0.95 f2*1.05 Fs/2]/Fs;
hbp_pm=firpm(128, f, [0 0 1 1 0 0]); %��������� �������� equiripple design filter 
wvtool(hbp_pm); %������������ ��������� ��� ����� ������ ��� ����������
s_pm=conv(s,hbp_pm);
figure; pwelch(s_pm,[],[],[],Fs);