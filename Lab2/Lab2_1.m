clear all; close all;
load sima; %�� ������ ���� �������� �� ���� s ��� �� Fs
figure; pwelch(s,[],[],[],Fs); %����� ���� �� �����������

%������� ���������� ��������� �, �� ��������� �������� Fs/8
H=[ones(1,Fs/8) zeros(1,Fs-Fs/4) ones(1,Fs/8)];

%��������� �������� �� ���������� �������������� Fourier
h=ifft(H, 'symmetric');
middle=length(h)/2; 
%�������� ��� ���������� ��������� ���� �� ���� ��� �����
h=[h(middle+1:end) h(1:middle)];
%� ��������� �������� ������������ �� ���� 32+1,64+1 ��� 128+1 ���������
h32=h(middle+1-16:middle+17);
h64=h(middle+1-32:middle+33);
h128=h(middle+1-64:middle+65);

wvtool(h32,h64,h128); %������������ ��������� ��� ����� ������ ��� ����������

%��������������� �������� Hamming ��� Kaiser ������ 64+1
wh=hamming(length(h64));
wk=kaiser(length(h64),5);
figure; plot(0:64,wk,'r',0:64,wh,'b'); grid;
h_hamming=h64.*wh';
h_kaiser=h64.*wk';
wvtool(h64, h_hamming, h_kaiser);

%����������� �� ���� ��� �� ������ ��� �� ����� ������
y_rect=conv(s,h64);
figure; pwelch(y_rect,[],[],[],Fs);
y_hamm=conv(s,h_hamming);
figure; pwelch(y_hamm,[],[],[],Fs);
y_kais=conv(s,h_kaiser);
figure; pwelch(y_kais,[],[],[],Fs);

%���������� Parks-MacClellan
hpm=firpm(64, [0 0.10 0.15 0.5]*2, [1 1 0 0]);
s_pm=conv(s,hpm);
figure; pwelch(s_pm,[],[],[],Fs);
sound(20*s); %������ �� ������ ����
sound(20*s_pm); %������ �� ������������� ����