clc;
clear all;
close all;

k=4;
M=50000;
nsamp=16;
EbNo=10;
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ��� ������ �������
% �������� ������� �������� {�1, �3, ... �(L-1)}. �� �����������
x=2*floor(L*rand(1,M))-L+1;
Px=(L^2-1)/3; % ��������� ����� �������
h=ones(1,nsamp); h=h/sqrt(h*h'); % ��������� �������� �������
 % ������ (����������� ������ ���������� ���������)
y=upsample(x,nsamp); % ��������� ��� ����� ������
y=conv(y,h); % �� ���� ������� ����
y=y(1:M*nsamp); % ������������ � ���� ��� ������ � ��������
ynoisy=y;%awgn(y,SNR,'measured'); % ��������� ����
for i=1:nsamp matched(i)=h(end-i+1); end
yrx=conv(ynoisy,matched);
z = yrx(nsamp:nsamp:M*nsamp); % Y������������ ��� ����� ���� �������� �
l=[-L+1:2:L-1];
for i=1:length(z)
  [m,j]=min(abs(l-z(i)));
  z(i)=l(j);
end
err=not(x==z);
errors=sum(err);