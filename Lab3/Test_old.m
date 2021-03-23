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
sum(x.^2)/length(x); % ���������� ����� ������� (��� ����������)
y=rectpulse(x,nsamp);
n=wgn(1,length(y),10*log10(Px)-SNR);
ynoisy=y+n; % ��������� ����
y=reshape(ynoisy,nsamp,length(ynoisy)/nsamp);
matched=ones(1,nsamp);
z=matched*y/nsamp;
l=[-L+1:2:L-1];
for i=1:length(z)
  [m,j]=min(abs(l-z(i)));
  z(i)=l(j);
end
err=not(x==z);
errors=sum(err);