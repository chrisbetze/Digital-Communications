clc;
clear all;
close all;

L=32;
k=log2(L);
Nbits=10000; %������ bits (����� ����������)
Nsymb=Nbits/k; %������ ��������
x=randi([0,1],[1,Nbits]); % ������ ������� ��������� 10000 bits

%������������ Gray
step=2;
mapping=[step/2; -step/2];
if(k>1)
 for j=2:k
 mapping=[mapping+2^(j-1)*step/2; -mapping-2^(j-1)*step/2];
 end
end
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb');
y=[];
for i=1:length(xsym)
 y=[y mapping(xsym(i)+1)];
end

nsamp=16; %����������� ���������������
delay=4; %group delay
filtorder=delay*nsamp*2; %���� �������
rolloff=0.25; %����������� ������
%��������� �������� �������
rNyquist=rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

%�������������� ��� �������� ������� rNyquist
y1=upsample(y,nsamp);
ytx=conv(y1,rNyquist); clear y1;
%����������� �������
yrx=conv(ytx,rNyquist);
%�������� ���� ������������ 
yrx=yrx(2*delay*nsamp+1:end-2*delay*nsamp);
%�������� yrx ��� �������� ��� �������
figure(1);
grid;
plot(yrx(1:10*nsamp)); hold;
stem([1:nsamp:nsamp*10],y(1:10),'filled');
% �� ����� ��� ������� ��� �����
figure(2);
pwelch(yrx,[],[],[],nsamp); 