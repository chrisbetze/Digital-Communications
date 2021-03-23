function errors=ask_errors_new_sin(k,M,nsamp,EbNo)
% � ��������� ���� ���������� ��� �������� ��� ���������������
% ���������� ������� L-ASK ��� ����� ��� ������ ��� ���������� ��������.
% ���������� ��� ������ ��� ���������� �������� (��� ��������� errors).
% k ����� � ������� ��� bits/�������, �������� L=2^k -- � ������� ���
% ������������ ������
% M ����� � ������� ��� ����������� �������� (����� ���������� L-ASK)
% nsamp � ������� ��� ��������� ��� ������� (oversampling ratio)
% EbNo ����� � ��������� ��������������� ����� Eb/No, �� db
%
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ��� ������ �������
% �������� ������� �������� {�1, �3, ... �(L-1)}. �� �����������
x=2*floor(L*rand(1,M))-L+1;
Px=(L^2-1)/3; % ��������� ����� �������

h=sin(pi*(1:nsamp)/nsamp);  h=h/sqrt(h*h');
y=upsample(x,nsamp); % ��������� ��� ����� ������
y=conv(y,h); % �� ���� ������� ����
y=y(1:M*nsamp); % ������������ � ���� ��� ������ � ��������

ynoisy=awgn(y,SNR,'measured'); % ��������� ����

for i=1:nsamp matched(i)=h(i); end
yrx=conv(ynoisy,matched);
z = yrx(nsamp:nsamp:M*nsamp); % Y������������ ��� ����� ���� �������� �
l=[-L+1:2:L-1];
for i=1:length(z)
  [m,j]=min(abs(l-z(i)));
  z(i)=l(j);
end
err=not(x==z);
errors=sum(err);
end