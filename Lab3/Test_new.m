clc;
clear all;
close all;

k=4;
M=50000;
nsamp=16;
EbNo=10;
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
% Διάνυσμα τυχαίων ακεραίων {±1, ±3, ... ±(L-1)}. Να επαληθευθεί
x=2*floor(L*rand(1,M))-L+1;
Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
h=ones(1,nsamp); h=h/sqrt(h*h'); % κρουστική απόκριση φίλτρου
 % πομπού (ορθογωνικός παλμός μοναδιαίας ενέργειας)
y=upsample(x,nsamp); % μετατροπή στο πυκνό πλέγμα
y=conv(y,h); % το προς εκπομπή σήμα
y=y(1:M*nsamp); % περικόπτεται η ουρά που αφήνει η συνέλιξη
ynoisy=y;%awgn(y,SNR,'measured'); % θορυβώδες σήμα
for i=1:nsamp matched(i)=h(end-i+1); end
yrx=conv(ynoisy,matched);
z = yrx(nsamp:nsamp:M*nsamp); % Yποδειγμάτιση στο τέλος κάθε περιόδου Τ
l=[-L+1:2:L-1];
for i=1:length(z)
  [m,j]=min(abs(l-z(i)));
  z(i)=l(j);
end
err=not(x==z);
errors=sum(err);