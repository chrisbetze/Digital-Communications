clear all;
close all;
clc;

k=mod(16067,2)+3;
L=2^k;
M=50000;
nsamp=16;
numBits=k*M;
EbNo=1:20;
Pe=((L-1)/L)*erfc(sqrt(3*k/(L^2-1)*(10.^(EbNo/10))));
theBER=Pe/k;
for i=1:20
    errors(i)=ask_errors(k,M,nsamp,i);
    simBER(i)=errors(i)/numBits;
end
figure(1);
hold on;
set(gca,'yscale','log');
semilogy(EbNo,simBER,'r+');
semilogy(EbNo,theBER,'b-');
title('BER of 16-ASK');
xlabel('Eb/No (dB)');
ylabel('BER');
hold off;