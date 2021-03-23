clc;
clear all;
close all;

L=32;
k=log2(L);
Nbits=10000; %πλήθος bits (μήκος ακολουθίας)
Nsymb=Nbits/k; %πλήθος συμβόλων
x=randi([0,1],[1,Nbits]); % τυχαία δυαδική ακολουθία 10000 bits

%κωδικοποίηση Gray
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

nsamp=16; %συντελεστής υπερδειγμάτισης
delay=4; %group delay
filtorder=delay*nsamp*2; %τάξη φίλτρου
rolloff=0.25; %συντελεστής πτώσης
%κρουστική απόκριση φίλτρου
rNyquist=rcosine(1,nsamp,'fir/sqrt',rolloff,delay);

%υπερδειγμάτιση και εφαρμογή φίλτρου rNyquist
y1=upsample(y,nsamp);
ytx=conv(y1,rNyquist); clear y1;
%φιλτράρισμα σήματος
yrx=conv(ytx,rNyquist);
%περικοπη λόγω καθυστέρησης 
yrx=yrx(2*delay*nsamp+1:end-2*delay*nsamp);
%σχεδίαση yrx και υπέρθεση της εισόδου
figure(1);
grid;
plot(yrx(1:10*nsamp)); hold;
stem([1:nsamp:nsamp*10],y(1:10),'filled');
% το φάσμα του σήματος στο δέκτη
figure(2);
pwelch(yrx,[],[],[],nsamp); 