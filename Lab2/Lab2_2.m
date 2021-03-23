% Bandpass FIR filter - closed-form and Kaiser windowing
clear all; close all;
load sima; %σημα s και συχνότητα Fs=8192Hz
figure; pwelch(s,[],[],[],Fs); %φάσμα πριν το φιλτράρισμα
f1=750; f2=950; %συχνότητες bandpass φίλτρου
Ts=1/Fs;
%ιδανική ζωνοπερατή συνάρτηση
H=[zeros(1,f1) ones(1,f2-f1) zeros(1,2*(Fs/2-f2)) ones(1,f2-f1) zeros(1,f1)];
hbp=ifft(H, 'symmetric'); %αντίστροφος μετασχηματισμός Fourier
%ολίσθηση της κρουστικής απόκρισης
middle=length(hbp)/2;
hbp=[hbp(middle+1:end) hbp(1:middle)];
%κρουστική απόκριση ορθογωνικού παραθύρου 256+1 δειγμάτων
h256=hbp(middle+1-128:middle+129);  
hbp_kaiser=h256.*kaiser(length(h256),5)'; %κρουστική απόκριση παραθύρου Kaiser
wvtool(h256,hbp_kaiser); %οπτικοποίηση παραθύρων στο πεδίο χρόνου και συχνότητας
%φιλτράρισμα με τα δύο φίλτρα
y_rect=conv(s,h256);
figure; pwelch(y_rect,[],[],[],Fs);
y_kai=conv(s,hbp_kaiser);
figure; pwelch(y_kai,[],[],[],Fs);

% Bandpass FIR filter – Equiripple design (Parks Mc Clellan)
f=2*[0 f1*0.95 f1*1.05 f2*0.95 f2*1.05 Fs/2]/Fs;
hbp_pm=firpm(128, f, [0 0 1 1 0 0]); %κρουστική απόκριση equiripple design filter 
wvtool(hbp_pm); %οπτικοποίηση παραθύρων στο πεδίο χρόνου και συχνότητας
s_pm=conv(s,hbp_pm);
figure; pwelch(s_pm,[],[],[],Fs);