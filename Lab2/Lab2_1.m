clear all; close all;
load sima; %το αρχείο αυτό περιέχει το σήμα s και τη Fs
figure; pwelch(s,[],[],[],Fs); %φάσμα πριν το φιλτράρισμα

%ιδανική βαθυπερατή συνάρτηση Η, με συχνότητα αποκοπής Fs/8
H=[ones(1,Fs/8) zeros(1,Fs-Fs/4) ones(1,Fs/8)];

%κρουστική απόκριση με αντίστροφο μετασχηματισμό Fourier
h=ifft(H, 'symmetric');
middle=length(h)/2; 
%ολίσθηση της κρουστικής απόκρισης κατά το μισό της μήκος
h=[h(middle+1:end) h(1:middle)];
%η κρουστική απόκριση περικόπτεται σε μήκη 32+1,64+1 και 128+1 δειγμάτων
h32=h(middle+1-16:middle+17);
h64=h(middle+1-32:middle+33);
h128=h(middle+1-64:middle+65);

wvtool(h32,h64,h128); %οπτοκοποίηση παραθύρων στο πεδίο χρόνου και συχνότητας

%καρασκευάζονται παράθυρα Hamming και Kaiser μήκους 64+1
wh=hamming(length(h64));
wk=kaiser(length(h64),5);
figure; plot(0:64,wk,'r',0:64,wh,'b'); grid;
h_hamming=h64.*wh';
h_kaiser=h64.*wk';
wvtool(h64, h_hamming, h_kaiser);

%φιλτράρουμε το σήμα μας με καθένα από τα τρεία φίλτρα
y_rect=conv(s,h64);
figure; pwelch(y_rect,[],[],[],Fs);
y_hamm=conv(s,h_hamming);
figure; pwelch(y_hamm,[],[],[],Fs);
y_kais=conv(s,h_kaiser);
figure; pwelch(y_kais,[],[],[],Fs);

%Βαθυπερατό Parks-MacClellan
hpm=firpm(64, [0 0.10 0.15 0.5]*2, [1 1 0 0]);
s_pm=conv(s,hpm);
figure; pwelch(s_pm,[],[],[],Fs);
sound(20*s); %ακούμε το αρχικό σημα
sound(20*s_pm); %ακούμε το φιλτραρισμένο σήμα