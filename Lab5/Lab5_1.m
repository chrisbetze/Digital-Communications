clc;
clear all;
close all;
% �������� mapping ��� ��� ������������ Gray M-QAM
% ����� �� ������ ���������� ������ �������, ��������� M=L^2
% l=log2(L): ������� bit ��� ��������� (inphase, quadrature)
M=64;
k=log2(M);
L=sqrt(M);
l=log2(L);
core=[1+1i;1-1i;-1+1i;-1-1i]; % ���������� ������������, M=4
mapping=core;
if(l>1)
 for j=1:l-1
 mapping=mapping+j*2*core(1);
 mapping=[mapping;conj(mapping)];
 mapping=[mapping;-conj(mapping)];
 end
end;
scatter(real(mapping), imag(mapping),'filled','d');
xlabel("In-Phase");
ylabel("Quadrature");
for i=1:M
    text(real(mapping(i))-l/8,imag(mapping(i))+l/8,num2str(de2bi(i-1,k,'left-msb')), 'FontSize', 8);
end