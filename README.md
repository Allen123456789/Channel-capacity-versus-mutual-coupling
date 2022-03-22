# MATLAB code: Channel-capacity-versus-mutual-coupling
% Reference
% https://www.researchgate.net/post/How-to-write-Monte-Carlo-simulation-code-in-Matlab-for-finding-outage-capacity-of-a-channel
% https://www.mathworks.com/help/antenna/ug/effect-of-mutual-coupling-on-mimo-communication.html#EffectOfMutualCouplingOnMIMOCommunicationExample-2
% https://www.mathworks.com/matlabcentral/answers/364786-where-i-can-find-helpercalculatecouplingmatrix-m-function


%2.14GHz Z matrix and S matrix from HFSS 
clear;
clc;

% couple_IFA Z matrix and S matrix antenna efficiency when f=2.14GHz
%Z=[48.18+i*(19.67) 28.88+i*(-2.69);28.88+i*(-2.69) 47.21+i*(19.68) ]; 
%S_11_linear=10^((-11.35)/20);
%S_22_linear=10^((-11.26)/20);
%S_21_linear=10^((-10.03)/20);
%S_12_linear=10^((-10.03)/20);
%total_efficiency1=0.87;
%total_efficiency2=0.87;

% self_decouple_h_1_6_C_1_348_1_567pF Z matrix and S matrix when f=2.14GHz
%Z=[51.81+i*(-8.42) 11.49+i*(-1.69);11.49+i*(-1.69) 33.35+i*(20.55) ];
%S_11_linear=10^((-22.41)/20);
%S_22_linear=10^((-9.75)/20);
%S_21_linear=10^((-17.43)/20);
%S_12_linear=10^((-17.43)/20);
%total_efficiency1=0.9;
%total_efficiency2=0.85;

% self_decouple_h_1_6_C_1_348_1_567pF Z matrix and S matrix when f=2.14GHz
Z=[51.81+i*(-8.42) 11.49+i*(-1.69);11.49+i*(-1.69) 33.35+i*(20.55) ];
S_11_linear=10^((-22.41)/20);
S_22_linear=10^((-9.75)/20);
S_21_linear=10^((-17.43)/20);
S_12_linear=10^((-17.43)/20);
S=[S_11_linear S_12_linear; S_21_linear S_22_linear];
total_efficiency1=0.9;
total_efficiency2=0.85;

% self_decouple_best_h_1_56_C_1_3_1_4pF Z matrix and S matrix antenna efficiency when f=2.14GHz
%Z=[56.54+i*(5.46) 3.36+i*(-0.14);3.36+i*(-0.14) 28.79+i*(-8.86) ];
%S_11_linear=10^((-21.99)/20);
%S_22_linear=10^((-10.72)/20);
%S_21_linear=10^((-28.53)/20);
%S_12_linear=10^((-28.53)/20);
%total_efficiency1=0.92;
%total_efficiency2=0.88;



total_efficiency=[total_efficiency1 0; 0 total_efficiency2];

A=Z; %transmit Z matrix
D=Z; %receive Z matrix
ZL=[50 0; 0 50];  %define ZL at each port 50 ohm

%calculate coupling matrix from the Z matrix
CRX=ZL*((ZL+D)\eye(2));
KRX=[0 CRX(1,2); CRX(2,1) 0]*(diag(diag(CRX))\eye(2));
CTX=A\eye(2);
KTX=diag(diag(CTX))\[0 CTX(1,2);CTX(2,1) 0];

%By "MIMO Capacity and Antenna Array Design" (26) (27) if the below is
%small than 0 the antenna coupling is harmful for the channel capacity
log2(det((eye(2)+KTX)*(eye(2)+KTX)'))+log2(det((eye(2)+KRX)*(eye(2)+KRX)')) 

%calculate ECC and from the S matrix this is a undervalue compared to the
%value calculated from HFSS (HFSS-toolkit-MIMO)
%term1 = conj(S(1,1)).*S(1,2);
%term2 = conj(S(2,1)).*S(2,2);
%numerator = abs(term1 + term2).^2;
%denominator = (1 - (abs(S(1,1)).^2 + abs(S(1,2)).^2)).*(1 - (abs(S(2,2)).^2 + abs(S(2,1)).^2));
%ECC_by_Sparameter=numerator./denominator;
%r_by_Sparameter = sqrt(numerator./denominator);

% couple_IFA Z matrix and S matrixwhen f=2.14GHz
%r_by_HFSS=0.039;   %phase?

% self_decouple_h_1_6_C_1_348_1_567pF Z matrix and S matrix when f=2.14GHz
r_by_HFSS=0.0039;   %phase?

% self_decouple_best_h_1_56_C_1_3_1_4pF Z matrix and S matrix when f=2.14GHz
%r_by_HFSS=0.000006;   %phase?



%先隨機生成一個phase
rand_phase=360*rand(1,1)/180*pi;

r_by_HFSS_random_phase=r_by_HFSS*(cos(rand_phase)+1i*sin(rand_phase));

%monti carlo simulation for channel capacity
SNR_dB=[0:30];
SNR_linear=10.^(SNR_dB/10);
nT=2; nR=2; % nTxnR Tx and Rx antennas
n=min(nT,nR);% determine the rank of channel matrix
I = eye(n);
N_trial= 50000*n; % Monte Carlo (Number of realization)
Ciid = zeros(N_trial,length(SNR_dB));
Cnc = zeros(N_trial,length(SNR_dB));
Cc = zeros(N_trial,length(SNR_dB));

R_bar=[1,r_by_HFSS_random_phase;conj(r_by_HFSS_random_phase),1];%arbitary receieve correlation matrix but it should be symmetrical matrix
[V,Diag,W] = eig(total_efficiency);
Diag_sqrt=Diag.^0.5;
total_efficiency_sqrt=V*Diag_sqrt*W;

R=total_efficiency_sqrt*R_bar*total_efficiency_sqrt;
[V,Diag,W] = eig(R);
Diag_sqrt=Diag.^0.5;
R_sqrt=V*Diag_sqrt*W;


for k=1:N_trial
Hw_real = randn(nT,nR);
Hw_img=randn(nT,nR);
var=1; %normal distribution的變異數為1 若不是normal distribution var要改 因為你另外有控制SNR，所以這裡VAR要設為1
Hw=sqrt(var/2)*(Hw_real+j*Hw_img);
HwHw=Hw*Hw';
Hnc=R_sqrt*Hw;
HncHnc=Hnc*Hnc';  % eigenvalue,phase
Hc=(eye(2)+KRX)*Hnc*(eye(2)+KTX);
HcHc=Hc*Hc';
for l=1:length(SNR_dB)
Ciid(k,l) = log2(real(det(I+(SNR_linear(l)/nT)*HwHw)));
Cnc(k,l) = log2(real(det(I+(SNR_linear(l)/nT)*HncHnc)));
Cc(k,l) = log2(real(det(I+(SNR_linear(l)/nT)*HcHc)));
end
end
Ciid_m=mean(Ciid);
Cnc_m=mean(Cnc);
Cc_m=mean(Cc);
% Plots
figure(1)
plot(SNR_dB(15:31),Ciid_m(15:31),'r-o', 'LineWidth',5);
hold on;
%plot(SNR_dB,Cnc_m,'g-x');
plot(SNR_dB(15:31),Cc_m(15:31),'b-*','LineWidth',5);


legend('C_{Rayleigh fading channel}','C_{Antenna coupling}','Location','NorthWest')
xlabel('SNR [dB]'); ylabel('Channel capccity [bps/Hz]');
xlim([15 30])
ylim([6 18])
grid on;
title('MIMO channel capacity nT=nR=2');
set(gca,'FontSize',15)

