function [z,PSA_Stochastic,exceed] =vec_fit_alt(x)
 %*ones(1,25)]

%x = mapvars(x);
% Rrup1 = [1,45,50];
% Mag = [5,7];
% Str_Frequency = [0.5,5,20];
Rrup1 = logspace(log10(1),log10(1000),30);
Str_Frequency = logspace(-1,log10(100),25);
Mag = 4;
%[4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2];

[MRF, MedianGMPEs, ~] = CalculateMedianNGA(Mag,Rrup1,Str_Frequency);
MedianGMPEs = reshape(MedianGMPEs,length(Rrup1)*length(Str_Frequency)*length(Mag),1);
AF1 = [0.001 1;0.008 1.003;0.023 1.01;0.040 1.017;0.061 1.026;
    0.108 1.047;0.234 1.069;0.345 1.084;0.508 1.101;1.090 1.135;
    1.370 1.143;1.690 1.148;1.970 1.150;2.420 1.151;101 1.15101]; %frequency, amp factor
%MedMinusSig = exp(log(MedianGMPEs)-log(1*Sigmas));
%MedPlusSig = exp(log(MedianGMPEs)+log(1*Sigmas));
Ro= 2.8; %2.72;%
Beta_s= 3.7;

Damp=0.05;

npionts=1024;
start_freq=0.1;
end_freq=100;

frinc=log10(end_freq/start_freq)/(npionts-1);
freq=start_freq*10.^((0:1:(npionts-1))*frinc);
AF = interp1(AF1(:,1),AF1(:,2),freq);
%round(x(:,1).*100/5)*5
stressdrop=round(x(:,1).*100/5)*5;  %repmat(x(:,1),1,length(unique(MRF(:,3))));
%stressdrop=repmat(stressdrop,1,length(unique(MRF(:,2))));
f_c=4.906*10^6*Beta_s*(stressdrop'./(10.^(1.5*MRF(:,1)+16.05))).^(1/3)';    % Corner frequency

Kappa=0.006.*ones(1,length(MRF(:,1)));%[Kappa1(indM1),Kappa2(indM2),Kappa3(indM3)];

Q = (x(:,13).*1000)*(freq.^(x(:,12)./10)); %bsxfun(@times,440,(bsxfun(@power,freq,0.47))); %

h = x(:,2).*10;

R=sqrt(MRF(:,2).^2+h.^2);

R1=50;
R2=125;
indR1=find(R<=R1);
indR2=find(R>R1 & R<=R2);
indR3=find(R>R2);

f = [Str_Frequency(1);Str_Frequency(4);Str_Frequency(7);Str_Frequency(10);Str_Frequency(13);Str_Frequency(17);Str_Frequency(21);Str_Frequency(25)];
%hz
b1in = [0.1 x(:,3)
    f(2) x(:,4)
    f(3) x(:,5)
    f(4) x(:,6)
    f(5) x(:,7)
    f(6) x(:,8)
    f(7) x(:,9)
    100 x(:,10)];
b1a = interp1(b1in(:,1),b1in(:,2),MRF(1:25,3),'makima');
b1 = repmat(b1a,length(Rrup1),1)';
%b11 = x(:,5)+x(:,6).*log10(f1)*ones(1,length(MRF(:,3))));
%disp(size(b11))
% b12= x(:,7);
% %disp(size(b12))
% b13 = x(:,8);
% %disp(size(b13))
% b1 = [b11(indf1),b12(indf2),b13(indf3)];
%disp(size(b1))
%GMSpredg = x(:,13);
%GMSpredg=repmat(GMSpredg,1,length(unique(MRF(:,2)))*length(unique(MRF(:,1))));
b2 = x(:,11);
b3 = -0.5;
Z1=R'.^b1;
Z2=(R1.^b1).*(R'/R1).^(b2);
Z3=((R1.^b1).*((R2/R1).^(b2))).*(R'/R2).^(b3);%
Z=[Z1(indR1),Z2(indR2),Z3(indR3)];

FAS=bsxfun(@times,...
    bsxfun(@rdivide,(10^-20*1/100*1/9.81*(0.781*pi/(Ro*Beta_s^3))*10.^(1.5*MRF(:,1)+16.05)*freq.^2),...
    (1+(bsxfun(@rdivide,freq,f_c').^2))),Z').*...
    bsxfun(@times,exp(bsxfun(@rdivide,-1*pi*R,(Q./freq*Beta_s))),...
    bsxfun(@times,exp(bsxfun(@times,-1*pi*Kappa',freq)),AF));

BTPATH=[0 0
    15 2.6
    35 17.5
    50 25.1
    125 25.1
    200 28.5
    392 46.0
    600 69.1];

Tgm1=interp1(BTPATH(:,1),BTPATH(:,2),R)+1./f_c';
Tgm2=69.1+0.111*(R-600)+1./f_c'; %EDITTED 7/14/2021
Tgm=[Tgm1(R<=600)',Tgm2(R>600)'];

PSA_Stochastic=FAS2SA(FAS,freq,Tgm,MRF(:,3), Damp)';

e1 = abs((PSA_Stochastic)-(MedianGMPEs))./abs(MedianGMPEs);
zzz = abs(log10(PSA_Stochastic)-log10(MedianGMPEs));
zz = zzz.^2;
z = sum(zz);
z = sqrt(sum(z));
% zzz = sum(e1);
% zz = sum(zzz);
% z = zz/length(PSA_Stochastic);

exceed = 0;
for i = 1:length(e1)
    if e1(i) >= 0.1
        exceed = exceed+1;
    end
end

end