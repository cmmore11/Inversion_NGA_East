function [SA]=FAS2SA(FAS,F,T_gm,freq, damping)

% INPUT
%   FAS  - Fourier Amplitude Spectrum
%   F    - Frquency array corresponding FAS (Hz)
%   T_gm - ground motion duration (sec)
%   freq - fundametnal frequency of SDOF oscilator (Hz)
%   damping - damping ratio of SDOF oscilator
% OUTPUT
%   SA   - Spectral acceleration (g)
%   PF   - Peak Factor 

  
    
    H=getH_SDF(damping,freq,F);  % get transfer function of SDOF oscilator
    
    Y=FAS.*H;  % FAS of acc response of oscilator
    Tn=1./freq ;  % fundamantal period
    T0=Tn/2/pi/damping;  %oscilator duration
    gamma=T_gm./Tn';  
  
    n_val=3; alpha=1./3.;
   
    T_rms=T_gm+T0'.*(gamma.^n_val./(gamma.^n_val+alpha)); % total duration
	
	
    moment0=getmoment(Y,F,0);
    moment2=getmoment(Y,F,2);
    moment4=getmoment(Y,F,4);
	
    Y_rms=(moment0./T_rms).^0.5;  % rms value
	
    Ne=1/pi*(moment4./moment2).^0.5.*T_gm;  % number of extrema
    % note that here use T_gm not T_rms (Refer to Boore and Joyer, 1984)
    
    xi=(moment2.^2./moment0./moment4).^0.5;  % number of width
    temp=(2*log(xi.*Ne)).^0.5;
    PF=temp+0.5772./temp;
	
% 	Fun=@(x) 1-(1-xi.*exp(-x.^2)).^Ne;
%     PF1 = 2^0.5*quadgk(Fun,0,Inf); 
    
    
    
    Ymax=PF.*Y_rms;   
   
   SA =Ymax;
	
 