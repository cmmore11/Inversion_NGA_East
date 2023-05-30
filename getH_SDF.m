function absH=getH_SDF(Damp,fn,freq)
% PURPOSE
%   get transfer function of a SDOF oscilator
%   acc response to acc excitation
% INPUT
%   Damp - damping ratio
%   fn - fundamental frequency of the SDOF (Hz)
%   freq - frequency array (Hz)
% OUTPUT
%   absH - Amplitude of the transfer function

% N=length(freq);
% H=zeros(N,1);
% for i=1:N
%     f=freq(i);
%     H(i)=(-fn^2)/((f^2-fn^2)-(2*f*fn*D*1i));
% end
% absH=abs(H);



 
    H =bsxfun(@rdivide,(-fn.^2),(bsxfun(@minus,freq.^2,fn.^2)-(2*bsxfun(@times,freq,fn)*Damp*1i)));


absH=abs(H);