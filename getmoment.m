function out = getmoment(Y,f,n)

%     N=length(Y);
%    S=zeros(N,1);
%     for i=1:N
%         S(i)=(2*pi*f(i))^n*Y(i)^2;
%     end
%     out=2*trapz(f,S);

 
        S=bsxfun(@times,(2*pi*f).^n,Y.^2);
   
    out=2*trapz(f,S');
