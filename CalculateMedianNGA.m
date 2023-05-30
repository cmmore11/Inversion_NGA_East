%%%%%%%%%%%%%
%M = magnitude
%Rrup = rupture distance (km)
%Structural frequency = frequencies (hz)
%change these three items to output the median from the NGA East GMM tool
%and develop the MRF matrix [magnitudes, distances, frequencies]
%%%%%%%%%%%%%
function [MRF, Medians, Sigmas] = CalculateMedianNGA(M,Rrup,Structural_Frequency)

counter=0;

%AF_hash_stew = [0.1	0.739910534;0.11	0.749992554;0.112	0.75247237;0.113	0.754205803;0.114	0.755363642;0.115	0.757280161;0.116	0.758619458;
%    0.117	0.760544238;0.118	0.761711808;0.119	0.763644433;0.12	0.765403602;0.125	0.773847077;0.13	0.782201421;0.135	0.789855707;0.14	0.797956655;
%    0.15	0.811047222;0.2	0.838209676;0.25	0.862538926;0.3	0.892889898;0.4	0.904942385;0.5	0.90833139;0.75	0.883568021;0.8	0.88696301;1	0.892898827;1.5	0.900525317;
%    2	0.8949298;3	0.897235419;4	0.889196529;5	0.896164759;7.5	0.911025856;10	0.918655583];

MedianGMMs = zeros(length(M)*length(Rrup)*length(Structural_Frequency),1);
SigmaGMMs = zeros(length(M)*length(Rrup)*length(Structural_Frequency),1);
MRF = zeros(length(M)*length(Rrup)*length(Structural_Frequency),3);
%AF = zeros(length(M)*length(Rrup)*length(Structural_Frequency),3);
for j=1:length(M)
    for i=1:length(Rrup)
        for k=1:length(Structural_Frequency)
            counter = counter+1;
            MRF(counter,1)= M(j);
            MRF(counter,2)= Rrup(i);
            MRF(counter,3)= Structural_Frequency(k);
            [P,S] =NGA_East_GMMtool_Matlab_R2003162(MRF(counter,1), MRF(counter,2),1./MRF(counter,3),"SingleStation",1.0);
            MedianGMMs(counter,1)= P;
            SigmaGMMs(counter,1) = S;
        end
    end
end

Medians = reshape(MedianGMMs,length(Rrup)*length(Structural_Frequency),length(M));
Sigmas = reshape(SigmaGMMs,length(Rrup)*length(Structural_Frequency),length(M));
            
            %AF(i) = interp1(AF_hash_stew(:,1),AF_hash_stew(:,2),1./MRF(i,3));
           
            %     inMagnitude = 5; % Range: 4-8.2
            %     inRrup = 111.2; % range: 0-1500km
            %     inPeriodList = [0.01,0.015,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7.5,10,"PGA","PGV"]; % range:0.01-10sec
            %     inSigmaType = "SingleStation"; % options: "Ergodic" "SingleStation"
            %     inNsigma = 1.0; 
            %[~,~,~,~,~,MedianGMPEs(counter)]=Func_Median_NGA(MRF(:,1), MRF(:,2),1/MRF(:,3));
% save('SigmaGMPEs','SigmaGMPEs')
% save('MedianGMPEs','MedianGMPEs')
% save('MRF','MRF')

end










