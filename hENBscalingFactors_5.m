%hENBscalingFactors Calculate ENB signal scaling factors 
%   [K1, K2, K3] = hENBscalingFactors(DIP2, DIP3, NOC, SINR, ENB1, ENB2, ENB3)
%   calculates the signal scaling factors to achieve the specified SINR,
%   DIP2 and DIP3 in a 3 cell scenario: serving cell and two interferers.
%
%   DIP2, DIP3 and SINR are specified in dB. NOC is specified in dBm/15kHz.
%
%   ENB1 is the ENB structure for the signal of interest, while ENB2 and
%   ENB3 are the ENB structures for the two interferers. These are needed
%   to calculate the power of the interferer signals.
%
%   K1, K2, K3 are the scaling factors to apply to the time domain
%   waveforms from the serving and interfering cells respectively to
%   achieve the specified SINR, DIP2 and DIP3.

%   Copyright 2014 The MathWorks, Inc.

function [K1, K2, K3, K4, K5] = hENBscalingFactors_5(DIP2, DIP3, DIP4, DIP5, NOC, SINR, enb1, enb2, enb3, enb4, enb5)

% DIP and Noc in linear scales    
dip2 = 10.^(DIP2/10);
dip3 = 10.^(DIP3/10);
dip4 = 10.^(DIP4/10);
dip5 = 10.^(DIP5/10);
Noc = 10.^((NOC-30)/10); % Convert to W/15kHz
            
% Average received power spectral densities for cells 2 and 3
% Solve linear system of equations
% dip2 = Ior2 / (Ior2 + Ior3 + Noc);
% dip3 = Ior3 / (Ior2 + Ior3 + Noc);
% Rewrite first in the form Ax = B with x = [Ior2; Ior3]
A = [dip2-1 dip2 dip2 dip2; ...
     dip3 dip3-1 dip3 dip3; ...
     dip4 dip4 dip4-1 dip4; ...
     dip5 dip5 dip5 dip5-1];
B = [-dip2*Noc; -dip3*Noc; -dip4*Noc; -dip5*Noc];
x = mldivide(A,B); % Solve Ax = B
Ior2 = x(1); % W/15kHz             
Ior3 = x(2); % W/15kHz 
Ior4 = x(3);
Ior5 = x(4);

% SINR = EsTot / NocTot = (Es1 + Es2) / (Noc1 + Noc2) 36.101 Sec. 8.1.1
NocTot = Ior2 + Ior3 + Ior4 + Ior5 + Noc; % W/15kHz 
sinr = 10.^((SINR-enb1.PDSCH.Rho)/10); % linear, compensate for the effect of Rho
EsTot = sinr*NocTot;        % W/15kHz 


% Amplitude scaling factors
% For Ior2 and Ior3 normalise by the power of the interfering waveforms. We
% could have approximated these and use the following scaling factor:
%   K2 = sqrt(Ior2*NoTxAnts/(10^(Rho/10)));
% However that's an approximation to the expected power of the interfering
% signal. Here we measure this power for 1 frame.
interfPower2 = powerETM3Interferer(enb2);
interfPower3 = powerETM3Interferer(enb3);
interfPower4 = powerETM3Interferer(enb4);
interfPower5 = powerETM3Interferer(enb5);

K2 = sqrt(Ior2/interfPower2);
K3 = sqrt(Ior3/interfPower3);
K4 = sqrt(Ior4/interfPower4);
K5 = sqrt(Ior5/interfPower5);

% K1 is not normalised because:
% - Rho has been taken into account in sinr = 10.^((SINR-enb1.PDSCH.Rho)/10);
% - The number of antennas is not necessary here, in K2 and K3 it is needed
%   because it affects the power of the transmitted signal (TxDiversity),
%   while EsTot is measured without taking into account the effect of
%   precoding (see Section 8.1.1 TS 36.101). Note that this is a difference
%   in the definition of Ior vs EsTot.
K1 = sqrt(EsTot);

end

function signalPower = powerETM3Interferer(enb)
% Measure the power of ETM3 interferer over 1 frame

nSubframes = 10; % number of subframes used to measure the power

% Generate required number of subframes (starting at 0) and measure the power
enb.TotSubframes = nSubframes;
enb.NSubframe = 0;
waveform = hTM3InterfModel(enb);
ofdmInfo = lteOFDMInfo(enb);
signalPower = (var(waveform)*double(ofdmInfo.Nfft)^2/(enb.NDLRB*12));

% Average from all antennas
signalPower = mean(signalPower);

end
 