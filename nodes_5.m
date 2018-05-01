 NFrames = 12;   % Number of frames

SINR = -1.4;    % SINR value in dB
DIP2 = -1.73;   % DIP in dB for cell 2
DIP3 = -8.66;   % DIP in dB for cell 3
DIP4 = -1.73;   % DIP in dB for cell 4
DIP5 = -8.66;   % DIP in dB for cell 5

Noc = -98; % dBm/15kHz average power spectral density\

% Set the random number generator seed
rng('default');

% Cell 1 eNodeB configuration according to R.46
enb1 = struct;
%enb1 = lteRMCDL('R.46'); 
enb1.NDLRB = 50;
enb1.CellRefP = 2;
enb1.NCellID = 0; 
enb1.CyclicPrefix = 'Normal';
enb1.CFI = 2;
enb1.PHICHDuration = 'Normal';
enb1.Ng = 'Sixth';
enb1.NFrame = 0;
enb1.TotSubframes = 1;
enb1.Windowing = 0;
enb1.DuplexMode = 'TDD';
enb1.SSC = 4;
enb1.TDDConfig = 1;
enb1.NPRSRB = 2;  
enb1.IPRS = 0; 
enb1.PRSPeriod = 'On'; 
enb1.Position = [2000,3000];
      

% PDSCH configuration substructure
enb1.PDSCH.TxScheme = 'TxDiversity'; % PDSCH transmission mode 2
enb1.PDSCH.Modulation = {'QPSK'};
enb1.PDSCH.NLayers = 2;
enb1.PDSCH.Rho = -3;
enb1.PDSCH.RNTI = 1;
enb1.PDSCH.RVSeq = [0 1 2 3];
enb1.PDSCH.RV = 0;
enb1.PDSCH.NHARQProcesses = 7;
enb1.PDSCH.NTurboDecIts = 5;
enb1.PDSCH.PRBSet = (0:49)';
% Table A.3.4.2.1-2, TS 36.101
enb1.PDSCH.TrBlkSizes = [5160 3880 0 0 5160 0 3880 0 0 5160];
% Table A.3.4.2.1-2, TS 36.101
enb1.PDSCH.CodedTrBlkSizes = [12528 10656 0 0 13200 0 10656 0 0 13200];
enb1.PDSCH.CSIMode = 'PUCCH 1-0';
enb1.PDSCH.PMIMode = 'Wideband';
enb1.PDSCH.CSI = 'On';

% PDSCH OCNG configuration
enb1.OCNGPDSCHEnable = 'On';             % Enable OCNG fill
enb1.OCNGPDSCHPower = -3;                % OCNG power same as PDSCH Rho
enb1.OCNGPDSCH.RNTI = 0;                 % Virtual UE RNTI
enb1.OCNGPDSCH.Modulation = 'QPSK';      % OCNG symbol modulation
enb1.OCNGPDSCH.TxScheme = 'TxDiversity'; % OCNG transmission mode 2\

info1 = lteOFDMInfo(enb1);

% Cell2
enb2 = enb1;
enb2.NCellID = 1;
enb2.Position = [2000,3500];
enb2.OCNGPDSCHEnable = 'Off';

% Cell 3
enb3 = enb1;
enb3.NCellID = 2;
enb3.Position = [2500,4000];
enb3.OCNGPDSCHEnable = 'Off';

% Cell 3
enb4 = enb1;
enb4.NCellID = 3;
enb4.Position = [2500,4500];
enb4.OCNGPDSCHEnable = 'Off';

% Cell 4
enb5 = enb1;
enb5.NCellID = 4;
enb5.Position = [3000,4500];
enb5.OCNGPDSCHEnable = 'Off';

hPositioningPlotPositions({enb1, enb2, enb3, enb4, enb5});

    
    
% eNodeB1 to UE propagation channel
channel1 = struct;                    % Channel config structure
channel1.Seed = 20;                   % Random channel seed
channel1.NRxAnts = 8;                 % 2 receive antennas
channel1.DelayProfile ='EVA';         % Delay profile
channel1.DopplerFreq = 5;            % Doppler frequency
channel1.MIMOCorrelation = 'Low';     % Multi-antenna correlation
channel1.NTerms = 16;                 % Oscillators used in fading model
channel1.ModelType = 'GMEDS';         % Rayleigh fading model type
channel1.InitPhase = 'Random';        % Random initial phases
channel1.NormalizePathGains = 'On';   % Normalize delay profile power
channel1.NormalizeTxAnts = 'On';      % Normalize for transmit antennas
% The channel sampling rate depends on the FFT size used in the OFDM
% modulator. This can be obtained using the function lteOFDMInfo.
ofdmInfo = lteOFDMInfo(enb1);
channel1.SamplingRate = ofdmInfo.SamplingRate;

% eNodeB2 (interference) to UE propagation channel
channel2 = channel1;
channel2.Seed = 122;                  % Random channel seed

% eNodeB3 (interference) to UE propagation channel
channel3 = channel1;
channel3.Seed = 36;                   % Random channel seed

% eNodeB4 (interference) to UE propagation channel
channel4 = channel1;
channel4.Seed = 120;                  % Random channel seed

% eNodeB5 (interference) to UE propagation channel
channel5 = channel1;
channel5.Seed = 40;                  % Random channel seed

% Uplink %

numSubframes = 10;  % Number of frames to simulate at each SNR
SNRdB = [-1.4]; % SNR points to simulate
frc.TotSubframes = 1;   % Total number of subframes to generate
frc.NCellID = 5;       % Cell identity
frc.RC = 'A4-3';        % FRC number

% Populate FRC configuration structure with default values for A4-3
frc = lteRMCUL(frc);

enb1.PUSCH = frc.PUSCH;

frc.PDSCH = enb1.PDSCH;

disp(frc.PUSCH);
% frc.Position = hPositioningPosition(0, 0);
chcfg.NRxAnts = 8;               % Number of receive antennas
chcfg.DelayProfile = 'ETU';      % Delay profile
chcfg.DopplerFreq = 5;          % Doppler frequency
chcfg.MIMOCorrelation = 'Low';   % MIMO correlation
chcfg.Seed = 91;                 % Channel seed
chcfg.NTerms = 16;               % Oscillators used in fading model
chcfg.ModelType = 'GMEDS';       % Rayleigh fading model type
chcfg.InitPhase = 'Random';      % Random initial phases
chcfg.NormalizePathGains = 'On'; % Normalize delay profile power
chcfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas

% Set channel model sampling rate
info = lteSCFDMAInfo(frc);
info.SamplingRate = info1.SamplingRate;
chcfg.SamplingRate = info.SamplingRate;


% End uplink setup %

cec = struct;                        % Channel estimation config structure
cec.PilotAverage = 'UserDefined';    % Type of pilot symbol averaging
cec.FreqWindow = 31;                 % Frequency window size
cec.TimeWindow = 23;                 % Time window size
cec.InterpType = 'Cubic';            % 2D interpolation type
cec.InterpWindow = 'Centered';       % Interpolation window type
cec.InterpWinSize = 1;               % Interpolation window size

% Channel noise setup
nocLin = 10.^(Noc/10)*(1e-3); % linear  in Watts
% Take into account FFT (OFDM) scaling
No = sqrt(nocLin/(2*double(ofdmInfo.Nfft)));

% Signal and interference amplitude scaling factors calculation. These
% ensure the SINR and DIP values specified are met
[K1, K2, K3, K4, K5] = hENBscalingFactors_5(DIP2, DIP3, DIP4, DIP5, Noc, SINR, enb1, enb2, enb3, enb4, enb5);

% Initialize state of all HARQ processes
harqProcesses = hInitHARQProcesses(enb1, enb1.PDSCH.NHARQProcesses);
% Initialize HARQ process IDs to 1 as the first non-zero transport block
% will always be transmitted using the first HARQ process. This will be
% updated with the full sequence output by lteRMCDLTool after the first
% call to the function
harqProcessSequence = 1;

% Set up variables for the main loop
lastOffset = 0;       % Initialize overall frame timing offset
frameOffset = 0;      % Initialize frame timing offset
blkCRC = [];          % Block CRC for all considered subframes
bitTput = [];         % Number of successfully received bits per subframe
txedTrBlkSizes = [];  % Number of transmitted bits per subframe
% Vector of total number of bits transmitted calculated for each subframe
runningMaxThPut = [];
% Vector storing the number of successfully received bits calculated for
% each subframe
runningSimThPut = [];

% Obtain the number of transmit antennas.
dims = lteDLResourceGridSize(enb1);
P = dims(3);



grid = [];
for nsf = 0:19
    enb1.NSubframe = mod(nsf,10);
    sfgrid = lteDLResourceGrid(enb1);       % Empty subframe
    sfgrid(ltePRSIndices(enb1)) = ltePRS(enb1);       % PRS REs
    sfgrid(ltePSSIndices(enb1)) = ltePSS(enb1);       % PSS REs
    sfgrid(lteSSSIndices(enb1)) = lteSSS(enb1);       % SSS REs
    sfgrid(lteCellRSIndices(enb1)) = lteCellRS(enb1); % Cell RS REs
    grid = [grid sfgrid]; %#ok<AGROW>
end
enb1.NSubframe = 0;
tx1 = lteOFDMModulate(enb1, grid);        % OFDM modulate


speedOfLight = 299792458.0; % Speed of light in m/s



[~, radius1] = cart2pol(enb1.Position(1), enb1.Position(2));
delay = radius1/speedOfLight;                  % Delay in seconds
sampleDelay1 = round(delay*info.SamplingRate); % Delay in samples



%  sumrx = length(tx1)+sampleDelay1;
%  rx = cell(1);
%     % Urban Macro LOS path loss per TR36.814
%     PLdB = hPositioningPathLoss(radius1, 2.1e9);
%     PL = 10^(PLdB/10);
% 
%     % Add delay, pad and attenuate
%     rx{1} = [zeros(sampleDelay1, 1); zeros(tx1, 1); ...
%                 zeros(0, 1)]/ sqrt(PL);
% 
%     % Sum waveforms from all eNodeBs
%     sumrx = sumrx + rx{1};
% 
% 
% % Plot received waveforms
% hPositioningPlotRx(enb1, rx{1}, info.SamplingRate);
% 
% % Assumed parameters for cell search
% searchcfg.CyclicPrefix = enb1.CyclicPrefix;
% searchcfg.DuplexMode = enb1.DuplexMode;
% searchcfg.NDLRB = enb1.NDLRB;
% 
% % Perform multi-cell search
% searchalg.MaxCellCount = 5;
% searchalg.SSSDetection = 'PostFFT';
% [cellIDs,offsets] = lteCellSearch(searchcfg,sumrx,searchalg);
% 
% % Set up configurations for each detected cell; cells are considered as
% % detected here if they meet a minimum RSRQ threshold Qqualmin
% Qqualmin = -20;
% RSRQdB = zeros(1,searchalg.MaxCellCount);
% rxcfg = cell(1,searchalg.MaxCellCount);
% for i = 1:searchalg.MaxCellCount
%     % Assumed parameters
%     rxcfg{i} = enb1;
%     % Use cell identity that was detected
%     rxcfg{i}.NCellID = cellIDs(i);
%     % Measure RSRQ
%     rxgrid = lteOFDMDemodulate(rxcfg{i},sumrx(1+offsets(i):end,:));
%     meas = hRSMeasurements(rxcfg{i},rxgrid);
%     RSRQdB(i) = meas.RSRQdB;
% end
% rxcfg(RSRQdB<Qqualmin) = [];
% Ndetected = numel(rxcfg);
% 
% 
% ref = cell(1, Ndetected);
% corr = cell(1, Ndetected);
% delayEst = zeros(1, Ndetected);
% for i = 1:Ndetected
%     % Generate reference PRS
%     sfgrid = lteDLResourceGrid(rxcfg{i});
%     sfgrid(ltePRSIndices(rxcfg{i})) = ltePRS(rxcfg{i});
%     ref{i} = lteOFDMModulate(rxcfg{i},sfgrid);
% 
%     % Correlate received signal with each reference PRS
%     c = abs(xcorr(sumrx,ref{i}));
% 
%     % Reduced length of correlation vector for positioning and plotting
%     c(1:length(sumrx)) = [];    % Remove meaningless result at beginning
%     corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
% 
%     % Delay estimate is at point of maximum correlation
%     delayEst(i) = find(corr{i}==max(corr{i}));
% end
% 
% % Plot correlation
% if (Ndetected>0)
%     hPositioningPlotCorr(rxcfg, corr, info.SamplingRate);
% end


% % Estimate time difference of arrival from each eNodeB
% tdoa = hPositioningTDOA(delayEst,info.SamplingRate);
% 
% % Plot hyperbolas
% figure(1);
% legendstr = feval(@(x)x.String,legend);
% enbs = [enb{:}];
% txCellIDs = [enbs.NCellID];
% for j = 1:Ndetected
%     for i = (j+1):Ndetected
%         dd = tdoa(i,j)*speedOfLight; % Delay distance
%         % establish the eNodeBs for which the delay distance
%         % is applicable by examining the detected cell identities
%         txi = find(txCellIDs==rxcfg{i}.NCellID);
%         txj = find(txCellIDs==rxcfg{j}.NCellID);
%         if (~isempty(txi) && ~isempty(txj))
%             % plot TDOA curve
%             [x, y] = hPositioningTDOACurve(enb{txi}.Position, ...
%                 enb{txj}.Position, dd);
%             plot(x, y, 'k:', 'LineWidth', 2);
%         end
%     end
% end
% legend(legendstr);


fprintf('\nSimulating %d frame(s)\n',NFrames);

% Main for loop: for all subframes
for subframeNo = 0:(NFrames*10-1)

    % Reinitialize channel seed for each subframe to increase variability
    channel1.Seed = 1+subframeNo;
    channel2.Seed = 1+subframeNo+(NFrames*10);
    channel3.Seed = 1+subframeNo+2*(NFrames*10);
    channel4.Seed = 1+subframeNo+3*(NFrames*10);
    channel5.Seed = 1+subframeNo+4*(NFrames*10);
    
    % Update subframe number
%     frc.NSubframe = subframeNo;
    enb1.NSubframe = subframeNo;
    enb2.NSubframe = subframeNo;
    enb3.NSubframe = subframeNo;
    enb4.NSubframe = subframeNo;
    enb5.NSubframe = subframeNo;
    
    duplexInfo = lteDuplexingInfo(enb1);
    duplexInfo2 = lteDuplexingInfo(enb1);

    if duplexInfo.NSymbolsDL ~= 0 % target only downlink subframes
%     if duplexInfo2.NSymbolsDL ~= 0 

        % Get HARQ process ID for the subframe from HARQ process sequence
        harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);

        % If there is a transport block scheduled in the current subframe
        % (indicated by a non-zero 'harqID'), perform transmission and
        % reception. Otherwise, continue to the next subframe.
        if harqID == 0
            continue;
        end

        % Update HARQ process: reset for successful transmission or set to
        % retransmit: This effectively gets the ACK/NACK
        harqProcesses(harqID) = hTxDiversityHARQScheduling( ...
            harqProcesses(harqID));

        % Update the RV value and RVSeq for correct waveform generation
        enb1.PDSCH.RV = harqProcesses(harqID).rvSeq ...
            (harqProcesses(harqID).rvIdx);
        enb1.PDSCH.RVSeq = harqProcesses(harqID).rvSeq ...
            (harqProcesses(harqID).rvIdx);

%         frc.PDSCH.RV = harqProcesses(harqID).rvSeq ...
%             (harqProcesses(harqID).rvIdx);
%         frc.PDSCH.RVSeq = harqProcesses(harqID).rvSeq ...
%             (harqProcesses(harqID).rvIdx);

        % Serving cell payload
        dlschTransportBlk = {harqProcesses(harqID).dlschTransportBlk};

        % Create transmit waveform and get the HARQ scheduling ID sequence
        % from 'enbOut' structure output which also contains the waveform
        % configuration and OFDM modulation parameters
        [txWaveform1,~,enbOut] = lteRMCDLTool(enb1,dlschTransportBlk);
%         [txWaveform1,~,ueOut] = lteRMCDLTool(frc,dlschTransportBlk);

        % Add 25 sample padding. This is to cover the range of delays
        % expected from channel modeling (a combination of implementation
        % delay and channel delay spread)
        txWaveform1 =  [txWaveform1; zeros(25, P)]; %#ok<AGROW>

        % Get the HARQ ID sequence from 'enbOut' for HARQ processing
        harqProcessSequence = enbOut.PDSCH.HARQProcessSequence;

        % Generate interferer model as per as per TS 36.101, B.5.2. The
        % function hTM3InterfModel generates the interferer transmit
        % signal. It internally sets transmission and modulation
        % schemes and generates the required payload bits.
        txWaveform2 = [hTM3InterfModel(enb2); zeros(25,P)];
        txWaveform3 = [hTM3InterfModel(enb3); zeros(25,P)];
        txWaveform4 = [hTM3InterfModel(enb4); zeros(25,P)];
        txWaveform5 = [hTM3InterfModel(enb5); zeros(25,P)];
        
        % Channel time for the present subframe
        channel1.InitTime = subframeNo/1000;
        channel2.InitTime = channel1.InitTime;
        channel3.InitTime = channel1.InitTime;
        channel4.InitTime = channel1.InitTime;
        channel5.InitTime = channel1.InitTime;
        
        % Filter the transmitted waveforms
        rxWaveform1 = lteFadingChannel(channel1,txWaveform1);
        rxWaveform2 = lteFadingChannel(channel2,txWaveform2);
        rxWaveform3 = lteFadingChannel(channel3,txWaveform3);
        rxWaveform4 = lteFadingChannel(channel4,txWaveform4);
        rxWaveform5 = lteFadingChannel(channel5,txWaveform5);
        % Noise generation
        noise = No*complex(randn(size(rxWaveform1)),...
                randn(size(rxWaveform1)));

        % Add AWGN to the received time domain waveform and scale for
        % required power
        rxWaveform = K1*rxWaveform1 + K2*rxWaveform2 + ...
                     K3*rxWaveform3 + noise + K4*rxWaveform5 + K5*rxWaveform5;

        % Receiver
        % Perform synchronization
        % An offset within the range of delays expected from the
        % channel modeling (a combination of implementation delay and
        % channel delay spread) indicates success
        if (mod(subframeNo,10)==0)
            frameOffset = lteDLFrameOffset(enb1,rxWaveform);
            if (frameOffset > 25)
                frameOffset = lastOffset;
            end
            lastOffset = frameOffset;
        end
        rxWaveform = rxWaveform(1+frameOffset:end,:);

%         if (mod(subframeNo,10)==0)
%             frameOffset = lteDLFrameOffset(frc,rxWaveform);
%             if (frameOffset > 25)
%                 frameOffset = lastOffset;
%             end
%             lastOffset = frameOffset;
%         end
%         rxWaveform = rxWaveform(1+frameOffset:end,:);


        % Scale rxWaveform by 1/K1 to avoid numerical issues with
        % channel decoding stages
        rxWaveform = (1/K1)*rxWaveform;

        % Perform OFDM demodulation on the received data to obtain
        % the resource grid
        rxSubframe = lteOFDMDemodulate(enb1,rxWaveform);

%         rxSubframe = lteOFDMDemodulate(frc,rxWaveform);

        % Channel estimation
        [estChannelGrid,noiseEst] = lteDLChannelEstimate(enb1,cec, ...
                                    rxSubframe);

%     [estChannelGrid,noiseEst] = lteDLChannelEstimate(frc,cec, ...
%                                     rxSubframe);

        % Perform equalization, deprecoding, layer demapping,
        % demodulation and descrambling on the received data using the
        % channel estimate.

        % Get PDSCH indices
        
        [pdschIndices,~] = ...
            ltePDSCHIndices(enb1, enb1.PDSCH, enb1.PDSCH.PRBSet);


%         [pdschIndices,~] = ...
%             ltePDSCHIndices(frc, frc.PDSCH, frc.PDSCH.PRBSet);

        % Get PDSCH resource elements. Scale the received subframe by
        % the PDSCH power factor Rho. The PDSCH is scaled by this
        % amount, while the cell reference symbols used for channel
        % estimation (used in the PDSCH decoding stage) are not.
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
            rxSubframe*(10^(-enb1.PDSCH.Rho/20)), estChannelGrid);
        
%         [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
%             rxSubframe*(10^(-frc.PDSCH.Rho/20)), estChannelGrid);

        % Decode PDSCH
        [dlschBits,~] = ltePDSCHDecode(enb1, enb1.PDSCH,...
            pdschRx, pdschHest, noiseEst);
        
%         [dlschBits,~] = ltePDSCHDecode(frc, frc.PDSCH,...
%             pdschRx, pdschHest, noiseEst);

        % Decode downlink shared channel (DL-SCH)
        [decbits,harqProcesses(harqID).crc,harqProcesses(harqID).decState] = ...
            lteDLSCHDecode(enb1,enb1.PDSCH,harqProcesses(harqID).trBlkSize, ...
            dlschBits{1},harqProcesses(harqID).decState);
        
%         [decbits,harqProcesses(harqID).crc,harqProcesses(harqID).decState] = ...
%             lteDLSCHDecode(frc,frc.PDSCH,harqProcesses(harqID).trBlkSize, ...
%             dlschBits{1},harqProcesses(harqID).decState);

        % Store values to calculate throughput
        blkCRC = [blkCRC harqProcesses(harqID).crc]; %#ok<AGROW>
        bitTput = [bitTput harqProcesses(harqID).trBlkSize.*(1- ...
            harqProcesses(harqID).crc)]; %#ok<AGROW>
        runningSimThPut = [runningSimThPut sum(bitTput,2)]; %#ok<AGROW>
        txedTrBlkSizes = [txedTrBlkSizes harqProcesses(harqID).trBlkSize]; %#ok<AGROW>
        runningMaxThPut = [runningMaxThPut sum(txedTrBlkSizes,2)]; %#ok<AGROW>

    end
end

maxThroughput = sum(txedTrBlkSizes); % Maximum possible throughput
simThroughput = sum(bitTput,2);      % Simulated throughput

% Display noise
% disp(noise)

% Display SNR
%rxWaveform = rxWaveform .* rxWaveform;
%SINR1 = snr(rxWaveform, noise);
%disp(SINR1);

% Display achieved throughput percentage
disp(['Achieved throughput ' num2str(simThroughput*100/maxThroughput) '%'])

% Plot running throughput
figure;plot(runningSimThPut*100./runningMaxThPut)
ylabel('Throughput (%)');
xlabel('Simulated subframe');
title('Throughput');

% Uplink run% 
% pFalse = zeros(size(SNRdB));  % Probability of false detection at each SNR
% pMissed = zeros(size(SNRdB)); % Probability of missed detection at each SNR
% 
% for nSNR = 1:length(SNRdB)
% 
%     % Initialize the random number generator stream
%     rng('default');
% 
%     % Extract SNR to test
%     SNR = 10^(SNRdB(nSNR)/20);
% 
%     % Scale noise to ensure the desired SNR after SC-FDMA demodulation
%     N = 1/(SNR*sqrt(double(info.Nfft)))/sqrt(2.0);
% 
%     offsetUsed = 0;
%     falseCount = 0;
%     falseTests = 0;
%     missCount = 0;
%     missTests = 0;
% 
%     for subframeNo = 0:(numSubframes-1)
% 
%         % Updating subframe number
%         frc.NSubframe = mod(subframeNo, 10);
% 
%         % Transmit ACK on every odd subframe
%         if (mod(subframeNo, 2)==0)
%             ACK = [];
%             falseTests = falseTests + 1;
%         else
%             ACK = 1;
%             missTests = missTests + 1;
%         end
% 
%         % Create random data to transmit
%         trblklen = frc.PUSCH.TrBlkSizes(frc.NSubframe+1);
%         trblk = randi([0 1], trblklen, 1);
% 
%         % Transmit a waveform with an additional 25 samples to cover the
%         % range of delays expected from the channel modeling (a
%         % combination of implementation delay and channel delay spread)
%         txWaveform = [lteRMCULTool(frc, trblk, [], [], ACK); zeros(25, 1)];
% 
%         % Pass waveform through fading channel model
%         chcfg.InitTime = subframeNo/1000;
%         rxWaveform = lteFadingChannel(chcfg, txWaveform);
% 
%         % Add noise
%         noise = N*complex(randn(size(rxWaveform)), ...
%             randn(size(rxWaveform)));
%         rxWaveform = rxWaveform + noise;
% 
%         % Synchronization
%         offset = lteULFrameOffset(frc, frc.PUSCH,rxWaveform);
%         if (offset < 25)
%             offsetUsed = offset;
%         end
% 
%         % SC-FDMA demodulation
%         rxSubframe = lteSCFDMADemodulate(frc, ...
%             rxWaveform(1+offsetUsed:end, :));
% 
%         % Channel Estimation
%         [estChannelGrid, noiseEst] = lteULChannelEstimate( ...
%             frc, frc.PUSCH, cec, rxSubframe);
% 
%         % PUSCH indices for given subframe
%         puschIndices = ltePUSCHIndices(frc,frc.PUSCH);
% 
%         % Minimum Mean Squared Error (MMSE) equalization
%         rxSymbols = lteEqualizeMMSE(rxSubframe(puschIndices), ...
%                         estChannelGrid(puschIndices), noiseEst);
% 
%         % Obtain UL-SCH coding information for current transport block and
%         % ACK, and concatenate this information with the PUSCH / UL-SCH
%         % configuration
%         frc.PUSCH = lteULSCHInfo(frc, frc.PUSCH, trblklen, ...
%             0, 0, 1, 'chsconcat');
% 
%         % Perform deprecoding, demodulation and descrambling on the
%         % received data
%         rxEncodedBits = ltePUSCHDecode(frc, frc.PUSCH, rxSymbols);
% 
%         % UL-SCH channel deinterleaving
%         [deinterleavedData, ccqi, cri, cack] = ...
%             lteULSCHDeinterleave(frc, frc.PUSCH, rxEncodedBits);
% 
%         % HARQ-ACK decoding
%         rxACK = lteACKDecode(frc.PUSCH, cack);
% 
%         % Detect false or missed HARQ-ACK
%         if (isempty(ACK) && ~isempty(rxACK))
%             falseCount = falseCount + 1;
%         end
%         if (~isempty(ACK) && ~isequal(ACK,rxACK))
%             missCount = missCount + 1;
%         end
% 
%     end
% 
%     % Calculate false or missed HARQ-ACK probability
%     pFalse(nSNR) = falseCount/falseTests;
%     pMissed(nSNR) = missCount/missTests;
% 
% end
% 
% hHARQACKResults(SNRdB, pFalse, pMissed);

