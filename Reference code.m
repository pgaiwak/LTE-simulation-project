NFrames = 10;   % Number of frames


SINR = -1.4;    % SINR value in dB
DIP2 = -1.73;   % DIP in dB for cell 2
DIP3 = -8.66;   % DIP in dB for cell 3

Noc = -98; % dBm/15kHz average power spectral density

% Set the random number generator seed
rng('default');

% Cell 1 eNodeB configuration according to R.46
enb1 = struct;
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
enb1.OCNGPDSCH.TxScheme = 'TxDiversity'; % OCNG transmission mode 2

% Cell2
enb2 = enb1;
enb2.NCellID = 1;
enb2.OCNGPDSCHEnable = 'Off';

% Cell 3
enb3 = enb1;
enb3.NCellID = 2;
enb3.OCNGPDSCHEnable = 'Off';

% eNodeB1 to UE propagation channel
channel1 = struct;                    % Channel config structure
channel1.Seed = 20;                   % Random channel seed
channel1.NRxAnts = 2;                 % 2 receive antennas
channel1.DelayProfile ='EVA';         % Delay profile
channel1.DopplerFreq = 70;            % Doppler frequency
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
[K1, K2, K3] = hENBscalingFactors(DIP2, DIP3, Noc, SINR, enb1, enb2, enb3);


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


fprintf('\nSimulating %d frame(s)\n',NFrames);

% Main for loop: for all subframes
for subframeNo = 0:(NFrames*10-1)

    % Reinitialize channel seed for each subframe to increase variability
    channel1.Seed = 1+subframeNo;
    channel2.Seed = 1+subframeNo+(NFrames*10);
    channel3.Seed = 1+subframeNo+2*(NFrames*10);

    % Update subframe number
    enb1.NSubframe = subframeNo;
    enb2.NSubframe = subframeNo;
    enb3.NSubframe = subframeNo;

    duplexInfo = lteDuplexingInfo(enb1);

    if duplexInfo.NSymbolsDL ~= 0 % target only downlink subframes

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

        % Serving cell payload
        dlschTransportBlk = {harqProcesses(harqID).dlschTransportBlk};

        % Create transmit waveform and get the HARQ scheduling ID sequence
        % from 'enbOut' structure output which also contains the waveform
        % configuration and OFDM modulation parameters
        [txWaveform1,~,enbOut] = lteRMCDLTool(enb1,dlschTransportBlk);

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

        % Channel time for the present subframe
        channel1.InitTime = subframeNo/1000;
        channel2.InitTime = channel1.InitTime;
        channel3.InitTime = channel1.InitTime;

        % Filter the transmitted waveforms
        rxWaveform1 = lteFadingChannel(channel1,txWaveform1);
        rxWaveform2 = lteFadingChannel(channel2,txWaveform2);
        rxWaveform3 = lteFadingChannel(channel3,txWaveform3);

        % Noise generation
        noise = No*complex(randn(size(rxWaveform1)),...
                randn(size(rxWaveform1)));

        % Add AWGN to the received time domain waveform and scale for
        % required power
        rxWaveform = K1*rxWaveform1 + K2*rxWaveform2 + ...
                     K3*rxWaveform3 + noise;

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

        % Scale rxWaveform by 1/K1 to avoid numerical issues with
        % channel decoding stages
        rxWaveform = (1/K1)*rxWaveform;

        % Perform OFDM demodulation on the received data to obtain
        % the resource grid
        rxSubframe = lteOFDMDemodulate(enb1,rxWaveform);

        % Channel estimation
        [estChannelGrid,noiseEst] = lteDLChannelEstimate(enb1,cec, ...
                                    rxSubframe);

        % Perform equalization, deprecoding, layer demapping,
        % demodulation and descrambling on the received data using the
        % channel estimate.

        % Get PDSCH indices
        [pdschIndices,~] = ...
            ltePDSCHIndices(enb1, enb1.PDSCH, enb1.PDSCH.PRBSet);

        % Get PDSCH resource elements. Scale the received subframe by
        % the PDSCH power factor Rho. The PDSCH is scaled by this
        % amount, while the cell reference symbols used for channel
        % estimation (used in the PDSCH decoding stage) are not.
        [pdschRx, pdschHest] = lteExtractResources(pdschIndices, ...
            rxSubframe*(10^(-enb1.PDSCH.Rho/20)), estChannelGrid);

        % Decode PDSCH
        [dlschBits,~] = ltePDSCHDecode(enb1, enb1.PDSCH,...
            pdschRx, pdschHest, noiseEst);

        % Decode downlink shared channel (DL-SCH)
        [decbits,harqProcesses(harqID).crc,harqProcesses(harqID).decState] = ...
            lteDLSCHDecode(enb1,enb1.PDSCH,harqProcesses(harqID).trBlkSize, ...
            dlschBits{1},harqProcesses(harqID).decState);

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

% Display achieved throughput percentage
disp(['Achieved throughput ' num2str(simThroughput*100/maxThroughput) '%'])

% Plot running throughput
figure;plot(runningSimThPut*100./runningMaxThPut)
ylabel('Throughput (%)');
xlabel('Simulated subframe');
title('Throughput');

