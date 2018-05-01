clc
clear all
close all

DIP2=[];
DIP3=[];
P=1000;
P1=25;
P2=P1;
N=10;
SINR=[];
enb_serv=[0,0];
enb_2=[5,0];
enb_3=[0,5];
for i=1:5
    ue=[i,i];
    x=[enb_serv;enb_2;enb_3;ue];
    dist=pdist(x,'euclidean');
    d1=dist(3);
    d2=dist(5);
    d3=dist(6);
    s=P/(d1^4);
    DIP2(i)=P1/(d2^4);
    DIP3(i)=P2/(d3^4);
    SINR(i)=s/(DIP2(i)+DIP3(i)+N);
end
SINR=10*log10(SINR);
DIP2=10*log10(DIP2);
DIP3=10*log10(DIP3);
%Plotting locations:
hold on
grid on
v=[0,0,5];
w=[5,0,0];
scatter(v,w,'filled');
a=[1,2,3,4,5];
b=a;
scatter(a,b,'filled','d');
line(a,b);
text(enb_2(1),enb_2(2),'enb_2')
text(enb_3(1),enb_3(2),'enb_3')
text(enb_serv(1),enb_serv(2),'serving enb');

for i=1:5
    txt=SINR(i);
    text(i,i,num2str(txt));
end

NFrames = 10;   % Number of frames

Noc = -98; % dBm/15kHz average power spectral density
% Set the random number generator seed
rng('default');
for idx = 1:5
    
    % Cell 1 eNodeB configuration according to R.46
    enb1 = struct;%A.1.
    enb1.NDLRB = 50;%A.2.
    enb1.CellRefP = 2;%A.3.
    enb1.NCellID = 1;%A.4.
    enb1.CyclicPrefix = 'Normal';%A.5.
    enb1.CFI = 2;%A.6.
    enb1.PHICHDuration = 'Normal';%A.7.
    enb1.Ng = 'Sixth';%A.8.
    enb1.NFrame = 0;%A.9.
    enb1.TotSubframes = 1;%A.10.
    enb1.Windowing = 0;%A.11.
    enb1.DuplexMode = 'TDD';%A.12.
    enb1.SSC = 4;%A.13
    enb1.TDDConfig = 1;%A.14.
    enb1.NPRSRB = 2;
    enb1.IPRS = 0;
    enb1.PRSPeriod = 'On';
    % enb1.Position = [1000, 1000];

    % PDSCH configuration substructure
    enb1.PDSCH.TxScheme = 'TxDiversity'; % PDSCH transmission mode 2 %B.1.
    enb1.PDSCH.Modulation = {'QPSK'};%B.2.
    enb1.PDSCH.NLayers = 2;%B.3.
    enb1.PDSCH.Rho = -3;%B.4.
    enb1.PDSCH.RNTI = 1;%B.5.
    enb1.PDSCH.RVSeq = [0 1 2 3];%B.6.
    enb1.PDSCH.RV = 0;%B.7.
    enb1.PDSCH.NHARQProcesses = 7;%B.8.
    enb1.PDSCH.NTurboDecIts = 5;%B.9.
    enb1.PDSCH.PRBSet = (0:49)';%B.10.
    % Table A.3.4.2.1-2, TS 36.101
    enb1.PDSCH.TrBlkSizes = [5160 3880 0 0 5160 0 3880 0 0 5160];%B.11.
    % Table A.3.4.2.1-2, TS 36.101
    enb1.PDSCH.CodedTrBlkSizes = [12528 10656 0 0 13200 0 10656 0 0 13200];%B.12.
    enb1.PDSCH.CSIMode = 'PUCCH 1-0';%B.13.
    enb1.PDSCH.PMIMode = 'Wideband';%B.14.
    enb1.PDSCH.CSI = 'On';%B.15.

    % PDSCH OCNG configuration
    enb1.OCNGPDSCHEnable = 'On';%B.16.              % Enable OCNG fill
    enb1.OCNGPDSCHPower = -3;   %B.17.              % OCNG power same as PDSCH Rho
    enb1.OCNGPDSCH.RNTI = 0;    %B.18.              % Virtual UE RNTI
    enb1.OCNGPDSCH.Modulation = 'QPSK';%B.19.       % OCNG symbol modulation
    enb1.OCNGPDSCH.TxScheme = 'TxDiversity'; %B.20. % OCNG transmission mode 2

    % Cell2
    enb2 = enb1;
    enb2.NCellID = 2;
    enb2.OCNGPDSCHEnable = 'Off'; 
    enb2.NPRSRB = 2;
    enb2.IPRS = 0;
    enb2.PRSPeriod = 'On';
    % enb2.Position = [1000, -1000];


    % Cell 3
    enb3 = enb1;
    enb3.NCellID = 3;
    enb3.OCNGPDSCHEnable = 'Off';
    enb3.NPRSRB = 2;
    enb3.IPRS = 0;
    enb3.PRSPeriod = 'On';
    % enb3.Position = [-1000, 1000];

    ofdmInfo = lteOFDMInfo(enb1);

    tx = cell(1,3);

    grid1 = [];
    for nsf = 0:19
        enb1.NSubframe = mod(nsf,10);
        sfgrid = lteDLResourceGrid(enb1);       % Empty subframe
        sfgrid(ltePRSIndices(enb1)) = ltePRS(enb1);       % PRS REs
        sfgrid(ltePSSIndices(enb1)) = ltePSS(enb1);       % PSS REs
        sfgrid(lteSSSIndices(enb1)) = lteSSS(enb1);       % SSS REs
        sfgrid(lteCellRSIndices(enb1)) = lteCellRS(enb1); % Cell RS REs
        grid1 = [grid1 sfgrid]; %#ok<AGROW>
    end
    enb1.NSubframe = 0;
    tx{1} = lteOFDMModulate(enb1, grid1);        % OFDM modulate

    grid2 = [];
    for nsf = 0:19
        enb2.NSubframe = mod(nsf,10);
        sfgrid = lteDLResourceGrid(enb2);       % Empty subframe
        sfgrid(ltePRSIndices(enb2)) = ltePRS(enb2);       % PRS REs
        sfgrid(ltePSSIndices(enb2)) = ltePSS(enb2);       % PSS REs
        sfgrid(lteSSSIndices(enb2)) = lteSSS(enb2);       % SSS REs
        sfgrid(lteCellRSIndices(enb2)) = lteCellRS(enb2); % Cell RS REs
        grid2 = [grid2 sfgrid]; %#ok<AGROW>
    end
    enb2.NSubframe = 0;
    tx{2} = lteOFDMModulate(enb2, grid2);        % OFDM modulate

    grid3 = [];
    for nsf = 0:19
        enb3.NSubframe = mod(nsf,10);
        sfgrid = lteDLResourceGrid(enb3);       % Empty subframe
        sfgrid(ltePRSIndices(enb3)) = ltePRS(enb3);       % PRS REs
        sfgrid(ltePSSIndices(enb3)) = ltePSS(enb3);       % PSS REs
        sfgrid(lteSSSIndices(enb3)) = lteSSS(enb3);       % SSS REs
        sfgrid(lteCellRSIndices(enb3)) = lteCellRS(enb3); % Cell RS REs
        grid3 = [grid3 sfgrid]; %#ok<AGROW>
    end
    enb3.NSubframe = 0;
    tx{3} = lteOFDMModulate(enb3, grid3);        % OFDM modulate

    speedOfLight = 299792458.0; % Speed of light in m/s

    %sampleDelay = zeros(1, 1);
    %radius = cell(1, 1);

    %[~, radius{1}] = cart2pol(enb1.Position(1), enb1.Position(2));
    %delay = radius{1}/speedOfLight;                  % Delay in seconds
    %sampleDelay(1) = round(delay*ofdmInfo.SamplingRate); % Delay in samples

    % eNodeB1 to UE propagation channel
    channel1 = struct;                    % Channel config structure
    channel1.Seed = 20;                   % Random channel seed
    channel1.NRxAnts = 2;                 % 2 receive antennas
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

    channel1.SamplingRate = ofdmInfo.SamplingRate;

    % eNodeB2 (interference) to UE propagation channel
    channel2 = channel1;
    channel2.Seed = 122;                  % Random channel seed

    % eNodeB3 (interference) to UE propagation channel
    channel3 = channel1;
    channel3.Seed = 36;                   % Random channel seed

    %hPositioningPlotPositions({enb1, enb2, enb3});

    numSubframes = 10;  % Number of frames to simulate at each SNR
    % SNRdB = [-1.4]; % SNR points to simulate
    frc.TotSubframes = 1;   % Total number of subframes to generate
    frc.NCellID = 5;       % Cell identity
    frc.RC = 'A4-3';        % FRC number

    % Populate FRC configuration structure with default values for A4-3
    frc = lteRMCUL(frc);
    frc.PUSCH.Modulation = {'QPSK'};
    enb1.PUSCH = frc.PUSCH;
    
    
    frc.PDSCH = enb1.PDSCH;

    disp(enb1.PUSCH);
    % frc.Position = hPositioningPosition(0, 0);
    chcfg.NRxAnts = 8;               % Number of receive antennas
    chcfg.DelayProfile = 'EVA';      % Delay profile
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
    chcfg.SamplingRate = info.SamplingRate;

    
    cec = struct;                        %C.1. % Channel estimation config structure
    cec.PilotAverage = 'UserDefined';    %C.2. % Type of pilot symbol averaging
    cec.FreqWindow = 31;                 %C.3. % Frequency window size
    cec.TimeWindow = 23;                 %C.4. % Time window size
    cec.InterpType = 'Cubic';            %C.5. % 2D interpolation type
    cec.InterpWindow = 'Centered';       %C.6. % Interpolation window type
    cec.InterpWinSize = 1;               %C.7. % Interpolation window size


    % Channel noise setup
    nocLin = 10.^(Noc/10)*(1e-3); % linear  in Watts
    % Take into account FFT (OFDM) scaling
    No = sqrt(nocLin/(2*double(ofdmInfo.Nfft)));

    % Signal and interference amplitude scaling factors calculation. These
    % ensure the SINR and DIP values specified are met
    [K1, K2, K3] = hENBscalingFactors(DIP2(idx), DIP3(idx), Noc, SINR(idx) , enb1, enb2, enb3);



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

end
