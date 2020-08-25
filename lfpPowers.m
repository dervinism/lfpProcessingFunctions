function [lfpPower, f] = lfpPowers(fileName, chN, sr, options)
% [lfpPower, f] = lfpPowers(fileName, chN, sr, options)
%
% Function estimates LFP signal power of various frequency bands,
% calculates the spectrogram, and detects LFP signal trace saturations
% (flat signal) in a given LFP binary file and recording channels of
% interest.
%
% Input: fileName - an LFP binary file name including the full path. In
%          some instances it may be a path to a file containing common
%          reference average (CAR). Or it may also be a cell array holding
%          paths to both. For more details see options.lfpCAR.
%        chN - number of channels in the LFP recording.
%        sr - sampling rate in Hz.
%        options - a structure variable with the following fields:
%          bandRange is a cell array with cells corresponding to different
%            band range frequencies. Below is an example of how to set up
%            one (also the default range, if the field is not set):
%
%                                 delta  alpha  h theta  theta   spindles   slow    beta    s gamma  f gamma   ripples/uf'};
%            options.bandRange = {[1 4]; [4 8];  [5 9];  [8 12]; [6.5 16]; [1 16]; [12 30]; [30 50]; [50 120];  [120 200]};
%
%          chunkSize is a scalar determining the number of recording data
%            samples to load in one instance in order not to overload the
%            computer memory. The default is 4500000.
%          srInterpInit is the initial interpolated sampling rate. All
%            analyses are carried out on data that is down-sampled to this
%            rate. The default is 1000 Hz.
%          srInterpFinal is the final interpolation rate. The output data
%            is further downsampled to this rate. The default is 10 Hz.
%          chOI is a vector with indices of channels of interest. The
%            default is 1.
%          deleteChans is a vector with indices of channels to be removed
%            when carrying out the common average referencing (CAR). The
%            default is none.
%          lfpCAR is a structure used to determine how the common average
%            reference should be used by the function. The following
%            options are available:
%              'none' - CAR is not used in any way (default).
%              'replace' - CAR is used instead of the LFP signal. In this
%                case the user has to supply the full path to the file
%                containing CAR as the first input to the function
%                (fileName).
%              'subtract' - subtracts CAR from the LFP signal.
%              'add' - add CAR to the LFP signal. In this case the user
%                must supply full path to both the lfp recording binary
%                file and the file containing CAR. The paths should be
%                provided via the first input variable in a form of a cell
%                array with the first cell value being the path to the LFP
%                binary file and the second cell being the path to the CAR
%                file.
%          transformFunc specifies the coefficients of a linear
%            transformation function of the LFP signal so that
%
%            tranformed LFP times = transformFunc.a + transformFunc.b * original LFP times
%
%            The default is no transformation. That is, a = 0 and b = 1.
%          powerCalcMethod specifies the method to calculate frequency band
%            power: 'wavelet' (based on wavelet transform of the raw
%            signal) or 'filter' (based on band-pass filtered signal). The
%            default is 'wavelet'.
%          saturationMethod is a structure used to determine the method for
%            finding LFP saturations. Options are the following:
%             'diff' - LFP zero rate of change method.
%             'hist1' - LFP extreme histogram values method. It looks for
%               three saturation values including one around zero.
%             'hist2' - LFP extreme histogram values method. It looks for
%               two saturation values not including one around zero. This
%               is the default method.
%          spectrogram should be set to true if in addition to frequency
%            band power measures you also want to obtain a spectrogram. The
%            default is false.
%          rippleDuration is a minimal duration of a single sharp
%            wave-ripples bout. This is the duration of the impulse that is
%            convolved with a Gaussian. The default is 0.05 seconds. For
%            more details see McGinley et al. (2015).
%          wGaussian is a Gaussian kernel width convolved with a ripple
%            bout in order to obtain the ripple rate. The default is 6 SD
%            (+/-3 SD). For more details consult McGinley et al. (2015).
%          sdGaussian is the duration of a single SD (see McGinley et al.,
%            2015). The default is 1.5 seconds.
%          saturationPlot should be true if you want to display the LFP
%            saturation detection graphs. The default option is false.
%
% Output: lfpPower - a structure variable with the following fields:
%           rippleRate is a cell array with different cells corresponding
%             to ripple rate vectors of different LFP channels oh interest.
%           meanRippleRate is a cell array of average ripple rates of
%             different LFP channels of interest in Hz.
%           theta2deltaRatio (a cell array of vectors).
%           slowPower (a cell array of vectors).
%           fastPower (a cell array of vectors).
%           ultraFastPower (a cell array of vectors).
%           LFPsaturations (a cell array of vectors).
%           nLFPsaturations is a total number of LFP saturations (a cell
%             array of scalars).
%           fLFPsaturations is a number of LFP saturation per minute on
%             average (a cell array fo scalars).
%           meanDurationLFPsaturations is the mean duration of LFP
%             saturations (a cell array fo scalars).
%           wtSpectrogram is a cell array of matrices with each matrix
%             corresponding to a channel of interest. The matrix dimensions
%             match the size of other output vectors on one side and the
%             number of spectral frequencies (fSpectrogram) on the other.
%           fSpectrogram is a vector of spectrogram frequencies.
%           time is the time vector corresponding to output variables.
%           options is the structure variable with options values that were
%             used by the function (same as input options).
%         f - a figure handle of the LFP rate of change graph.
%
% In order to visualise the spectrogram, adapt the following code example:
%   helperCWTTimeFreqPlot(lfpPower.wtSpectrogram{1}, lfpPower.time,...
%     lfpPower.fSpectrogram, 'surf', 'Spectrogram for LFP channel #1',...
%     'Seconds', 'Hz');
%   set(gca, 'YScale', 'log')
%
% References: McGinley MJ, David SV, McCormick DA. Cortical Membrane
%             Potential Signature of Optimal States for Sensory Signal
%             Detection. Neuron. 2015;87(1):179-192.
%             doi:10.1016/j.neuron.2015.05.038


%% Test the input variables
if iscell(fileName)
  medianFile = fileName{2};
  fileName = fileName{1};
end

% Default options
options.bandNames =   {'delta'; 'alpha'; 'h theta'; 'theta'; 'spindles'; 'slow'; 'beta'; 's gamma'; 'f gamma'; 'ripples/uf'};
if ~isfield(options, 'bandRange')
  options.bandRange = { [1 4];   [4 8];    [5 9];    [8 12];  [6.5 16];  [1 16]; [12 30]; [30 50];   [50 120];   [120 200]}; % Hz
end

if ~isfield(options, 'chunkSize')
  options.chunkSize = 4500000; % number of samples to read at a time
end

if ~isfield(options, 'srInterpInit')
  options.srInterpInit = 1000; % Hz
end
if ~isfield(options, 'srInterpFinal')
  options.srInterpFinal = 10; % Hz
end

if ~isfield(options, 'chOI')
  options.chOI = 1;
end
if ~isfield(options, 'deleteChans')
  options.deleteChans = [];
end
if ~isfield(options, 'lfpCAR')
  options.lfpCAR = 'none';
end
if ~isfield(options, 'transformFunc')
  options.transformFunc.a = 0;
  options.transformFunc.b = 1;
end
if ~isfield(options, 'powerCalcMethod')
  options.powerCalcMethod = 'wavelet';
end
if ~isfield(options, 'saturationMethod')
  options.saturationMethod = 'hist2';
end

if ~isfield(options, 'spectrogram')
  options.spectrogram = 0;
end

if ~isfield(options, 'rippleDuration')
  options.rippleDuration = 0.05; % s
end
if ~isfield(options, 'wGaussian')
  options.wGaussian = 6; % Gassian kernel size in SD
end
if ~isfield(options, 'sdGaussian')
  options.sdGaussian = 1.5; % SD size in seconds
end

if ~isfield(options, 'saturationPlot')
  options.saturationPlot = false;
end

bandNames = options.bandNames;
bandRange = options.bandRange;
chunkSize = options.chunkSize;
srInterpInit = options.srInterpInit;
srInterpFinal = options.srInterpFinal;
chOI = options.chOI;
deleteChans = options.deleteChans;
deleteChans = unique(deleteChans);
lfpCAR = options.lfpCAR;
transformFunc = options.transformFunc;
powerCalcMethod = options.powerCalcMethod;
saturationMethod = options.saturationMethod;
spectrogram = options.spectrogram;
rippleDuration = round(options.rippleDuration*srInterpInit);
wGaussian = options.wGaussian;
sdGaussian = options.sdGaussian;
saturationPlot = options.saturationPlot;


%% intialise storage variables
LFPbandPower = {};
wtSpectrogram = {};
LFPsaturations = {};
nLFPsaturations = {};
fLFPsaturations = {};
meanDurationLFPsaturations = {};

rippleRate = {};
meanRippleRate = {};
theta2deltaRatio = {};
slowPower = {};
fastPower = {};
ultraFastPower = {};
f = {};


%% Load and interpolate the data; then apply wavelet transform to estimate LFP band power
% Load
if ~strcmp(lfpCAR, 'replace')
  fid = fopen(fileName, 'r');
  d = dir(fileName);
  nSampsTotal = d.bytes/chN/2;
else
  load(medianFile); %#ok<*LOAD>
  nSampsTotal = numel(medianTrace);
end
nChunksTotal = ceil(nSampsTotal/chunkSize);

chunkInd = 1;
sampleCount = 0;
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  if ~strcmp(lfpCAR, 'replace')
    dat = fread(fid, [chN chunkSize], '*int16');
  else
    if sampleCount+chunkSize > nSampsTotal
      chunkSize = nSampsTotal - sampleCount;
    end
    if chunkSize == 0
      dat = [];
    else
      dat = zeros(chN,chunkSize);
    end
  end
  if ~isempty(dat)
    if strcmp(lfpCAR, 'replace')
      datSamples = size(dat,2);
      if addMedian && (exist('medianTrace','var') && ~isempty(medianTrace))
        median2Add = int16(repmat(medianTrace(sampleCount+1:sampleCount+datSamples),size(dat,1),1));
      else
        median2Add = int16(repmat(zeros(size(sampleCount+1:sampleCount+datSamples)),size(dat,1),1));
      end
      if addMedian && (exist('channelMedian','var') && ~isempty(channelMedian)) %#ok<*USENS,*NODEF>
        median2Add = median2Add + repmat(channelMedian(:,chunkInd),1,size(median2Add,2));
      end
      if numel(swapOrder_temp) < nChans
        median2Add = [median2Add; int16(zeros(1,size(median2Add,2)))];
      end
      dat = dat+median2Add;
      sampleCount = sampleCount + datSamples;
    elseif strcmp(lfpCAR, 'subtract')
      if deleteChans
        chans2include = ones(1,size(dat,1));
        chans2include(deleteChans) = zeros(1,numel(deleteChans));
        chm = zeros(size(dat,1),1);
        chm(logical(chans2include)) = median(dat(logical(chans2include),:),2);
        dat = bsxfun(@minus, dat, int16(chm)); % subtract median of each channel
        tm = int16(median(dat(logical(chans2include),:),1));
      else
        chm = median(dat,2);
        dat = bsxfun(@minus, dat, chm); % subtract median of each channel
        tm = median(dat,1);
      end
      dat = bsxfun(@minus, dat, tm); % subtract median of each time point
    elseif strcmp(lfpCAR, 'add')
      load(medianFile); %#ok<*LOAD>
      datSamples = size(dat,2);
      if (exist('medianTrace','var') && ~isempty(medianTrace))
        median2Add = int16(repmat(medianTrace(sampleCount+1:sampleCount+datSamples),size(channelMedian,1),1));
      else
        median2Add = int16(repmat(zeros(size(sampleCount+1:sampleCount+datSamples)),size(channelMedian,1),1));
      end
      if (exist('channelMedian','var') && ~isempty(channelMedian)) %#ok<*USENS,*NODEF>
        median2Add = median2Add + repmat(channelMedian(:,chunkInd),1,size(median2Add,2));
      end
      dat(1:size(median2Add,1),:) = dat(1:size(median2Add,1),:)+median2Add;
      sampleCount = sampleCount + datSamples;
    end
    
    % Interpolate
    originalTimes = 1/sr:1/sr:size(dat,2)/sr;
    interpTimes = 1/srInterpInit:1/srInterpInit:size(dat,2)/sr;
    interpLFP = interp1(originalTimes, double(dat'), interpTimes)';
    
    % Wavelet transform
    for iCh = 1:numel(chOI)
      LFPbandPowerChunk = zeros(numel(bandRange),numel(interpTimes));
      if strcmpi(powerCalcMethod, 'wavelet')
        fb = cwtfilterbank('SignalLength',size(interpLFP,2),'SamplingFrequency',srInterpInit,...
          'FrequencyLimits',[1 200],'WaveletParameters',[3 16],'VoicesPerOctave',20);
        [wt,fr] = cwt(interpLFP(chOI(iCh),:)-mean(interpLFP(chOI(iCh),:)),'FilterBank',fb); % Continuous wavelet transform
        for iBand = 1:numel(bandRange)
          for iLim = 1:2
            [~, LFPbandLimits{iBand}(iLim)] = min(abs(fr-bandRange{iBand}(iLim))); %#ok<*AGROW>
          end
          LFPbandPowerChunk(iBand,:) = sum(abs(wt(LFPbandLimits{iBand}(2):LFPbandLimits{iBand}(1),:)).^2,1);
        end
      elseif strcmpi(powerCalcMethod, 'filter')
        for iBand = 1:numel(bandRange)
          LFPbandPowerChunk(iBand,:) = leaveFrequencies(interpLFP(chOI(iCh),:)-mean(interpLFP(chOI(iCh),:)),...
            srInterpInit, bandRange{iBand}(1), bandRange{iBand}(2)).^2;
        end
      end
      if chunkInd == 1
        LFPbandPower{iCh} = LFPbandPowerChunk;
      else
        LFPbandPower{iCh} = [LFPbandPower{iCh} LFPbandPowerChunk];
      end
      
      % Spectrogram
      if spectrogram
        fb = cwtfilterbank('SignalLength',size(interpLFP,2),'SamplingFrequency',srInterpInit,...
          'FrequencyLimits',[1 200],'WaveletParameters',[3 16],'VoicesPerOctave',10);
        [wtSpectrogramChunk, fSpectrogram] = cwt(interpLFP(chOI(iCh),:),'FilterBank',fb); % Continuous wavelet transform
        %helperCWTTimeFreqPlot(wt,interpTimes,f,'surf','Spectrogram for CA1 channel','Seconds','Hz');
        %set(gca, 'YScale', 'log')
        if chunkInd == 1
          wtSpectrogram{iCh} = abs(wtSpectrogramChunk).^2;
        else
          wtSpectrogram{iCh} = [wtSpectrogram{iCh} abs(wtSpectrogramChunk).^2];
        end
      end
      
      % Concatenate interpolated LFP traces
      if chunkInd == 1
        LFPsaturations{iCh} = interpLFP(chOI(iCh),:);
      else
        LFPsaturations{iCh} = [LFPsaturations{iCh} interpLFP(chOI(iCh),:)];
      end
    end
  else
    break
  end
  chunkInd = chunkInd+1;
end

% Detect LFP trace saturations
for iCh = 1:numel(chOI)
  [LFPsaturations{iCh}, interpTimes, nLFPsaturations{iCh}, fLFPsaturations{iCh}, meanDurationLFPsaturations{iCh},...
    f{iCh}] = detectLFPsaturations(LFPsaturations{iCh}, interpTimes(2)-interpTimes(1), saturationMethod, saturationPlot, fileName, chOI(iCh));
end


%% Estimate ripple rate
% Ripple rate calculations are based on McGinley et al. (2015)
for iCh = 1:numel(chOI)
  
  % Descriptive measures
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'ripples/uf')
      iRipples = iName;
      break
    end
  end
  %       figure; plot(interpTimes, LFPbandPower{iCh}(iRipples,:)); hold on
  medianPower = median(LFPbandPower{iCh}(iRipples,:),'omitnan');
  stdPower = std(LFPbandPower{iCh}(iRipples,:),'omitnan');
  %       stdPowerTop = medianPower+1.96*stdPower;
  %       stdPowerBottom = medianPower-1.96*stdPower;
  stdPowerTop = medianPower+5*stdPower;
  stdPowerBottom = medianPower-5*stdPower; %#ok<*NASGU>
  %plot(medianPower*ones(1,size(LFPbandPower{iCh},2)))
  %plot(stdPowerTop*ones(1,size(LFPbandPower{iCh},2)))
  %plot(stdPowerBottom*ones(1,size(LFPbandPower{iCh},2)))
  
  % Detect initial ripples
  [ripplePowerPeaks, locations] = findpeaks(LFPbandPower{iCh}(iRipples,:));
  locations = locations(ripplePowerPeaks > stdPowerTop);
  totalRipplesInit = numel(locations);
  %       ripplePowerPeaks = ripplePowerPeaks(ripplePowerPeaks > stdPowerTop);
  %       plot(interpTimes(locations),ripplePowerPeaks, 'r.', 'MarkerSize',5)
  
  % Mark ripple initiation and termination
  ripplePowerPeaksInit = zeros(1,size(LFPbandPower{iCh},2));
  ripplePowerPeaksInit(locations) = 1;
  ripplePowerPeaks = ripplePowerPeaksInit;
  for dtRipple = 1:round(rippleDuration/2)
    ripplePowerPeaks = ripplePowerPeaks + [ripplePowerPeaksInit(1+dtRipple:end) zeros(1,dtRipple)];
    ripplePowerPeaks = ripplePowerPeaks + [zeros(1,dtRipple) ripplePowerPeaksInit(1:end-dtRipple)];
  end
  ripplePowerPeaks(ripplePowerPeaks > 0) = 1;
  locations = 1:size(LFPbandPower{iCh},2);
  locations = locations(logical(ripplePowerPeaks));
  ripplePowerPeaksInit = ripplePowerPeaks;
  %       ripplePowerPeaks = ripplePowerPeaks(ripplePowerPeaks > 0);
  %       plot(interpTimes(locations),ripplePowerPeaks, 'k.', 'MarkerSize',10); hold on
  
  % Total ripple event count
  ripplePowerPeaks2 = findpeaks(ripplePowerPeaksInit);
  totalRipples = numel(ripplePowerPeaks2);
  
  % Convolve with Gaussian
  w = gausswin(wGaussian*sdGaussian*srInterpInit, (wGaussian*sdGaussian*srInterpInit-1)/(2*sdGaussian*srInterpInit));
  w = w/sum(w);
  ripplePowerPeaks = filtfilt(w,1,ripplePowerPeaksInit);
  ripplePowerPeaks = ripplePowerPeaks/sum(ripplePowerPeaks)*sum(LFPbandPower{iCh}(iRipples,:));
  %plot(ripplePowerPeaks)
  %rippleRateInit = ripplePowerPeaks/mean(ripplePowerPeaks)*((totalRipplesInit/size(LFPbandPower{iCh},2))*srInterpInit); % Hz
  rippleRate{iCh} = ripplePowerPeaks/mean(ripplePowerPeaks)*((totalRipples/size(LFPbandPower{iCh},2))*srInterpInit); % Hz
  meanRippleRate{iCh} = mean(rippleRate{iCh});
  %interpTimes = 1/srInterpInit:1/srInterpInit:size(LFPbandPower{iCh},2)/srInterpInit;
  %       plot(interpTimes,rippleRate{iCh});
  
  % Down-sample again
  %interpTimesFinal = 1/srInterpFinal:1/srInterpFinal:size(LFPbandPower{iCh},2)/srInterpInit;
  %rippleRate{iCh} = interp1(interpTimes, rippleRate{iCh}, interpTimesFinal)';
  %interpTimesFinal = 1/srInterpFinal:1/srInterpFinal:size(LFPbandPower{iCh},2)/srInterpInit;
  %plot(interpTimesFinal,rippleRate{iCh}); hold off
end


%% Calculate theta/delta ratio, slow, fast, and ultra-fast activity power
% Again based on McGinley et al. (2015)
for iCh = 1:numel(chOI)
  % Theta/delta ratio
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'delta')
      iDelta = iName;
    end
    if strcmpi(bandNames{iName}, 'h theta')
      iTheta = iName;
    end
  end
  theta2deltaRatio{iCh} = LFPbandPower{iCh}(iTheta,:)./LFPbandPower{iCh}(iDelta,:);
  %d = designFilterLP(1.5, 2, 0.5, 65, srInterpInit);
  %t2dRatioFilt{iCh} = filtfilt(d,theta2deltaRatio{iCh});
  %t2dRatioFilt{iCh}(t2dRatioFilt{iCh} < 1e-9) = 1e-9;
  %figure; plot(interpTimes,theta2deltaRatio{iCh}); hold on
  %plot(interpTimes(2000:end-2000),t2dRatioFilt{iCh}(2000:end-2000));
  %figure; plot(interpTimes,t2dRatioFilt{iCh}); hold on
  %t2dRatioFilt{iCh} = interp1(interpTimes, t2dRatioFilt{iCh}, interpTimesFinal)';
  %plot(interpTimesFinal,t2dRatioFilt{iCh}); hold off
  
  % Slow activity power
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'delta')
      iDelta = iName;
    end
    if strcmpi(bandNames{iName}, 'alpha')
      iAlpha = iName;
    end
    if strcmpi(bandNames{iName}, 'theta')
      iTheta = iName;
    end
  end
  slowPower{iCh} = LFPbandPower{iCh}(iDelta,:) + LFPbandPower{iCh}(iAlpha,:) + LFPbandPower{iCh}(iTheta,:);
  %slowPowerFilt{iCh} = filtfilt(d,slowPower{iCh});
  %slowPowerFilt{iCh}(slowPowerFilt{iCh} < 1e-9) = 1e-9;
  %figure; plot(interpTimes,slowPower{iCh}); hold on
  %plot(interpTimes(2000:end-2000),slowPowerFilt{iCh}(2000:end-2000));
  %figure; plot(interpTimes,slowPowerFilt{iCh}); hold on
  %slowPowerFilt{iCh} = interp1(interpTimes, slowPowerFilt{iCh}, interpTimesFinal)';
  %plot(interpTimesFinal,slowPowerFilt{iCh}); hold off
  
  % Fast activity power
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'beta')
      iBeta = iName;
    end
    if strcmpi(bandNames{iName}, 's gamma')
      iSlowGamma = iName;
    end
    if strcmpi(bandNames{iName}, 'f Gamma')
      iFastGamma = iName;
    end
  end
  fastPower{iCh} = LFPbandPower{iCh}(iBeta,:) + LFPbandPower{iCh}(iSlowGamma,:) + LFPbandPower{iCh}(iFastGamma,:);
  %fastPowerFilt{iCh} = filtfilt(d,fastPower{iCh});
  %fastPowerFilt{iCh}(fastPowerFilt{iCh} < 1e-9) = 1e-9;
  %figure; plot(interpTimes,fastPower{iCh}); hold on
  %plot(interpTimes(2000:end-2000),fastPowerFilt{iCh}(2000:end-2000));
  %figure; plot(interpTimes,fastPowerFilt{iCh}); hold on
  %fastPowerFilt{iCh} = interp1(interpTimes, fastPowerFilt{iCh}, interpTimesFinal)';
  %plot(interpTimesFinal,fastPowerFilt{iCh}); hold off
  
  % Ultra-fast activity power
  for iName = 1:numel(bandNames)
    if strcmpi(bandNames{iName}, 'ripples/uf')
      iRipples = iName;
    end
  end
  ultraFastPower{iCh} = LFPbandPower{iCh}(iRipples,:);
  %ultraFastPowerFilt{iCh} = filtfilt(d,ultraFastPower{iCh});
  %ultraFastPowerFilt{iCh}(ultraFastPowerFilt{iCh} < 1e-9) = 1e-9;
  %figure; plot(interpTimes,ultraFastPower{iCh}); hold on
  %plot(interpTimes(2000:end-2000),ultraFastPowerFilt{iCh}(2000:end-2000));
  %figure; plot(interpTimes,ultraFastPowerFilt{iCh}); hold on
  %ultraFastPowerFilt{iCh} = interp1(interpTimes, ultraFastPowerFilt{iCh}, interpTimesFinal)';
  %plot(interpTimesFinal,ultraFastPowerFilt{iCh}); hold off
end


%% Transform the LFP signal
if transformFunc.a || transformFunc.b ~= 1
  for iCh = 1:numel(chOI)
    rippleRate{iCh} = syncLFP(rippleRate{iCh}, 1/srInterpInit, transformFunc);
    theta2deltaRatio{iCh} = syncLFP(theta2deltaRatio{iCh}, 1/srInterpInit, transformFunc);
    slowPowerFilt{iCh} = syncLFP(slowPowerFilt{iCh}, 1/srInterpInit, transformFunc);
    fastPowerFilt{iCh} = syncLFP(fastPowerFilt{iCh}, 1/srInterpInit, transformFunc);
    ultraFastPowerFilt{iCh} = syncLFP(ultraFastPowerFilt{iCh}, 1/srInterpInit, transformFunc);
    if spectrogram
      wtSpectrogram{iCh} = syncLFP(wtSpectrogram{iCh}, 1/srInterpInit, transformFunc);
    end
    LFPsaturations{iCh} = round(syncLFP(LFPsaturations{iCh}, 1/srInterpInit, transformFunc));
  end
end


%% Filter and further down-sample LFP frequency band power measures
interpTimes = 1/srInterpInit:1/srInterpInit:numel(rippleRate{iCh})/srInterpInit;
interpTimesFinal = 1/srInterpFinal:1/srInterpFinal:interpTimes(end);
for iCh = 1:numel(chOI)
  d = designFilterLP(1.5, 2, 0.5, 65, srInterpInit);
  
  rippleRate{iCh} = interp1(interpTimes, rippleRate{iCh}, interpTimesFinal);
  %plot(interpTimesFinal,rippleRate{iCh}); hold off
  
  t2dRatioFilt = filtfilt(d, theta2deltaRatio{iCh});
  t2dRatioFilt = interp1(interpTimes, t2dRatioFilt, interpTimesFinal);
  t2dRatioFilt(t2dRatioFilt < 1e-9) = 1e-9;
  
  slowPowerFilt = filtfilt(d, slowPower{iCh});
  slowPowerFilt = interp1(interpTimes, slowPowerFilt, interpTimesFinal);
  slowPowerFilt(slowPowerFilt < 1e-9) = 1e-9;
  
  fastPowerFilt = filtfilt(d, fastPower{iCh});
  fastPowerFilt = interp1(interpTimes, fastPowerFilt, interpTimesFinal);
  fastPowerFilt(fastPowerFilt < 1e-9) = 1e-9;
  
  ultraFastPowerFilt = filtfilt(d, ultraFastPower{iCh});
  ultraFastPowerFilt = interp1(interpTimes, ultraFastPowerFilt, interpTimesFinal);
  ultraFastPowerFilt(ultraFastPowerFilt < 1e-9) = 1e-9;
  
  if spectrogram
    wtSpectrogram{iCh} = interp1(interpTimes, wtSpectrogram{iCh}', interpTimesFinal)';
  end
  
  LFPsaturations{iCh} = interp1(interpTimes, LFPsaturations{iCh}, interpTimesFinal);
  
  lfpPower.rippleRate{iCh} = single(rippleRate{iCh});
  lfpPower.meanRippleRate{iCh} = meanRippleRate{iCh};
  lfpPower.theta2deltaRatio{iCh} = single(t2dRatioFilt);
  lfpPower.slowPower{iCh} = single(slowPowerFilt);
  lfpPower.fastPower{iCh} = single(fastPowerFilt);
  lfpPower.ultraFastPower{iCh} = single(ultraFastPowerFilt);
  lfpPower.LFPsaturations{iCh} = single(LFPsaturations{iCh});
end


%% Assign remaining output variables
lfpPower.nLFPsaturations = nLFPsaturations;
lfpPower.fLFPsaturations = fLFPsaturations;
lfpPower.meanDurationLFPsaturations = meanDurationLFPsaturations;
if spectrogram
  lfpPower.wtSpectrogram = wtSpectrogram;
  lfpPower.fSpectrogram = fSpectrogram;
  %helperCWTTimeFreqPlot(wtSpectrogram{1}, timeFinal, fSpectrogram, 'surf', 'Spectrogram for CA1 channel', 'Seconds', 'Hz');
  %set(gca, 'YScale', 'log')
end
lfpPower.time = single(interpTimesFinal);
lfpPower.options = options;