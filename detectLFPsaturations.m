function [LFPsaturations, time, nSaturations, fSaturations, meanSatDuration, f] = detectLFPsaturations(lfp, dt, method, methodPlot, SDfraction, prefix, ch)
% [LFPsaturations, time, nSaturations, fSaturations, meanSatDuration, f] = detectLFPsaturations(lfp, dt, method, methodPlot, SDfraction, prefix, ch)
% Function detects local field potential (LFP) signal saturations (periods
% with flat signal) given a single channel of LFP recording or the median
% value of all LFP channels (common average reference). Three different
% detection algorithms are available to chose from: diff, hist1, and hist2
% (default).
% Inputs: lfp - a single channel LFP signal vector (voltage).
%         dt - signal time step (s).
%         method - a structure to determine the method for finding LFP
%           saturations. Options are the following:
%             'diff' - LFP zero rate of change method.
%             'hist1' - LFP extreme histogram values method. It looks for
%               three saturation values including one around zero.
%             'hist2' - LFP extreme histogram values method. It looks for
%               two saturation values not including one around zero. This
%               is the default method.
%             'combined' - combined hist2 and diff method.
%         methodPlot - if true, draws a graph showing the rate of change of
%           the signal(method = 'diff'). If method = 'hist1' or 'hist2',
%           will draw the LFP values histogram and the LFP trace. The
%           default option is false.
%         SDfraction - fraction of the standard deviation window around
%           saturation voltage value used for saturation detection if hist1
%           or hist2 methods are used (the default value is 0.05 uV) or
%           fraction of the standard deviation window around 0 rate of
%           change value if the diff detection method is used (the default
%           value is 0.25 (uV/s). If combined method is used, one has to
%           specify both values as a two element vector. In this case the
%           default is [0.05 0.25].
%         prefix - figure name prefix (typically recording ID or similar).
%         ch - LFP channel number. it will be converted to figure name
%              suffix.
% Output: LFPsaturations - a vector of ones and zeros marking LFP
%                          saturations by ones. Vector length is the same
%                          as the input LFP vector.
%         time - a corresponding time vector.
%         nSaturations - total number of saturations.
%         nSaturations - saturation frequency per minute.
%         meanSatDuration - average duration of saturations in seconds.
%         f - a figure handle to the graph of the LFP signal rate of
%             change.

%% Initialise variables
if nargin < 7
  ch = 0;
end
if nargin < 6
  prefix = '';
end
if nargin < 4
  methodPlot = false;
end
if nargin < 3
  method = 'hist2';
end

time = dt:dt:dt*numel(lfp);

%% Detect saturations
if strcmp(method, 'hist1') || strcmp(method, 'hist2') || strcmp(method, 'combined') % histogam and combined methods
  nBins = 1000;
  if nargin < 5 || isempty(SDfraction)
    SDfraction1 = 0.05;
  elseif numel(SDfraction) == 2
    SDfraction1 = SDfraction(1);
  else
    SDfraction1 = SDfraction;
  end
  extremeProportions = 6;
  minCount = 50;
  minDuration1 = 0.001;
  minDuration2 = 0.05;
  nPeaks = 3;
  dPeaks = 3;
  [lfpHisto, voltage] = hist(lfp,nBins); %#ok<*HIST>
  if max(abs(voltage)) < 400 || max(abs(voltage)) > 2*400
    maxValue = 0.95*max(abs(voltage));
  else
    maxValue = max([390 min([400 0.8*max(abs(voltage))])]);
  end
  % Peak 1
  [count1, peak1Loc] = findpeaks(lfpHisto(1:floor(nBins/extremeProportions)),...
    'MinPeakHeight',minCount, 'MinPeakProminence',minCount, 'MinPeakDistance',dPeaks);
  pks2include = abs(voltage(peak1Loc)) >= maxValue & count1 > 2*lfpHisto(floor(nBins/extremeProportions));
  count1 = count1(logical(pks2include));
  peak1Loc = peak1Loc(logical(pks2include));
  if ~sum(peak1Loc == 1) && lfpHisto(1) > lfpHisto(2)...
      && lfpHisto(1) >= minCount && abs(voltage(1)) >= maxValue
    count1 = [count1 lfpHisto(1)];
    peak1Loc = [peak1Loc 1];
  end
  [count1, inds] = maxk(count1,nPeaks);
  peak1Loc = peak1Loc(inds);
  if ~isempty(count1)
    peak1 = voltage(peak1Loc);
  else
    peak1 = [];
  end
  % Peak 2
  if strcmp(method, 'hist1')
    [count2, peak2Loc] = max(lfpHisto(round(nBins/2)-10+1:round(nBins/2)+10));
    if count2 < minCount
      count2 = [];
      peak2Loc = [];
      peak2 = [];
    else
      peak2Loc = round(nBins/2) - 10 + peak2Loc;
      peak2 = voltage(peak2Loc);
    end
  end
  % Peak 3
  [count3, peak3Loc] = findpeaks(lfpHisto(end-floor(nBins/extremeProportions)+1:end),...
    'MinPeakHeight',minCount, 'MinPeakProminence',minCount, 'MinPeakDistance',dPeaks);
  pks2include = abs(abs(voltage(end-floor(nBins/extremeProportions)+peak3Loc)) >= maxValue & count3 > 2*lfpHisto(end-floor(nBins/extremeProportions)+1));
  count3 = count3(logical(pks2include));
  peak3Loc = peak3Loc(logical(pks2include));
  if ~sum(peak3Loc == 1) && lfpHisto(end) > lfpHisto(end-1)...
      && lfpHisto(end) >= minCount && abs(voltage(end)) >= maxValue
    count3 = [count3 lfpHisto(end)];
    peak3Loc = [peak3Loc numel(lfpHisto(end-floor(nBins/extremeProportions)+1:end))];
  end
  [count3, inds] = maxk(count3,nPeaks);
  peak3Loc = peak3Loc(inds);
  if ~isempty(count3)
    peak3Loc = nBins - floor(nBins/extremeProportions) + peak3Loc;
    peak3 = voltage(peak3Loc);
  else
    peak3 = [];
  end
  SD = std(lfp);
  SDfraction1 = SDfraction1*SD;
  
  % Peak 1
  LFPsaturations1 = zeros(size(lfp)); nSaturations1 = 0; fSaturations1 = 0;
  if ~isempty(peak1)
    for iPk = 1:numel(peak1)
      [LFPsaturations1temp, nSaturations1temp, fSaturations1temp] = saturationAlgorithm(lfp, time, peak1(iPk), SDfraction1, minDuration1);
      LFPsaturations1 = LFPsaturations1 + LFPsaturations1temp;
      nSaturations1 = nSaturations1 + nSaturations1temp;
      fSaturations1 = fSaturations1 + fSaturations1temp;
    end
  end
  % Peak 2
  if strcmp(method, 'hist1') && ~isempty(peak2)
    [LFPsaturations2, nSaturations2, fSaturations2] = saturationAlgorithm(lfp, time, peak2, SDfraction1, minDuration2);
  else
    LFPsaturations2 = zeros(size(lfp)); nSaturations2 = 0; fSaturations2 = 0;
  end
  % Peak 3
  LFPsaturations3 = zeros(size(lfp)); nSaturations3 = 0; fSaturations3 = 0;
  if ~isempty(peak3)
    for iPk = 1:numel(peak3)
      [LFPsaturations3temp, nSaturations3temp, fSaturations3temp] = saturationAlgorithm(lfp, time, peak3(iPk), SDfraction1, minDuration1);
      LFPsaturations3 = LFPsaturations3 + LFPsaturations3temp;
      nSaturations3 = nSaturations3 + nSaturations3temp;
      fSaturations3 = fSaturations3 + fSaturations3temp;
    end
  end
  % Extreme deviations
  if isempty(peak1)
    LFPsaturations4 = zeros(size(lfp)); nSaturations4 = 0; fSaturations4 = 0;
  else
    LFPsaturations = zeros(size(lfp));
    LFPsaturations(LFPsaturations1 | LFPsaturations2 | LFPsaturations3) = 1;
    [LFPsaturations4, nSaturations4, fSaturations4] = deviationAlgorithm(-lfp, time, max([maxValue max(-peak1)]), LFPsaturations);
  end
  if isempty(peak3)
    LFPsaturations5 = zeros(size(lfp)); nSaturations5 = 0; fSaturations5 = 0;
  else
    LFPsaturations = zeros(size(lfp));
    LFPsaturations(LFPsaturations1 | LFPsaturations2 | LFPsaturations3 | LFPsaturations4) = 1;
    [LFPsaturations5, nSaturations5, fSaturations5] = deviationAlgorithm(lfp, time, max([maxValue max(peak3)]), LFPsaturations);
  end
  if strcmp(method, 'combined')
    if nargin < 5 || isempty(SDfraction)
      SDfraction2 = 0.25;
    else
      SDfraction2 = SDfraction;
    end
    [LFPsaturations6, nSaturations6, fSaturations6, ~, diffLFP, thrDiff] = diffWrap(lfp, time, SDfraction2);
  else
    LFPsaturations6 = zeros(size(lfp)); nSaturations6 = 0; fSaturations6 = 0;
  end
  LFPsaturations = zeros(size(lfp));
  LFPsaturations(LFPsaturations1 | LFPsaturations2 | LFPsaturations3 | LFPsaturations4 | LFPsaturations5 | LFPsaturations6) = 1;
  nSaturations = nSaturations1 + nSaturations2 + nSaturations3 + nSaturations4 + nSaturations5 + nSaturations6;
  fSaturations = fSaturations1 + fSaturations2 + fSaturations3 + fSaturations4 + fSaturations5 + fSaturations6;
  meanSatDuration = (sum(LFPsaturations)*dt)/nSaturations;
elseif strcmp(method, 'diff') % rate of change method
  if nargin < 5 || isempty(SDfraction)
    SDfraction2 = 0.25;
  end
  [LFPsaturations, nSaturations, fSaturations, meanSatDuration, diffLFP, thrDiff] = diffWrap(lfp, time, SDfraction2);
end

%% Draw graphs
if methodPlot
  if sum(LFPsaturations)
    saturationTimes = time(logical(LFPsaturations));
  else
    saturationTimes = [];
  end
  if strcmp(method, 'hist1') || strcmp(method, 'hist2') || strcmp(method, 'combined')
    f(1) = figure; hold on
    plot(voltage, lfpHisto);
    if strcmp(method, 'hist1') && (~isempty(count1) || ~isempty(count2) || ~isempty(count3))
      plot(voltage([peak1Loc peak2Loc peak3Loc]), [count1 count2 count3], 'r.', 'MarkerSize',10);
    elseif (strcmp(method, 'hist2') || strcmp(method, 'combined')) && (~isempty(count1) || ~isempty(count3))
      plot(voltage([peak1Loc peak3Loc]), [count1 count3], 'r.', 'MarkerSize',10);
    end
    hold off
    xlabel('LFP (\muV)')
    ylabel('Count')
    title('LFP histogram');
    figName = [prefix '_LFP_histogram_ch' num2str(ch)];
    set(f, 'Name',figName);
    
    f(2) = figure; hold on
    plot(time,lfp);
    % Peak 1
    if ~isempty(count1)
      for iPk = 1:numel(peak1)
        p1 = plot([0 time(end)], [peak1(iPk) peak1(iPk)], 'r'); %#ok<*NASGU>
        p2 = plot([0 time(end)], [peak1(iPk)-SDfraction1 peak1(iPk)-SDfraction1], 'g');
        plot([0 time(end)], [peak1(iPk)+SDfraction1 peak1(iPk)+SDfraction1], 'g');
      end
    else
      p1 = []; p2 = [];
    end
    % Peak 2
    if strcmp(method, 'hist1') && ~isempty(count2)
      p1 = plot([0 time(end)], [peak2 peak2], 'r');
      p2 = plot([0 time(end)], [peak2-SDfraction1 peak2-SDfraction1], 'g');
      plot([0 time(end)], [peak2+SDfraction1 peak2+SDfraction1], 'g');
    else
      p1 = []; p2 = [];
    end
    % Peak 3
    if ~isempty(count3)
      for iPk = 1:numel(peak3)
        p1 = plot([0 time(end)], [peak3(iPk) peak3(iPk)], 'r');
        p2 = plot([0 time(end)], [peak3(iPk)-SDfraction1 peak3(iPk)-SDfraction1], 'g');
        plot([0 time(end)], [peak3(iPk)+SDfraction1 peak3(iPk)+SDfraction1], 'g');
      end
    else
      p1 = []; p2 = [];
    end
    if sum(LFPsaturations)
      p3 = plot(saturationTimes, zeros(size(saturationTimes)), 'r.', 'MarkerSize',10);
    end
    hold off
    if sum(LFPsaturations)
      if isempty(p1)
        legend(p3, {'LFP saturations'})
      else
        legend([p1 p2 p3], {'saturation values', ['+/-' num2str(SDfraction1/SD) '*SD'], 'LFP saturations'})
      end
    end
    xlabel('Time (s)')
    ylabel('LFP (\muV)')
    title('LFP saturations');
    figName = [prefix '_LFP_saturations_ch' num2str(ch)];
    set(f, 'Name',figName);
    if strcmp(method, 'combined')
      f = diffPlot(diffLFP, time, SDfraction2, thrDiff, saturationTimes, prefix, ch);
    end
  elseif strcmp(method, 'diff')
    f = diffPlot(diffLFP, time, SDfraction2, thrDiff, saturationTimes, prefix, ch);
  end
else
  f = [];
end



%% Local functions
function [LFPsaturations, nSaturations, fSaturations, meanSatDuration] = saturationAlgorithm(signal, time, saturationValue, thresholdDeviation, minDuration)

absDiff = abs(signal-saturationValue);
absDiff(absDiff < thresholdDeviation) = 0;
absDiff(absDiff > 0) = 1;
cumsumDiff = cumsum(absDiff);
saturations = unique(cumsumDiff);
histDiff = hist(cumsumDiff, saturations);
dt = time(2) - time(1);
threshold = round(minDuration/dt) + 1;
histDiff(histDiff < threshold) = 0;
saturations = saturations(logical(histDiff));
nSaturations = numel(saturations);
fSaturations = nSaturations/((time(end) - time(1))/60);
[~,overlapInds] = ismember(cumsumDiff,saturations);
LFPsaturations = zeros(size(signal));
LFPsaturations(logical(overlapInds)) = 1;
meanSatDuration = (sum(LFPsaturations)*dt)/nSaturations;


function [LFPsaturations, nSaturations, fSaturations, meanSatDuration, diffLFP, thrDiff] = diffWrap(lfp, time, SDfraction)

diffLFP = diff(lfp);
diffLFP = [diffLFP diffLFP(end)];
stdDiff = std(diffLFP);
if numel(SDfraction) == 2
  SDfraction2 = SDfraction(2);
else
  SDfraction2 = SDfraction;
end
thrDiff = SDfraction2*stdDiff;
minDuration = 0.1;
[LFPsaturations, nSaturations, fSaturations, meanSatDuration] = saturationAlgorithm(diffLFP, time, 0, thrDiff, minDuration);


function [deviations, nDeviations, fDeviations, meanDevDuration] = deviationAlgorithm(signal, time, thresholdDeviation, exclude)

deviations = zeros(size(signal));
deviations(signal >= thresholdDeviation) = 1;
deviations(logical(exclude)) = 0;
nDeviations = numel(findpeaks(deviations));
fDeviations = nDeviations/((time(end) - time(1))/60);
dt = time(2) - time(1);
meanDevDuration = (sum(deviations)*dt)/nDeviations;


function f = diffPlot(diffLFP, time, SDfraction2, thrDiff, saturationTimes, prefix, ch)

f = figure; plot(time,diffLFP); hold on
p1 = plot([time(1) time(end)], [0 0], 'r');
p2 = plot([time(1) time(end)], [thrDiff thrDiff], 'g');
plot([time(1) time(end)], [-thrDiff -thrDiff], 'g')
p3 = plot(saturationTimes, zeros(size(saturationTimes)), 'r.', 'MarkerSize',10);
hold off
if numel(SDfraction2) == 2
  SDfraction2 = SDfraction2(2);
end
legend([p1 p2 p3], {'0 rate of change', ['+/-' num2str(SDfraction2) '*SD'], 'LFP saturations'})
xlabel('Time (s)')
ylabel('Rate of change (units/sec)')
title('LFP rate of change');
figName = [prefix '_LFP_rate_of_change_ch' num2str(ch)];
set(f, 'Name',figName);