function [LFPsaturations, time, nSaturations, fSaturations, meanSatDuration, f] = detectLFPsaturations(lfp, dt, method, methodPlot, prefix, ch)
% [LFPsaturations, time, nSaturations, fSaturations, meanSatDuration, f] = detectLFPsaturations(lfp, dt, method, methodPlot, prefix, ch)
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
%         methodPlot - if true, draws a graph showing the rate of change of
%           the signal(method = 'diff'). If method = 'hist1' or 'hist2',
%           will draw the LFP values histogram and the LFP trace. The
%           default option is false.
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
if nargin < 6
  ch = 0;
end
if nargin < 5
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
if strcmp(method, 'hist1') || strcmp(method, 'hist2') % histogam method
  nBins = 1000;
  SDfraction = 0.05;
  [lfpHisto, voltage] = hist(lfp,nBins); %#ok<*HIST>
  [count1, peak1Loc] = max(lfpHisto(1:floor(nBins/10)));
  peak1 = voltage(peak1Loc);
  [count2, peak2Loc] = max(lfpHisto(round(nBins/2)-10+1:round(nBins/2)+10));
  peak2Loc = round(nBins/2) - 10 + peak2Loc;
  peak2 = voltage(peak2Loc);
  [count3, peak3Loc] = max(lfpHisto(end-floor(nBins/10)+1:end));
  peak3Loc = nBins - floor(nBins/10) + peak3Loc;
  peak3 = voltage(peak3Loc);
  SD = std(lfp);
  SDfraction = SDfraction*SD;
  minDuration1 = 0.002;
  minDuration2 = 0.1;
  
  [LFPsaturations1, nSaturations1, fSaturations1] = detectionAlgorithm(lfp, time, peak1, SDfraction, minDuration1);
  if strcmp(method, 'hist1')
    [LFPsaturations2, nSaturations2, fSaturations2] = detectionAlgorithm(lfp, time, peak2, SDfraction, minDuration2);
  else
    LFPsaturations2 = zeros(size(lfp)); nSaturations2 = 0; fSaturations2 = 0;
  end
  [LFPsaturations3, nSaturations3, fSaturations3] = detectionAlgorithm(lfp, time, peak3, SDfraction, minDuration1);
  LFPsaturations = zeros(size(lfp));
  LFPsaturations(LFPsaturations1 | LFPsaturations2 | LFPsaturations3) = 1;
  nSaturations = nSaturations1 + nSaturations2 + nSaturations3;
  fSaturations = fSaturations1 + fSaturations2 + fSaturations3;
  meanSatDuration = (sum(LFPsaturations)*dt)/nSaturations;
elseif strcmp(method, 'diff') % rate of change method
  diffLFP = diff(lfp);
  diffLFP = [diffLFP diffLFP(end)];
  stdDiff = std(diffLFP);
  thrDiff = 0.25*stdDiff;
  minDuration = 0.1;
  [LFPsaturations, nSaturations, fSaturations, meanSatDuration] = detectionAlgorithm(diffLFP, time, 0, thrDiff, minDuration);
end

%% Draw graphs
if methodPlot
  saturationTimes = time(logical(LFPsaturations));
  if strcmp(method, 'hist1') || strcmp(method, 'hist2')
    f(1) = figure; hold on
    plot(voltage, lfpHisto);
    if strcmp(method, 'hist1')
      plot(voltage([peak1Loc peak2Loc peak3Loc]), [count1 count2 count3], 'r.', 'MarkerSize',10);
    elseif strcmp(method, 'hist2')
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
    p1 = plot([0 time(end)], [peak1 peak1], 'r');
    p2 = plot([0 time(end)], [peak1-SDfraction peak1-SDfraction], 'g');
    plot([0 time(end)], [peak1+SDfraction peak1+SDfraction], 'g');
    if strcmp(method, 'hist1')
      plot([0 time(end)], [peak2 peak2], 'r');
      plot([0 time(end)], [peak2-SDfraction peak2-SDfraction], 'g');
      plot([0 time(end)], [peak2+SDfraction peak2+SDfraction], 'g');
    end
    plot([0 time(end)], [peak3 peak3], 'r');
    plot([0 time(end)], [peak3-SDfraction peak3-SDfraction], 'g');
    plot([0 time(end)], [peak3+SDfraction peak3+SDfraction], 'g');
    p3 = plot(saturationTimes, zeros(size(saturationTimes)), 'r.', 'MarkerSize',10);
    hold off
    legend([p1 p2 p3], {'saturation values','+/-0.05*SD','LFP saturations'})
    xlabel('Time (s)')
    ylabel('LFP (\muV)')
    title('LFP saturations');
    figName = [prefix '_LFP_saturations_ch' num2str(ch)];
    set(f, 'Name',figName);
  elseif strcmp(method, 'diff')
    f = figure; plot(time,diffLFP); hold on
    p1 = plot([time(1) time(end)], [0 0], 'r');
    p2 = plot([time(1) time(end)], [thrDiff thrDiff], 'g');
    plot([time(1) time(end)], [-thrDiff -thrDiff], 'g')
    p3 = plot(saturationTimes, zeros(size(saturationTimes)), 'r.', 'MarkerSize',10);
    hold off
    legend([p1 p2 p3], {'0 rate of change','0.25*SD','LFP saturations'})
    xlabel('Time (s)')
    ylabel('Rate of change (units/sec)')
    title('LFP rate of change');
    figName = [prefix '_LFP_rate_of_change_ch' num2str(ch)];
    set(f, 'Name',figName);
  end
else
  f = [];
end



%% Local functions
function [LFPsaturations, nSaturations, fSaturations, meanSatDuration] = detectionAlgorithm(signal, time, saturationValue, thresholdDeviation, minDuration)

absDiff = abs(signal-saturationValue);
absDiff(absDiff < thresholdDeviation) = 0;
absDiff(absDiff > 0) = 1;
cumsumDiff = cumsum(absDiff);
saturations = unique(cumsumDiff);
histDiff = hist(cumsumDiff, saturations);
dt = time(2) - time(1);
threshold = minDuration/dt;
histDiff(histDiff < threshold) = 0;
saturations = saturations(logical(histDiff));
nSaturations = numel(saturations);
fSaturations = nSaturations/((time(end) - time(1))/60);
LFPsaturations = zeros(size(signal));
for iSat = 1:numel(saturations)
  LFPsaturations(cumsumDiff == saturations(iSat)) = 1;
end
meanSatDuration = (sum(LFPsaturations)*dt)/nSaturations;