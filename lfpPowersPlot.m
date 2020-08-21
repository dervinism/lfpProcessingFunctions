function [lfpPower, f, g] = lfpPowersPlot(fileName, chN, sr, options)
% [lfpPower, f, g] = lfpPowersPlot(fileName, chN, sr, options)
%
% Function estimates LFP signal power of various frequency bands,
% calculates the spectrogram, and detects LFP signal trace saturations
% (flat signal) in a given LFP binary file and recording channels of
% interest. It also produces a graph with LFP frequency band measures.
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


%% Estimate LFP frequency band power measures
[lfpPower, g] = lfpPowers(fileName, chN, sr, options);


%% Plot the normalised LFP frequency band power graph
for iCh = 1:numel(options.chOI)
  f{iCh} = figure; %#ok<*AGROW>
  plot(lfpPower.time,lfpPower.rippleRate{iCh}./mean(lfpPower.rippleRate{iCh}));
  hold on
  plot(lfpPower.time,lfpPower.theta2deltaRatio{iCh}./mean(lfpPower.theta2deltaRatio{iCh}))
  plot(lfpPower.time,lfpPower.slowPower{iCh}./mean(lfpPower.slowPower{iCh}))
  plot(lfpPower.time,lfpPower.fastPower{iCh}./mean(lfpPower.fastPower{iCh}))
  plot(lfpPower.time,lfpPower.ultraFastPower{iCh}./mean(lfpPower.ultraFastPower{iCh}))
  plot(lfpPower.time(logical(lfpPower.LFPsaturations{iCh})), zeros(size(lfpPower.time(logical(lfpPower.LFPsaturations{iCh})))), 'r.', 'MarkerSize',10)
  hold off
  legend('Ripple rate','Theta2delta ratio','Slow power','Fast power','Ultra fast power','LFP saturations')
  xlabel('Time (s)')
  ylabel('Normalised signal')
  title(['LFP frequency band power measures: Trace saturation rate of ' num2str(lfpPower.fLFPsaturations{iCh}) ' min^-^1 and mean duration of '...
    num2str(lfpPower.meanDurationLFPsaturations{iCh}) ' s']);
  figName = ['LFP_frequency_band_power_measures_ch' num2str(options.chOI(iCh))];
  set(f{iCh}, 'Name',figName);
end