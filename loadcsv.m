function [lfp, time, chIDs] = loadcsv(fileName)
% [lfp, time, chIDs] = loadcsv(fileName)
%
% Function loads a csv file with a header row, a time column, and
% subsequent lfp channel columns.
% Input: fileName - file name string
% Output: lfp - a row matrix of LFP traces corresponding to probe recording
%               channels.
%         time - a time vector.
%         chIDs - channel IDs corresponding to lfp rows.

% Load the probe file
T = readtable(fileName);

% Extract relevant data
time = str2double(table2array(T(2:end,1)))';
chIDs = table2array(T(1,2:end))';
lfp = table2array(T(2:end,2:end))';