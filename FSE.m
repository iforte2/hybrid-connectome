function [rssc] = FSE(SC,time_series,lambda,B)

%ensures time series has no zero time points
[ROIs,len] = size(time_series');
time = time_series;
time( ~any(time,2), : ) = [];

%{
temptime = nonzeros(time);
[a,b] = size(temptime);
timenz = reshape(temptime,[a/ROIs,ROIs]);
%}

%convert to z-scores
ts = zscore(time);
originalData = ts';

%binarize about zero
binarizedData = sign(originalData);

%compute constrained pseudolikelihood
[rssc] = pseudo(binarizedData,SC,lambda,B);

end
