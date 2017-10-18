function [bi,fatb,avsc]=tslburstiness(tsl,varargin)
% ** function [bi,fatb,avsc]=tslburstiness(tsl,varargin)
% computes 'burstiness index' of neuronal firing as proposed by Wagenaar et
% al. (J.Neurosci 25:680-688, 2005). Its value varies between zero (tonic
% firing) and one (extremely bursty firing). Two additional parameters are
% produced; the fraction of time bins with nonzero activity and the
% geometric mean of the spike counts in these bins.
% 
% Input parameters listed below except tsl are optional and must be
% specified as parameter/value pairs, e.g. as in
%          tlsburstiness(tsl,'binW',200);
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT          DESCRIPTION
% tsl            array                 time stamp list of events (ms)
% recTime        array, []             start and stop time of recording (ms)
%                                       if empty, will be set to first and
%                                       last time stamp
% binW           scalar, 100           width of bins
% cutoff         scalar, 0.15          cutoff value: 0.15 means that events
%                                       in upper 15% of bins will be
%                                       counted
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT           DESCRIPTION
% bi             scalar                 burstiness index
% fatb           scalar                 fraction of active time bins
% avsc           scalar                 geometric mean of spike counts in
%                                       active bins


% ----- default values & varargin -----

recTime=[];
binW=100;
cutoff=.15;

pvpmod(varargin);

if isempty(tsl) || any(~isfinite(tsl))
  warning('tsl is empty or contains nonfinite values -  setting bi to nan')
  bi=nan;
  return
end
  

% if recTime was not specified set according to tsl
if isempty(recTime)
  recTime=tsl([1 end]);
end

if cutoff<=0 || cutoff>=1
  error('check cutoff - must be a value between 0 and 1')
end


% do it:
% - total number of events
nEv=numel(tsl);
% - bins
bin=recTime(1):binW:recTime(2);
% - histogram
n=histc(tsl,bin);
% - threshold (number of events corresponding to percentile defined by
%   cutoff)
thresh=prctile(n,(1-cutoff)*100);
% identify bins above thresh and compute fraction of events in these
fEvInBurst=sum(n(n>thresh))/nEv;
% bi is this fraction, minus its expected value in the case of tonic
% firing, divided by the expected fraction of spikes outside said bins
bi=(fEvInBurst-cutoff)/(1-cutoff);
% fraction of active time bins: fraction of bins with more than zero spikes
fatb=sum(n>0)/numel(bin);
% geometric mean of spike counts in active bins
avsc=geomean(n(n>0));