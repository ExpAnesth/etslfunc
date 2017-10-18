function [pepfr,bin]=tslpepfr(refTsl,tsl,varargin)
% ** function [pepfr,bin]=tslpepfr(refTsl,tsl,varargin)
%    computes peri-event PEAK firing rates from time stamp lists: 
%    for each event listed in refTsl the peak rate of events in tsl 
%    within the peri-event window defined below is determined.
%                      ** time unit is ms for all variables **
%
%                    >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT        DESCRIPTION
% refTsl         tsl                 time stamp list containing the reference events 
% tsl            tsl                 time stamp list containing the events to be counted
% interval       array, [-500 1000]  the time interval around the occurrence of a reference event 
%                                     in which peak firing rate in tsl will be counted computed.
%                                      t=0 is the occurrence of the reference event.
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME            TYPE/DEFAULT    DESCRIPTION
% pepfr            2d array       


% defaults
interval=[-500 1000];

pvpmod(varargin);

% the size of things
nTs=size(tsl,1);
nRefTs=size(refTsl,1);

% array collecting data: 1st col peak fr, 2nd col time (within interval
% defined above)
pepfr=nan*zeros(nRefTs,2);

if isempty(tsl),
  warn('tsl is empty - outputting nans');
else   
  totTime=tsl(end,1); 
  % omit events from refTsl which are too close to beginning or end of recording
  badIx=refTsl(:,1)+interval(1)<0.0 | refTsl(:,1)+interval(2)>totTime;
  refTsl(badIx,:)=nan;

  % find entries with interval(2) > tsl-refTsl > interval(1), interval(1) < 0, interval(2) > 0
  for i=1:nRefTs,
    tdiff=tsl(:,1)-refTsl(i,1);  
    tdiff=tdiff(tdiff>=interval(1) & tdiff<=interval(2));  
    if length(tdiff)>=2 
      % def peak fr: inverse of shortest isi (unit: Hz)
      [tmpMin,tmpMinT]=min(diff(tdiff));
      pepfr(i,1)=1000/tmpMin; 
      % def T_peakfr: middle of shortest isi
      pepfr(i,2)=mean(tdiff(tmpMinT+[0 1]));
    end
  end
end
