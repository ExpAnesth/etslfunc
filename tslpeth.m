function [peth,bin,varargout]=tslpeth(refTsl,tsl,varargin)
% ** function [peth,bin,varargout]=tslpeth(refTsl,tsl,varargin)
%    computes peri-event time histogram from time stamp lists. 
%    Works similarly to tslcc; the differences are:
%      - peth will be computed for each reference event individually
%      - events belonging to next reference event can be excluded
%
%                      ** time unit is ms for all variables **
%
%                    >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT        DESCRIPTION
% refTsl         tsl                 time stamp list containing the reference events 
% tsl            tsl                 time stamp list containing the events to be counted
% interval       array, [-500 1000]  the time interval around the occurrence of a reference event 
%                                     in which events in tsl will be counted (=added to the histogram).
%                                     t=0 is the occurrence of the reference event.
% intvDistance   scalar, nan         if 'interval' is so long that it extends into 
%                                     the interval of neighboring reference events
%                                     this parameter determines the minimal distance 
%                                     events in tsl must have to the next reference event 
%                                     (e.g. a value of 50 says that all events of the current 
%                                     ref event less than 50 ms before the next ref event 
%                                     shall be omitted). Set to nan if no minimal
%                                     distance is desired
% binw           scalar, 2           the bin width in ms
% convHz         scalar, 0           if nonzero output will be converted to
%                                     events per second (Hz)
%                    <<< OUTPUT VARIABLES <<<
%
% NAME            TYPE/DEFAULT    DESCRIPTION
% peth            2d array        the peth, one column per event in refTsl 
% bin             1d array        the BORDERS of the bins of the peth
% varargout{1}    cell array      time stamps in interval, time shifted
%                                 (one cell per reference event)


% defaults
interval=[-500 1000];
intvDistance=nan;
binw=2;
verbose=0;
convHz=0;
% modify according to input
pvpmod(varargin);

if verbose
  disp(['**** ' mfilename ':']);
end

% the size of things
nTs=size(tsl,1);
nRefTs=size(refTsl,1);

% time window 
bin=(interval(1):binw:interval(2))';
% array collecting data
peth=zeros(numel(bin),nRefTs);
% pre-generate the list of possibly dynamic interval limits 
intvEnd=repmat(interval(2),size(refTsl));
if isfinite(intvDistance)
  if nRefTs>1
    intvEnd=[min(refTsl(2:end)-refTsl(1:end-1)-intvDistance,intvEnd(2:end)); inf];
  end
end

% cell array of local tsls
if nargout>2
  doCollectTsl=true;
  tslColl=cell(1,nRefTs);
else
  doCollectTsl=false;
end
  
if isempty(tsl),
  peth=zeros(length(bin),1)*nan;
  warning('tsl is empty - outputting nans');
else   
  refTotTime=refTsl(end,1);
  totTime=tsl(end,1); 
  timeDiff=abs(totTime-refTotTime);
  % omit events from refTsl which are too close to beginning or end of recording
  badIx= refTsl(:,1)+interval(1)<0.0 | refTsl(:,1)+interval(2)>totTime;
  refTsl(badIx,:)=nan;
  peth(:,badIx)=nan;
  % find entries with intvEnd > tsl-refTsl > interval(1), interval(1) < 0, interval(2) > 0
  for i=1:nRefTs,
    tdiff=tsl(:,1)-refTsl(i,1);  
    tdiff=tdiff(tdiff>=interval(1) & tdiff<=intvEnd(i));  
    if doCollectTsl
      % stuff tsl in cell
      tslColl{i}=tdiff;
    end
    % reference events with zero events are included due to preallocation
    if ~isempty(tdiff)
      tmpCount=histc(tdiff,bin);
      peth(:,i)=tmpCount; 
    end
  end
  % convert to Hz (1/(spikes per bin))
  if convHz
    peth=peth/binw*1000;
  end
end

if nargout>2
  varargout{1}=tslColl;
end