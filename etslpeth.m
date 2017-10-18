function [peth,bin,varargout]=etslpeth(refEtsl,tsl,varargin)
% ** function [peth,bin,varargout]=etslpeth(refEtsl,tsl,varargin)
%    computes peri-event time histogram from time stamp lists. 
%            *** time unit is ms for all variables ***
%                    >>> INPUT VARIABLES >>>
% NAME           TYPE/DEFAULT        DESCRIPTION
% refEtsl        etsl                extended time stamp list containing  
%                                     the reference events
% tsl            tsl                 time stamp list containing the events 
%                                     to be counted
% interval       array, [-500 1000]  time interval around the occurrence of
%                                     a reference event in which events in
%                                     tsl will be counted (=added to the
%                                     histogram). t=0 is the start of the
%                                     reference event.
% binw           scalar, 2           bin width in ms (regularly spaced
%                                     bins)
% bin            array               bins (left bin borders)
%                                     ** NOTE: if this input argument is
%                                     specified it overrides both 'interval'
%                                     and 'binw', and output will be
%                                     converted to Hz (see below) because
%                                     bins of unequal width are allowed
% convHz         scalar, 0           if nonzero, output will be converted 
%                                     to events per second (Hz)
% iRestrict      {scalar,char},     interval restriction: 
%                {'none',0}          - if set to e.g. {'cur',100} only
%                                    events occurring up to 100 ms after the
%                                    END of the CURRENT ref event are collected;
%                                    the remaining bins are set to nan.
%                                    - if set to e.g. {'next',-100} only
%                                    events occurring up to 100 ms
%                                    BEFORE the NEXT ref event are
%                                    collected; the remaining bins are set
%                                    to nan
%                                    - set to {'none',x}, x an arbitrary
%                                    scalar, if neither restriction of the 
%                                    interval is intended
%
%                    <<< OUTPUT VARIABLES <<<
% NAME            TYPE/DEFAULT    DESCRIPTION
% peth            2d array        the peth, one column per event in refEtsl 
% bin             1d array        the BORDERS of the bins of the peth
% varargout{1}    1d array        bin widths
% varargout{2}    struct          a struct containing additional parameters
%                                  that may be of interest: 
%                                  - fractionEvInBurst, the fraction of spx
%                                    in bursts as defined in etsl
%                                  - evRateInBurst, average spx rate (Hz) 
%                                    within bursts



% defaults
interval=[-500 1000];
iRestrict={'none',0} ;
binw=2;
bin=[];
convHz=0;
% modify according to input
pvpmod(varargin);

etslconst;
% counting spx in bursts: spx clearly belonging to a burst may precede the
% LFP by a few ms; the value below specifies the offsets in ms by which
% right and left borders of the interval in which spx are to be collected
% shall be extended (t=0 = beginning of burst)
offsIntvSpxCt=[-10 10];

% the size of things
nTs=size(tsl,1);
nRefTs=size(refEtsl,1);

if ~isempty(bin)
  % ** user-defined bins 
  [tmpr,tmpc]=size(bin);
  if ~any([tmpr tmpc]==1)
    error('input arg ''bin'' must be a single-column (or single-row) array');
  end
  % make sure values are ascending
  if any(diff(bin)<=0)
    error('values in input arg ''bin'' must be ascending');
  end
  % column array
  bin=bin(:);
  % extend by width of last bin (which may result in biased counts in case
  % of heterogeneous bin widths) 
  bin=[bin; bin(end)+diff(bin(end-1:end))];
  % overwrite interval
  interval=bin([1 end]); 
  convHz=1;
else
  % bins (more precisely, bin borders)  ** note the extension of the
  % interval to the right because of histc
  bin=(interval(1):binw:interval(2)+binw)';
end

% 'compute' bin width - the last bin will be eliminated further below so
% that this array of bin widths matches the array of bin borders in terms
% of number of elements
binw=diff(bin);
if nargout>2
  varargout{1}=binw;
end
if nargout>3
  pethPar.fractionEvInBurst=nan;
  pethPar.evRateInBurst=nan;
  varargout{2}=pethPar;
end
% array collecting data
peth=zeros(numel(bin),nRefTs);
% number of spx residing in bursts as defined in etsl
spxInBurst=0;

% pre-generate the list of possibly dynamic interval limits 
if ~iscell(iRestrict)
  error('input argument ''iRestrict'' must be a cell array');
else
  % - the default case, no restriction 
  intvEnd=repmat(interval(2),nRefTs,1)+binw(end);
  switch iRestrict{1}
    case 'cur'
      intvEnd=refEtsl(:,etslc.durCol)+iRestrict{2};
    case 'next'
      if nRefTs>1
        intvEnd=[min(refEtsl(2:end,etslc.tsCol)-refEtsl(1:end-1,etslc.tsCol)+iRestrict{2},intvEnd(2:end)); inf];
      end
    case 'none'
    otherwise
  end
end

if isempty(tsl)
  % **note the -1 
  peth=zeros(length(bin)-1,1)*nan;
  warning('tsl is empty - outputting nans');
else   
  % time stamps representing the ends of bursts
  eob=sum(refEtsl(:,[etslc.tsCol etslc.durCol]),2); 
  % amplitude column will be used for counting number of spx per reference
  % event
  refEtsl(:,etslc.amplCol)=nan;
  % find events from refEtsl which are too close to beginning or end of
  % recording
  badIx=refEtsl(:,1)+interval(1)<0.0 | refEtsl(:,1)+intvEnd>eob(end);
  goodIx=find(~badIx);
  % time interval in which spx are counted for computation of fraction of
  % spx in bursts: from beginning of first non-bad burst to last burst (no
  % matter whether it is 'bad' or not)
  if ~isempty(goodIx)
    tmpIntv=refEtsl([goodIx(1) end],etslc.tsCol);
  else
    tmpIntv=[nan nan];
  end
  % omit 'bad' events (in above sense)
  refEtsl(badIx,:)=nan;
  peth(:,badIx)=nan;
  % find entries with intvEnd > tsl-refEtsl > interval(1), interval(1) < 0, interval(2) > 0
  for g=1:nRefTs
    % 'local time' time stamps
    tdiff=tsl(:,1)-refEtsl(g,1);  
    if ~isempty(tdiff)
      % number of spx residing in burst
      refEtsl(g,etslc.amplCol)=numel(find(tdiff>offsIntvSpxCt(1) & tdiff<=refEtsl(g,etslc.durCol)+offsIntvSpxCt(2)));
    end
    % the local time stamps falling into PETH bins 
    tdiff=tdiff(tdiff>=interval(1) & tdiff<=intvEnd(g));  
    % index to 'allowed' bins
    tmpIx=bin>intvEnd(g);
    if ~isempty(tdiff)
      tmpCount=histc(tdiff,bin);
      peth(:,g)=tmpCount; 
    end
    % reference events with zero events are included due to preallocation,
    % but intvEnd should nonetheless be represented via NaNs
    peth(tmpIx,g)=nan;
  end
  % get rid of last bin which contains only border counts anyways
  peth(end,:)=[];
  % convert to Hz (1/(spikes per bin))
  if convHz
    peth=peth./repmat(binw,[1 nRefTs])*1000;
  end
  % now kick last bin
  bin(end)=[];
  if nargout>3
    % if additional analysis parameters were requested provide them
    % - fraction of spx in bursts:
    pethPar.fractionEvInBurst=nansum(refEtsl(:,etslc.amplCol))/...
    numel(find(tsl(:,1)>=tmpIntv(1) & tsl(:,1)<tmpIntv(end)));
    % - spx rate (Hz) within bursts: compute for each burst individually,
    % then average
    pethPar.evRateInBurst=nanmean(refEtsl(:,etslc.amplCol)./...
      refEtsl(:,etslc.durCol))*1000;
    varargout{2}=pethPar;
  end
end