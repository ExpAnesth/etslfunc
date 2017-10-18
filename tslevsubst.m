function [d,sIx,eIx]=tslevsubst(d,si,tsl,sIntv,eIntv,varargin)
% ** function [d,sIx,eIx]=tslevsubst(d,si,tsl,sIntv,eIntv,varargin)
%    substitutes events as defined by time stamp lists (e.g. spikes) in
%    continuous raw data d by inconspicuous data values (e.g. mean of
%    surrounding points).
%
%                    >>> INPUT VARIABLES >>>
% NAME       TYPE/DEFAULT           DESCRIPTION
% d          1d                     sampled data, time runs along column
% tsl        1d arr                 time stamp list  
% si         a) scalar              sampling interval in us
%            b) string 'idx'        means that ALL timing information is
%                                    specified in points (as opposed to ms)
% sIntv      2-element arr          peri-event substitution interval (ms 
%                                    or points)
% eIntv      2-element arr          peri-event estimation interval (ms or 
%               OR                   points) from which substitute data 
%                                    will be estimated
%            scalar                 substitute value (all points in sIntv
%                                    will be set to this value
% doSubst    logical, true          if true, events will be substituted
%                                    (set to false if you need only sIx and
%                                    eIx as output) 
%
%                  <<< OUTPUT VARIABLES <<<
% NAME        TYPE/DEFAULT            DESCRIPTION
% d           same as input           cleaned raw data
% sIx         n by 2 array            start and stop indexes (columns 1 
%                                      and 2, respectively) to substitution
%                                      intervals
% eIx         n by 2 array            start and stop indexes (columns 1
%                                       and 2, respectively) to estimation
%                                       intervals

warning('tslevsubst is outdated - use etslexcsubst');

doSubst=true;
pvpmod(varargin);

etslconst;

[n1 n2]=size(d);
if n2>1
  error('d must be a single-column array')
end
if size(tsl,2)>1
  warning('tsl has more than one column, eliminating all but the first');
  tsl=tsl(:,1);
end

% tsl has to be sorted
tsl=sort(tsl);

if diff(sIntv)<=0
  error('check sIntvC');
end

if numel(eIntv)==1
  isSubstConstant=true;
  substVal=eIntv;
  % set eIntv to sIntv so the code below will run through without glitch
  eIntv=sIntv;
else
  isSubstConstant=false;  
  substVal=nan;
end

if ~isSubstConstant
  if sIntv(1)<eIntv(1) || sIntv(2)>eIntv(2)
    error('sIntv must be contained in eIntv');
  end
end
sIx=zeros(0,2);
eIx=zeros(0,2);

if ischar(si)
  if ~strcmpi(si,'idx')
    error('input variable si must be either a scalar or string ''idx''');
  end
else
  % ** convert to ticks (easier to calculate with) **
  sIntv=cont2discrete(sIntv,si*.001,'intv',1)-1;
  eIntv=cont2discrete(eIntv,si*.001,'intv',1)-1;
  tsl=cont2discrete(tsl,si*.001);
end

% kick out ts completely outside of borders
tsl(tsl<1 | tsl>n1)=[];
nTs=numel(tsl);

if nTs
  % ** there may be sequences of overlapping events, in other words,
  % bursts, so detect them and from here on deal with an etls
  
  % if ts are separated by less than this distance combine them into bursts
  % - this way we also ensure that neither of the two extrapolation flanks
  % extend into the intervals to be substituted of neighboring events
  maxIEI=max(abs(eIntv))+1;
  % turn off the 'considering single events as bursts' warning
  warning('off', 'etslburstf:singleEvBu')
  etsl=etslburstf(tsl,maxIEI,'minNEvPerBurst',1,'recLen',n1);
  warning('on', 'all')
  nTs=size(etsl,1);
end

if nTs
  stasto=cumsum(etsl(:,[etslc.tsCol etslc.durCol]),2);
  % substitution interval
  sIx=repmat(sIntv,nTs,1)+stasto;
  % extrapolation interval
  eIx=repmat(eIntv,nTs,1)+stasto;
  % correct border points (by definition, only the first and last event
  % can cause trouble)
  if eIx(1,1)<1
    eIx(1,1)=1;
    if sIx(1,1)<1
      % if sIx is also out of bounds set 2nd data point as the first one to
      % be substituted and substitute the missing pre-event part by the
      % mean of the post-event one as a best guess
      sIx(1,1)=2;
      d(1)=mean(d(sIx(1,2)+1:eIx(1,2)));
    end
  end
  % correspondingly for last event
  if eIx(end,2)>n1
    eIx(end,2)=n1;
    if sIx(end,2)>n1
      sIx(end,2)=n1-1;
      d(end)=mean(d(eIx(end,1):sIx(end,1)));
    end
  end
  % now substitute
  if doSubst
    if isSubstConstant
      for g=1:nTs
        d(sIx(g,1):sIx(g,2),1)=substVal;
      end
    else
      for g=1:nTs
        x=[eIx(g,1):sIx(g,1)-1   1+sIx(g,2):eIx(g,2)]';
        d(sIx(g,1):sIx(g,2),1)=interp1(x,d(x),sIx(g,1):sIx(g,2));
      end
    end
  end
else
  disp('tsl is empty or devoid of suitable events - nothing to clean');
end