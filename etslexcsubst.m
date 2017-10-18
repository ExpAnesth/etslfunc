function [d,sIx,eIx]=etslexcsubst(d,si,tsl,sIntv,eIntv,varargin)
% ** function [d,sIx,eIx]=etslexcsubst(d,si,tsl,sIntv,eIntv,varargin)
%    substitutes events as defined by (extended) time stamp lists (e.g.
%    stimulation artifacts, spikes in intracellular data) in time series
%    d by inconspicuous data values (e.g. fixed value or by extending
%    linear fit of surrounding points).
%
%                    >>> INPUT VARIABLES >>>
% NAME       TYPE/DEFAULT           DESCRIPTION
% d          1d                     sampled data, time runs along column
% tsl        1d arr                 time stamp list of artifacts (best used
%                                    if artifacts have constant width)
%            2d arr                 *extended* time stamp list (if arti-
%                                    facts vary in width)
% si         a) scalar              sampling interval in us
%            b) string 'idx'        means that ALL timing information is
%                                    specified in points (as opposed to ms)
% sIntv      2-element arr          peri-event substitution interval (ms 
%                                    or points; see NOTE 1 below)
% eIntv      2-element arr          peri-event estimation interval (ms or 
%                                    points) from which substitute data
%                                    will be estimated; will be ignored if
%                                    sVal is specified as input; also see
%                                    NOTE 1 and NOTE 2 below
% sVal       scalar, []             fixed substitute value: all points in 
%                                    sIntv will be set to this value; eIntv
%                                    will be ignored
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
%                                      and 2, respectively) to estimation
%                                      intervals
% 
% ** NOTE 1: etslexcsubst was overhauled in May 2014 to allow elimination
% of artifacts of varying width. Besides time stamp lists (tsl) it now
% accepts *extended* time stamp lists (etsl) which are designed for
% instances of neuronal activity with a varying duration ('bursts') and
% therefore contain information on the duration (width) of the artifact.
% Here's the important bit: time intervals sIntv and eIntv mean different
% things depending on whether a tsl or an etsl is specified as an input: by
% definition, tsl lack information on the duration of events, so it is
% assumed to be zero. Therefore, input variable sIntv must encompass the
% real width of the events (e.g. [-0.1 1.1] would be appropriate for
% artifacts of 1 ms length, because in this case an extension of the
% elimination window by 0.1 ms on either side of the artifact would make
% reasonably sure it got wiped out completely. Remember that in detected
% events zero is the time of the threshold crossing). If an etsl is
% specified as input, events do have a duration. In this case, sIntv of
% e.g. [-0.3 0.5] would mean that the artifact would be eliminated in a
% window 0.3 ms before its start up to a point 0.5 ms after its end. The
% same logic applies if zero-duration events occur so close to each other
% that they must be combined into bursts.
% ** NOTE 2: eIntv, the interval from which substitute data will be
% estimated, should in theory not overlap with sIntv, because the point of
% function tslev


doSubst=true;
sVal=[];
pvpmod(varargin);

etslconst;

[n1 n2]=size(d);
if n2>1
  error('d must be a single-column array')
end

etslconst;
if size(tsl,2)>1
  % extended tsl as input - retain only ts and dur columns
  etsl=tsl(:,[etslc.tsCol etslc.durCol]);
  isFixedLen=false;
else
  % tsl as input - make into etsl and rename
  etsl=[tsl tsl*nan];
  isFixedLen=true;
end
% etsl has to be sorted
etsl=sortrows(etsl,1);
nTs=size(etsl,1);

% check sIntv
if diff(sIntv)<=0
  error('values in sIntv must be in ascending order');
end

if ~isempty(sVal)
  isSubstConstant=true;
  % although eIntv is not needed in this case, set it to sIntv so the code
  % below will run through without glitch
  eIntv=sIntv;
else
  isSubstConstant=false;  
  sVal=nan;
  % now check eIntv:
  if diff(eIntv)<=0
    error('values in eIntv must be in ascending order');
  end
  % - make sure that at least a portion of eIntv is outside sIntv
  if sIntv(1)<=eIntv(1) && sIntv(2)>=eIntv(2)
    error('eIntv is completely contained in sIntv');
  end
end

if ischar(si)
  if ~strcmpi(si,'idx')
    error('input variable si must be either a scalar or string ''idx''');
  end
else
  % ** convert to ticks (easier to calculate with) **
  sIntv=cont2discrete(sIntv,si*.001,'intv',1)-1;
  eIntv=cont2discrete(eIntv,si*.001,'intv',1)-1;
  etsl=cont2discrete(etsl,si*.001,'intv',0);
end

if nTs
  % kick out start ts outside of borders
  killIx=etsl(:,1)<1 | etsl(:,1)>n1;
  if any(killIx)
    warning('etslexcsubst:timestamp_outside_data',[int2str(numel(find(killIx))),...
      ' events/artifacts as listed in etsl are outside of borders of data d']);
    etsl(killIx,:)=[];
    nTs=size(etsl,1);
  end
end

sIx=zeros(0,2);
eIx=zeros(0,2);

if nTs & isFixedLen
  % in case we're dealing with fixed-length events: detect sequences of
  % overlapping events (=bursts) 
  % - if ts are separated by less than a minimal distance combine them into
  % bursts. A reasonable value is the duration of an artifact as specified
  % in sIntv
  maxIEI=diff(sIntv);
  % turn off the 'considering single events as bursts' warning
  warning('off', 'etslburstf:singleEvBu')
  etsl=etslburstf(etsl(:,1),maxIEI,'minNEvPerBurst',1,'recLen',n1);
  warning('on', 'all')
  % retain only ts and dur columns
  etsl=etsl(:,[etslc.tsCol etslc.durCol]);
  nTs=size(etsl,1);
end

if nTs
  startStopList=cumsum(etsl(:,[etslc.tsCol etslc.durCol]),2);
  % substitution interval
  sIx=repmat(sIntv,nTs,1)+startStopList;
  % extrapolation interval (potentially including points not to be used)
  eIx=repmat(eIntv,nTs,1)+startStopList;
  % because sIntv and eIntv are independent correct border points for them
  % separately
  % - sIntv:
  sIx(sIx(:,1)<1,1)=1;
  sIx(sIx(end,2)>n1,2)=n1;
  % - eIntv: we have to check both borders of eIntv
  ix1=eIx<1;
  ix2=eIx>n1;
  % adjust eIx accordingly 
  eIx(ix1(:,1),1)=1;
  eIx(ix2(:,2),2)=n1;
  % these are the (unlikely) cases in wich the complete eIntv is outside
  % data range, so account for that
  impossibleIx=all(ix1,2) | all(ix2,2);
  if any(impossibleIx)
    warning('etslexcsubst:eIntv_outside_data',[int2str(numel(find(impossibleIx))),...
      ' events/artifacts are too close to data borders to be substituted with the current choice of eIntv']);
    eIx(impossibleIx,:)=nan;
  end
  % now substitute
  if doSubst
    % - by constant value
    if isSubstConstant
      for g=1:nTs
        d(sIx(g,1):sIx(g,2),1)=sVal;
      end
    % - by extending linear fit of pre-interval
    else
      warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
      for g=1:nTs
        if ~impossibleIx(g)
          % indexes in sIx which are to be used for estimation (they are
          % outside of eIx)
          x=setdiff(eIx(g,1):eIx(g,2),sIx(g,1):sIx(g,2))';
          % determine linear fit through estimation interval
          fitPar=polyfit(x,d(x),1);
          % extend fit
          d(sIx(g,1):sIx(g,2),1)=fitPar(2)+(sIx(g,1):sIx(g,2))*fitPar(1);
        end
      end
      warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');
    end
  end
else
  disp('tsl is empty or devoid of suitable events - nothing to clean');
end