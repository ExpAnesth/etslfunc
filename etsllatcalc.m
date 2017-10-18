function lat=etsllatcalc(intv,varargin)
% ** function lat=etsllatcalc(intv,varargin)
% calculates the latencies between events listed in extended time 
% stamp lists (etsl). This function works similar to a peth-generator, 
% the major difference being that for each reference event only 
% the closest event (in time) in each of the other lists is picked.
%
%                    >>> INPUT VARIABLES >>>
% NAME      TYPE/DEFAULT     DESCRIPTION
% varargin  etsls            the time stamp lists. The first of these is 
%                            the reference for all other lists (termed 
%                            'matching', see below)
% intv      2-element array  the time interval around the occurrence of a 
%                            reference event (t=0) in which events in the 
%                            matching time stamp list(s) will be regarded as 
%                            candidate events (relative to which the latency 
%                            will be calculated). [-100 200] for example means 
%                            that only events at most 100 ms before and 200 ms 
%                            after the occurrence of the reference event will 
%                            be taken into consideration. Among these, the event 
%                            with the minimal distance (latency) will be chosen. 
%                            If for a given reference event no event in the 
%                            matching time stamp list occurs within this interval, 
%                            its latency value is set to NaN.
%
%                    <<< OUTPUT VARIABLES <<<
% NAME       TYPE/DEFAULT    DESCRIPTION
% lat        2D-array        first column: time stamps of the reference tsl
%                            second+ columns: latencies
%                            Negative values of the latency mean that the 
%                            reference event occurred before the matching event.

% improvements:
% - vectorized 
% - reverse sign of latency


etslconst;
disp(['**** ' mfilename ':']);
numOfLists=size(varargin,2);
if numOfLists<2, error('more than one time stamp list needed'); end
if intv(1)>0.0 | intv(2)<0.0, error('check intv'); end

refetsl=varargin{1};
refEvIdx=1:length(refetsl);
tmp=refEvIdx(end);
if size(refetsl,2)<etslc.tagCol
  disp('reference tsl lacks column for tags, no check for rejected events possible');
else
  disp('scanning reference list for rejected events..');
  gutso=find(refetsl(:,etslc.tagCol)~=etslc.reject);
  if isempty(gutso), error('only rejected events in list'); end
  refEvIdx=intersect(refEvIdx,gutso); 
  disp(['out of a total of ' int2str(tmp) ', ' int2str(tmp-length(refEvIdx)) ' rejected events were excluded from latency analysis']);
end
% output
% 1st col: time stamp column of reference etsl
% 2nd and higher cols: latencies 
lat=[refetsl(:,1) repmat(nan,[size(refetsl,1) numOfLists-1])];
% initialize to NaN
minLat=lat(:,2);
idxIdx=lat(:,2);
for h=2:numOfLists,
  cmpetsl=varargin{h};
  cmpEvIdx=1:length(cmpetsl);
  tmp=cmpEvIdx(end);
  if size(cmpetsl,2)<etslc.tagCol
    disp(['ts list #' num2str(h-1) ' lacks column for tags, no check for rejected events possible']);
  else
    disp(['scanning list #' num2str(h-1) ' for rejected events..']);
    gutso=find(cmpetsl(:,etslc.tagCol)~=etslc.reject);
    % this might be improved - put out nans in the respective column
    if isempty(gutso), error('only rejected events in list'); end
    cmpEvIdx=intersect(cmpEvIdx,gutso); 
    disp(['out of a total of ' int2str(tmp) ', ' int2str(tmp-length(cmpEvIdx)) ' rejected events were excluded from latency analysis']);
  end
  % reinitialize tmpIdx and minLat
  idxIdx(:)=NaN;
  minLat(:)=NaN;
  % for each single event in reference list..
  for i=1:length(refEvIdx),
    % ..compute latencies to all events in comparison list..
    tdiff=refetsl(refEvIdx(i),etslc.tsCol)-cmpetsl(cmpEvIdx,etslc.tsCol);
    % ..and store latency and index of event which comes closest to it 
    [wattdenn,wodenn]=min(abs(tdiff));
    idxIdx(i)=wodenn;
    minLat(i)=tdiff(wodenn);
  end
  % identify those events in reference list which have one and the same event in 
  % comparison list as nearest..
  idx1=find(diff(idxIdx)==0);
  if ~isempty(idx1)
    for j=1:length(idx1)
      idx2=find(idxIdx==idxIdx(idx1(j)));
      [nix,idx3]=min(abs(minLat(idx2)));
      idx4=setdiff(idx2,idx2(idx3));
      % .. and set to nan all except the one with the smallest latency
      minLat(idx4)=NaN;
      idxIdx(idx4)=NaN;
    end  
  end  
  clear idx1 idx2 idx3 idx4;
  % now make to nan those events which are not within acceptable range
  nanIdx=find(((minLat<0.0) & minLat<=intv(1)) | (minLat>=0.0 & minLat>=intv(2)));
  minLat(nanIdx)=NaN;
  idxIdx(nanIdx)=NaN;
  % minLat now is OK
  lat(:,h)=minLat;
end     

