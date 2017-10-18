function tArr=etslsubdiv(etsl,silentEtsl,intv)
% ** function tArr=etslsubdiv(etsl,silentEtsl,intv)
% subdivides a time interval in which bursts and silent periods occur into
% smaller subintervals. Each subinterval will get three new (relative) time
% frames:
% i) pre-burst position: beginning of subinterval minus beginning of burst
% (if the subinterval 'touches' the burst with its end it is not counted as
% a pre-burst interval); values are by definition all negative
% ii) post-burst position: beginning of subinterval minus end of burst;
% values are by definition all positive
% iii) within-burst position: beginning of subinterval minus beginning of
% burst; values may be negative or positive
%
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT          DESCRIPTION
% etsl
% silentEtsl
% intv             scalar                length (ms) of elementary analysis
%                                        interval
%                         OR
%                  2 column array        start and stop times of analysis
%                                        intervals (ms, may overlap)
% ´
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT           DESCRIPTION
% tArr             4 columns, in following order:
%                  - pre-burst time of subintervals
%                  - post-burst time of subintervals
%                  - within-burst time of subintervals
%                  - length of preceding burst



% --- initialization, predefined vars
etslconst;
nCol=4;


% --- checks of input
if isempty(etsl) && isempty(silentEtsl)
  warning('both time stamp lists are empty - output will be empty');
  tArr=zeros(0,nCol);
  return
end

if numel(intv)==1
  error('automatic generation of analysis intervals not yet implemented, sorry');
else
  [n1 n2]=size(intv);
  if n2~=2
    error('input arg intv must be a scalar or a 2 column arr');
  end
  if n1<1
    error('input arg intv contains zero lines');
  end
  if any(diff(intv,1,2)<=0)
    error('at least one subinterval as defined in intv is zero or negative');
  end
  if any(diff(intv,1,1)<=0)
    error('subintervals as defined in intv must be sorted');
  end
  % finally, preallocate tArr 
  tArr=repmat(nan,n1,nCol);
end

% start by picking only columns of interest
etsl=etsl(:,[etslc.tsCol etslc.durCol]);
silentEtsl=silentEtsl(:,[etslc.tsCol etslc.durCol]);
% start of subintervals
siStartTsl=intv(:,1);
% end of subintervals
siEndTsl=intv(:,2);
nBu=size(etsl,1);
nSp=size(silentEtsl,1);

% continue depending on number of events:
if nBu==0 && nSp==1
  % one full silent period & no burst
  preIx=find(siStartTsl>silentEtsl(1,1) & siEndTsl<sum(silentEtsl(1,[1 2])));
  postIx=preIx;
  withinIx=find(siEndTsl>=sum(silentEtsl(1,[1 2])));
  % -- column 1: pre-burst time (ms)
  tArr(preIx,1)=siStartTsl(preIx)-sum(silentEtsl(1,[1 2]));
  % -- column 2: post-burst time (ms)
  tArr(postIx,2)=siStartTsl(postIx)-silentEtsl(1,1);
  % -- column 3: within-burst time (ms)
  tArr(withinIx,3)=siStartTsl(withinIx)-sum(silentEtsl(1,[1 2]));
  % -- column 4: remains undefined
  % -- (burst length undefined)
elseif nBu==1 && nSp==0
  % one full burst
  postIx=find(siStartTsl>sum(etsl(1,[1 2])));
  preIx=find(siEndTsl<etsl(1,1));
  withinIx=find(siStartTsl<=sum(etsl(1,[1 2])) & siEndTsl>=etsl(1,1));
  % -- column 1: pre-burst time (ms)
  tArr(preIx,1)=siStartTsl(preIx)-etsl(1,1);
  % -- column 2: post-burst time (ms)
  tArr(postIx,2)=siStartTsl(postIx)-sum(etsl(1,[1 2]));
  % -- column 3: within-burst time (ms)
  tArr(withinIx,3)=siStartTsl(withinIx)-etsl(1,1);  
  % -- column 4: pre burst len (ms)
  tArr(postIx,4)=etsl(1,2);
else
  % do it burst by burst
  for buIx=1:nBu
    if buIx==1
      if etsl(1,1)<silentEtsl(1,1)
        % if first burst comes before first silent period take all
        % pre-burst sweeps
        preIx=find(siEndTsl<etsl(buIx,1));
      else
        % all sweeps from beginning of first silent period
        preIx=find(siEndTsl<etsl(buIx,1) & siStartTsl>silentEtsl(1,1));
      end
      % post-burst sweeps:
      if buIx==nBu
        % if this is the only burst...
        postIx=find(siStartTsl>sum(etsl(buIx,[1 2])));
      else
        postIx=find(siStartTsl>sum(etsl(buIx,[1 2])) & siEndTsl<etsl(buIx+1,1));
      end
    elseif buIx==nBu
      if etsl(end,1)>silentEtsl(end,1)
        % if last burst comes after last silent period take all post-
        % burst ones
        postIx=find(siStartTsl>sum(etsl(buIx,[1 2])));
      else
        % take all up to end of last silent period
        postIx=find(siStartTsl>sum(etsl(buIx,[1 2])) & siEndTsl<sum(silentEtsl(end,[1 2])));
      end
      % pre-burst sweeps as usual
      preIx=find(siStartTsl>sum(etsl(buIx-1,[1 2])) & siEndTsl<etsl(buIx,1));
    else
      postIx=find(siStartTsl>sum(etsl(buIx,[1 2])) & siEndTsl<etsl(buIx+1,1));
      preIx=find(siStartTsl>sum(etsl(buIx-1,[1 2])) & siEndTsl<etsl(buIx,1));
    end
    withinIx=find(siStartTsl<=sum(etsl(buIx,[1 2])) & siEndTsl>=etsl(buIx,1));

    % -- column 1: pre-burst time (ms)
    tArr(preIx,1)=siStartTsl(preIx)-etsl(buIx,1);
    % -- column 2: post-burst time (ms)
    tArr(postIx,2)=siStartTsl(postIx)-sum(etsl(buIx,[1 2]));
    % -- column 3: within-burst time (ms, same reference frame as pre-burst time)
    tArr(withinIx,3)=siStartTsl(withinIx)-etsl(buIx,1);
    % -- column 4: pre burst len (ms)
    tArr(postIx,4)=etsl(buIx,2);
  end
end

% just to be on the safe side...(within- and peri-burst coverage must be
% mutually exclusive)
tmp=~isnan(tArr);
if any((tmp(:,1) | tmp(:,2)) & tmp(:,3))
  error('overlap of within- and peri-burst coverage');
end