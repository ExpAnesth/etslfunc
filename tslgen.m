function tsl=tslgen(pdf,varargin)
% ** function tsl=tslgen(pdf,varargin)
% generates artificial time stamp lists (tsl). 
% There are two - mutually exclusive - ways to set the length of the tsl: 
% a) by specifying the number of time stamps. This is the simple option
% b) by specifying 
%                         >>> INPUT VARIABLES >>>
% NAME             TYPE,DEFAULT                 DESCRIPTION
% pdf              2d-array                     the probability density function from which the ISI are drawn.
%                                               1st column: inter-event-time (ms)
%                                               2nd column: probability density
% nts              scalar,[]                    the length of the tsl: number of time stamps
% iLen             scalar, []                   max length of the tsl in terms of time (ms) - this option overrides 'nts'
% si               scalar, 100                  the sampling interval of the ts
% verbose          scalar, 1                    if nonzero, all sorts of information about the process of tsl generation will be splashed on the screen
%
%                         <<< OUTPUT VARIABLES <<<
% NAME             TYPE/DEFAULT          DESCRIPTION
% tsl              column array          the time stamp list - note that the first ts is always at time 0
%
%
% 
%  ------- non-standard matlab functions used:
% 
% --------------------------------------------
% default values
%

% to do
% - rand
% - cutting down tsl (de-nan)
% - introduce some jitter in spike times


nts=[];
iLen=[];
verbose=1;

pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;
% check & normalize pdf
pdfSz=size(pdf);
if (pdfSz(1)>=1 & pdfSz(2)==2)
  % no negative probabilities or isi, please
  if find(pdf<0.0)
    error('check pdf - negative value found');
  end
  pdfBinW=unique(diff(pdf(:,1)));
  if ~isequal(size(pdfBinW),[1 1])
    error('abscissa step size of pdf must be unique');
  end
  % 3rd col of pdf=cumulative prob
  pdf(:,3)=cumsum(pdf(:,2)/pdfSz(1));
  if pdf(end,3)~=1.0
    disp('pdf is not normalized - doing this now');
    pdf(:,[2 3])=pdf(:,[2 3])/pdf(end,3);
    % mean of (isi*prob(isi)) corresponds to the expected (mean) isi
    isiMn=mean(pdf(:,2).*pdf(:,1));
    % the expected event rate in Hz
    frMn=1/isiMn*1000;
  end
  if verbose,
    figure(99), clf
    plot(pdf(:,1),pdf(:,[2 3])); 
    niceyuax; xlabel('IEI (ms)');ylabel('prob | normalized freq');
    ultext(['mean event rate of normalized pdf: ' num2str(frMn)]);
  end
else
  error('pdf must be an array of 2 columns')
end


% determine whether length of tsl is given in number of events or in terms of max time
if ~isempty(nts)
  if ~isempty(iLen)
    error('length of tsl must be specified EITHER in terms of number of events (nts) OR in terms of time (ilen), NOT both');
  else
    % the 'expected' number of events is in fact not estimated but predetermined in this case
    expNumTs=nts;
  end
elseif ~isempty(iLen)
  % estimate the expected number of events..
  expNumTs=frMn*iLen/1000;
  % ..and multiply by a factor ensuring that it will really be sufficiently high
  expNumTs=round(expNumTs*1.1*((expNumTs+10)/expNumTs));
end
  

% generate

tmp=rand(expNumTs,1);
tsl=nan*zeros(expNumTs,1);
for i=1:expNumTs
  [y,idx]=min(abs(tmp(i)-pdf(:,3)));
  tsl(i)=pdf(idx,1);
end
tsl=cumsum(tsl);
% cut down to given length
if ~isempty(iLen)
  tmp=find(tsl>iLen);
  if isempty(tmp)
    warning('internal problem: tsl does probably not contain enough ts');
  else
    tsl(tmp)=[];
  end
end

if verbose
  disp(['expected event rate: ' num2str(frMn)]);
  disp(['actual event rate  : ' num2str(1000*length(tsl)/tsl(end))]);
end

% mean of the pdf must correspond to the instantaneous ISI






