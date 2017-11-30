function adr=etslAnomalyDetection(d,transIA,transAI,varargin)
% ** 
%                    >>> INPUT VARIABLES >>>
%
% NAME         TYPE/DEFAULT      DESCRIPTION
% d            1-column array    raw data (time series)
% transIA      col array         list of transitions i->a
% transAI      col array         list of transitions a->i
% funHandle    fHandle, {@mean,@std}      
% detMethod    char,'kmeans'     method: 'Gaussian','kmeans','explicit'
% thresh       array             threshold(s) for explicit method
% pThresh      scalar, .01       threshold for Gaussian method
% axH          axis handle       handle to axis for plot (no plot will be
%                                produced if empty or not specified)
%                    <<< OUTPUT VARIABLES <<<
%
% NAME         TYPE/DEFAULT      DESCRIPTION
% adr          struct            .funVal: function values
%                                .isAnomTransIA: result of classification
%                                .isAnomTransAI: result of classification
% **********************************************************************
% * ALL INPUT ARGUMENTS MUST BE SPECIFIED IN IDENTICAL TIME UNITS      *
% * OUTPUT ARGUMENTS WILL HAVE SAME TIME UNITS                         *
% **********************************************************************
%             - work in progress - 

% to do:
% rational choice of b

% ------ defaults

funHandle={@mean,@std};
detMethod='kmeans';
thresh=[];
pThresh=.01;
axH=[];
pvpmod(varargin,{'funHandle','detMethod','thresh','pThresh','axH'});

% convert to cell
if ~iscell(funHandle)
  funHandle={funHandle};
end
numQuantifier=numel(funHandle);
etslconst;

% ------ checks/prep work
[n1,n2,n3]=size(d);
if n2>1 || n3>1
  error('d must be a single column array');
end

% reject too short lists 
if numel(transIA)<2 || numel(transAI)<2
  warning('too few events - quantifying events but not running anomaly detection')
  adr.isAnomTransIA=false(size(transIA));
  adr.isAnomTransAI=false(size(transAI));
  % adr.funVal=nan;
  detMethod='none';
end

if strcmp(detMethod,'explicit')
  if numel(thresh) ~= numQuantifier
    error('in the ''explicit'' method the number of thresholds must equal the number of function handles')
  end
end

% ------ THE WORKS
% find combinations of transitions amounting to full events (that is, the
% i-a transitions followed by a-i transitions
transIAIx=true(size(transIA));
transAIIx=true(size(transAI));
% sequence starts with a-i transition: get rid of it
if transAI(1)<transIA(1)
  transAIIx(1)=false;
end
% sequence ends with a i-a transition: get rid of it
if transIA(end)>transAI(end)
  transIAIx(end)=false;
end
% now convert to linear indexes 
transIAIx=find(transIAIx);
transAIIx=find(transAIIx);

% check to make sure
numEv=unique([length(transIAIx),length(transAIIx)]);
if numel(numEv)>1
  error('internal: transIAIx and transAIIx don''t match')
end

% preallocate logical arrays pointing to anomalous events
adr.isAnomTransIA=false(size(transIA));
adr.isAnomTransAI=false(size(transAI));
% preallocate array for quantifiers (function values)
adr.funVal=nan(numEv,numQuantifier);

% loop over events
for k=1:numEv
  dExc=d(transIA(transIAIx(k)):transAI(transAIIx(k)));
  for g=1:numQuantifier
    adr.funVal(k,g)=funHandle{g}(dExc);
  end
end

% default: none is anomalous
badIx=false(size(adr.funVal,1),1);
switch detMethod
  case 'kmeans'
    numCluster=max(2,min(8,numEv-1));
    [cluIx,cluCentroid]=kmeans(adr.funVal,numCluster);
    % identify cluster with smallest mean: this is the poor one
    [~,badCluIx]=min(cluCentroid);
    badCluIx=unique(badCluIx,'stable');
    if numel(badCluIx)==1
      badIx=cluIx==badCluIx;
    else
      disp('anomalous samples cannot be identified because none of the clusters can be identified as the bad one')
    end
    
  case 'Gaussian'
    numCluster=2;
    try
      gmModel=fitgmdist(adr.funVal,numCluster,'Options',statset('Display','final','MaxIter',50,'TolFun',0.001));
      % identify cluster with larger means: these are the good ones
      [mainCluMn,mainCluIx]=max(gmModel.mu);
      mainCluIx=unique(mainCluIx,'stable');
    catch
      mainCluIx=[];
    end
    
    if numel(mainCluIx)==1
      % compute posterior prob of belonging to different clusters
      probCluster=posterior(gmModel,adr.funVal);
      % note that 'bad' bursts are those with a low probability of belonging to
      % the main cluster AND with values BELOW those of the main cluster
      badIx=probCluster(:,mainCluIx)<=pThresh & any(adr.funVal<mainCluMn,2);
    elseif isempty(mainCluIx)
      disp('fitting Gaussian Mixture Model failed')
    else
      disp('anomalous samples cannot be identified because none of the Gaussians can be identified as the main one')
    end
    
  case 'explicit'
    badIx=false(size(adr.funVal,1),1);
    % note definition of bad=smaller than threshold in ANY dimension!
    badIx=badIx | any(adr.funVal<thresh(:)',2);
    
  case 'none'
    disp('no anomaly detection requested')
    
  otherwise
    error('illegal method')
end
    

% linear indexes to transitions belonging to anomalous bursts
adr.isAnomTransIA(transIAIx(badIx))=true;
adr.isAnomTransAI(transAIIx(badIx))=true;

% plot
if ~isempty(axH) && isgraphics(axH)
  axes(axH)
  if numQuantifier>=2
    % currently, 2D scatterplot
    if numQuantifier==2
      markerSize=10;
    else
      markerSize=adr.funVal(:,3)*30/max(adr.funVal(:,3));
    end
    cla
    scatter(adr.funVal(:,1),adr.funVal(:,2),markerSize,[double(badIx)*.8 zeros(numEv,2)]);
    nicexyax
    grid on
    xlabel(char(funHandle{1}))
    ylabel(char(funHandle{2}))
  else
    markerSize=10;
    histogram(adr.funVal(:,1),40);
    axis tight
    yl=get(gca,'ylim');
    hold on
    scatter(adr.funVal(:,1),rand(numEv,1)*yl(2),markerSize,[double(badIx)*.8 zeros(numEv,2)]);
    xlabel(char(funHandle{1}))
    hold off
  end
end
