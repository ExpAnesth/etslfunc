function varargout=utility_effectSize_peth(ds,varargin)
% ** function varargout=utility_effectSize_peth(ds,varargin) compares two
% populations of peri-event time histograms as computed by tslbatch or
% perievdeal. Bin-by-bin effect size measures are computed, and the following
% plots generated:
% i) average peths 
% ii) average peths scaled by their maxima 
% iii) the ratio of the unscaled, averaged peths
% iv) effect sizes including 95% confidence intervals
% v) burst prevalence (duration) during PE interval
% vi) in case of paired data: averaged normalized PETH differences
% (peth1-peth2)/(peth1+peth2)
% The function needs input structure ds (see
% template_utility_effectSize_peth.m)
% varagin, if specified, will be plugged into mes.m
% OUTPUTS: 
% - stats (the output of mes.m) to which raw, averaged and difference
% peths are appended as fields
% *** NOTE: the function had been designed for tslbatch, but was later
% adapted for perievdeal. In order to be able to deal with both data
% formats, the function features two different subfunctions responsible for
% reading data from the two analysis routines, and some code which is
% awkward or inefficient in terms of data handling but necessary for
% compatibility


% -------------------------------------------------------------------------
%                   I. LOAD, EXTRACT & RESHAPE DATA  
% -------------------------------------------------------------------------
% average electrodes within experiments??
ds.doAvEl=ischar(ds.curRecSite) && strcmpi(ds.curRecSite,'av');
if ds.doAvEl
  % if electrodes are to be averaged, this is the internal code
  ds.curRecSite=-1;
end

% data to be loaded from separate files?
% §§ include dirs
if iscell(ds.dFn)
  if numel(ds.dFn)==1
    doSeparateDataSets=false;
    ds.dFn=ds.dFn{1};
  elseif numel(ds.dFn)==2
    doSeparateDataSets=true;
    if ds.pair
      error('comparison of paired data sets in different files not allowed');
    end
    if numel(unique(ds.curIpVal))~=1
      butt=questdlg('Two different data sets will be compared at different values of ds.curIpVal - sure this is what you want?','Your valued attention, please');
      if ~strcmpi(butt,'yes');
        error('user break');
      end
    end
    % split second data set from first
    ds2=ds;
    ds2.dFn=ds2.dFn{2};
    ds2.dDir=ds2.dDir{2};
    ds2.curIpVal=ds2.curIpVal(2);
    ds.dFn=ds.dFn{1};
    ds.dDir=ds.dDir{1};
    ds.curIpVal=ds.curIpVal(1);
  else
    error('too many data files specified (one or two allowed)');
  end
else
  doSeparateDataSets=false;
end

% find out whether we're dealing with results files of tslbatch or of
% perievdeal (the first file must suffice)
s=whos('-file',[ds.dDir ds.dFn]);
if any(strcmp({s.name},'dataInfo'))
  loadFun=@loadPerievdealData;
else
  loadFun=@loadTslbatchData;
end

% load (first) data file
[recSite,expChanName1,pethBin,pethMn,pethN,binW,columnIx]=loadFun(ds);

% extract first set of data
peth1=permute(pethMn(:,columnIx(1),:),[1 3 2]);
ix1=any(isfinite(peth1));
pethN1=permute(pethN(:,columnIx(1),:),[1 3 2]);
recSite1=permute(recSite(1,columnIx(1),:),[1 3 2]);
% pick only data obtained from specific recording site
ix1= ix1 & ismember(recSite1,ds.curRecSite);
ix1=find(ix1);

% load second data file and/or extract second set of data
if doSeparateDataSets
  [recSite,expChanName2,pethBin2,pethMn,pethN,binW,columnIx2]=loadFun(ds2);
  columnIx=[columnIx columnIx2];
  if ~isequal(pethBin,pethBin2)
    error('data sets with different bins');
  end
end

peth2=permute(pethMn(:,columnIx(2),:),[1 3 2]);
ix2=any(isfinite(peth2));
pethN2=permute(pethN(:,columnIx(2),:),[1 3 2]);
recSite2=permute(recSite(1,columnIx(2),:),[1 3 2]);

if ds.pair
  ix2=ix2 & ismember(recSite2,ds.curRecSite);
  ix2=intersect(find(ix2),ix1);
  ix1=ix2;
  nExp=numel(ix2);
else
  ix2=ix2 & ismember(recSite2,ds.curRecSite);
  ix2=find(ix2);
  nExp=[numel(ix1) numel(ix2)];
end

clear recSite pethBin2 columnIx* pethMn pethN
% -------------------------------------------------------------------------
%                   III. COMPUTE & PLOT
% -------------------------------------------------------------------------

if ds.indvPlot
  % --- overview of individual experiments fig
  if ds.pair
    % paired data: plot together in one axis
    figure;
    nRow=ceil(sqrt(nExp));
    nCol=ceil(nExp/nRow);
    for h=1:nExp
      subplot(nRow,nCol,h),
      hold on
      set(gca,'fontsize',7);
      % control
      plot(pethBin+binW/2,peth1(:,ix1(h)),'k');
      % applic
      plot(pethBin+binW/2,peth2(:,ix1(h)),'r');
      niceyuax;
      smarttext(['ID: ' int2str(ix1(h))],.9,.92,'fontsize',8);
      if iscell(expChanName1)
        title(expChanName1{ix1(h)},'interpreter','none');
      else
        title(expChanName1(ix1(h)).expID,'interpreter','none');
      end
      % mark if indexed for omission
      if ismember(ix1(h),union(ds.killer(1).ix,ds.killer(2).ix))
        set(gca,'color',[.6 .6 .6])
      end
    end
  else
    % unpaired data
    for g=1:2
      figure;
      nRow=ceil(sqrt(nExp(g)));
      nCol=ceil(nExp(g)/nRow);
      peth=eval(['peth' int2str(g) ';']);
      ix=eval(['ix' int2str(g) ';']);
      if g==2 && doSeparateDataSets
        expChanName=expChanName2;
      else
        expChanName=expChanName1;
      end
      for h=1:nExp(g)
        subplot(nRow,nCol,h),
        set(gca,'fontsize',7);
        plot(pethBin+binW/2,peth(:,ix(h)),'k');
        niceyuax;
        smarttext(['ID: ' int2str(ix(h))],.9,.92,'fontsize',8);
        if iscell(expChanName)
          title(expChanName{ix(h)},'interpreter','none');
        else
          title(expChanName(ix(h)).expID,'interpreter','none');
        end
        % mark if indexed for omission
        if ismember(ix(h),ds.killer(g).ix)
          set(gca,'color',[.6 .6 .6])
        end
      end
    end
  end
  clear peth ix expChanName nRow nCol h g 
end

% *** now kill selected IDs
if ds.pair
  ix2=setdiff(ix2,union(ds.killer(1).ix,ds.killer(2).ix));
  ix1=ix2;
  nExp=numel(ix2);
else
  ix1=setdiff(ix1,ds.killer(1).ix);  
  ix2=setdiff(ix2,ds.killer(2).ix);
  nExp=[numel(ix1) numel(ix2)];
end

% compute averages (for output struct and plots)
mnPeth1=ds.avFuncH(peth1(:,ix1),2);
mnPeth2=ds.avFuncH(peth2(:,ix2),2);

% compute difference PETHs (only for paired data)
if ds.pair
  diffPETH=(peth1(:,ix1)-peth2(:,ix2))./(peth1(:,ix1)+peth2(:,ix2));
  
  diffPETH=peth1(:,ix1)-peth2(:,ix2);
  
  %   % normalize: divide each bin of the difference peth by sum of both
  %   diffPETH=diffPETH./repmat(nansum(peth2(:,ix2)+peth1(:,ix1)),size(peth1,1),1);
    
  % diffPETH=diffPETH./repmat(nansum(abs(diffPETH)),size(peth1,1),1);
  % diffPETH=diffPETH./peth1(:,ix1);
  
  mnDiffPETH=nanmedian(diffPETH,2);
  stdDiffPETH=prctile(diffPETH,[25 75],2);
else
  diffPETH=[];
end

% *** compute effect size
% transpose peth2 and peth1: experiments down the columns (paired or not),
% time along rows (columns must match anyways)
stats=mes(peth2(:,ix2)',peth1(:,ix1)',...
  ds.esType,'isDep',ds.pair,'missVal','pairwise','doPlot',false,varargin{:});
es=stats.(ds.esType);
ci=stats.([ds.esType 'Ci']);


if ds.doPlot
  % *** main plot
  tmpScrSz=get(0,'Screensize');
  tmpScrSz([1 2])=tmpScrSz([3 4])*.05;
  tmpScrSz([3 4])=tmpScrSz([3 4]).*[.5 .85];
  fh=figure; orient tall, clf
  set(fh,'position',tmpScrSz)
  % --- plot PETH
  subplot(4,2,1), hold on
  % first set
  ph=stairs(pethBin,mnPeth1);
  set(ph,'color',ds.pCol(1,:),'linewidth',1.5);
  % second set
  ph=stairs(pethBin,mnPeth2);
  set(ph,'color',ds.pCol(2,:),'linewidth',1.5);
  if isfinite(ds.pethYLim)
    set(gca,'ylim',[0 ds.pethYLim]);
  else
    nicexy0ax(30);
  end
  set(gca,'xlim',ds.intv);
  ylabel('frequency (Hz)');
  title(['indep par levels=' num2str(ds.curIpVal) '; n=' int2str(nExp)]);
  
  % --- plot of PETH, normalized to peak
  mnPeth1_norm=mnPeth1/max(mnPeth1);
  mnPeth2_norm=mnPeth2/max(mnPeth2);
  subplot(4,2,2), hold on
  % first set
  ph=stairs(pethBin,mnPeth1_norm);
  set(ph,'color',ds.pCol(1,:),'linewidth',1.5);
  % second set
  ph=stairs(pethBin,mnPeth2_norm);
  set(ph,'color',ds.pCol(2,:),'linewidth',1.5);
  nicexy0ax(30);
  set(gca,'xlim',ds.intv);
  ylabel('frequency (normalized)');
  
  % --- effect size plot
  subplot(4,2,3), hold on
  switch ds.mesPlotStyle
    case 'errorbar'
      ph=plot(pethBin+binW/2,es','ko');
      set(ph,'markerfacecolor','k');
      lh=line([(pethBin+binW/2)'; (pethBin+binW/2)'],ci);
      set(lh,'color','k')
    case 'boundedline'
      [lh,ph]=boundedline(pethBin+binW/2,es',abs(ci-es)','k');
    case 'area'
      ah=area(pethBin+binW/2,es',ds.NullEffectVal,'FaceColor','k');
      lh=line(pethBin+binW/2,ci','linestyle','-','color',[.5 .5 .5]);
    otherwise
      error('illegal choice of ds.mesPlotStyle')
  end
  if all(isfinite(ds.esYLim))
    set(gca,'ylim',ds.esYLim);
  else
    nicexyax(30);
  end
  set(gca,'xlim',ds.intv);
  set(gca,'ygrid','on')
  line(ds.intv,ds.NullEffectVal*[1 1],'linestyle','--','color','k','linewidth',get(gca,'defaultlinelinewidth')*1.5);
  xlabel('time (ms)');
  ylabel('effect size');
  
  % --- plot of burst prevalence (duration) during PE interval
  subplot(4,2,5)
  hold on
  ph=plot(pethBin,mean(pethN1(:,ix1),2));
  set(ph,'color',ds.pCol(1,:));
  ph=plot(pethBin,mean(pethN2(:,ix2),2));
  set(ph,'color',ds.pCol(2,:));
  set(gca,'xlim',ds.intv,'ylim',[0 1.05]);
  xlabel('time (ms)');
  ylabel('burst prevalence');
  
  if ds.pair
    % --- plot of averaged normalized PETH differences
    subplot(4,2,6), hold on
    ph=plot(pethBin+binW/2,mnDiffPETH,'ko');
    set(ph,'markerfacecolor','k');
    % lh=line([(pethBin+binW/2)'; (pethBin+binW/2)'],(mnDiffPETH*[1 1])'+(stdDiffPETH*[-1 1])');
    lh=line([(pethBin+binW/2)'; (pethBin+binW/2)'],(stdDiffPETH)');
    set(lh,'color','k')
    nicexyax(30);
    set(gca,'xlim',ds.intv);
    set(gca,'ygrid','on')
    line(ds.intv,[0 0],'linestyle','--','color','k','linewidth',get(gca,'defaultlinelinewidth')*1.5);
    xlabel('time (ms)');
    ylabel('av. norm. diff. PETH');
  end
  
  % there are folks who prefer cumulative histograms, so let's plot them
  subplot(4,2,7), hold on
  % first set
  ph=stairs(pethBin,cumsum(mnPeth1));
  set(ph,'color',ds.pCol(1,:),'linewidth',1.5);
  % second set
  ph=stairs(pethBin,cumsum(mnPeth2));
  set(ph,'color',ds.pCol(2,:),'linewidth',1.5);
  nicexy0ax(30);
  set(gca,'xlim',ds.intv);
  xlabel('time (ms)')
  ylabel('cum. frequency (Hz)');
  
  
  if ~isempty(ds.printas)
    % shave any file extension off ds.exportFn, if any
    [~,ds.exportFn]=fileparts(ds.exportFn);
    print(ds.printas,'-r400',[ds.exportDir ds.exportFn]);
  end
  
  % figdefault
end


% outputs:
if nargout>0
  % attach bins & peths to stats
  stats.pethBin=pethBin;
  stats.pethBinW=binW;
  stats.mnPeth1=mnPeth1;
  stats.mnPeth2=mnPeth2;
  stats.peth1=peth1;
  stats.peth2=peth2;
  stats.diffPeth=diffPETH;
  varargout{1}=stats;
end

% ======================= LOCAL FUNCTIONS =================================

function [recSite,expChanName,pethBin,pethMn,pethN,binW,columnIx]=loadPerievdealData(ds)

% load parameters needed independent of ds.typeEvent  
load([ds.dDir ds.dFn],'indepPa*','fileInfo','wp_s','peth');

switch lower(ds.typeEvent)
  case 'start'
    pethBin=wp_s.pethBin;
    % the only two possibilities
    fina={'nil_pe_bu','nil_pe_bu_pc'};
    ix=isfield(peth,fina);
    if all(ix)
      choice=questdlg('Which shall be the reference channel (pc=principal channel)?', ...
        'PETH Menu', ...
        fina{:},fina{1});
      fina={choice};
      ix=1;
    end
    pethMnCell=peth.(fina{ix}).spxHistMn;
    pethNCell=peth.(fina{ix}).spxHistN;
    
  case 'end'
    pethBin=wp_s.spPethBin;
    % the only two possibilities
    fina={'nil_pe_sp','nil_pe_sp_pc'};
    ix=isfield(peth,fina);
    if all(ix)
      choice=questdlg('What would you like on your birthday cake?', ...
        'Cake Menu', ...
        fina{:},fina{1});
      fina={choice};
      ix=1;
    end
    pethMnCell=peth.(fina{ix}).spxHistMn;
    pethNCell=peth.(fina{ix}).spxHistN;

  otherwise
    error('illegal value for ds.typeEvent (choose ''start'' or ''end''');
end

% size of things
[n1 n2 n3]=size(pethMnCell);
nBin=numel(pethMnCell{1});

% index to columns (=independent variable)
[nix,columnIx,tmpIx]=intersect(wp_s.indepParLevel,ds.curIpVal(:));
if ~isequal(nix,unique(ds.curIpVal(:)))
  errordlg({'variable ''ds.curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(wp_s.indepParLevel)})
  return
end

% revert sorting by intersect, if any
columnIx=columnIx(tmpIx);
pethMn=[];
pethN=[];
recSite=[];
% here, loadPerievdealData deliberately differs from loadTslbatchData
expChanName=fileInfo(1,1,:);
if ds.doAvEl
  % preallocate output vars depending on choice of channels etc
  pethMn=nan([nBin n2 n3]);
  pethN=nan([nBin n2 n3]);
  % pick peths columnwise from cell arrays, average, embed
  for g=1:n2*n3
    pethMn(:,g)=nanmean(cat(2,pethMnCell{:,g}),2);
    pethN(:,g)=nanmean(cat(2,pethNCell{:,g}),2);
  end
  % set values: internal code for averages
  recSite=repmat(ds.curRecSite,[1 n2 n3]);
else
  % loop thru layers of pethMnCell, cat in third dimension
  for g=1:n1
    pethMn=cat(3,pethMn,reshape(cat(2,pethMnCell{g,:,:}),[nBin n2 n3]));
    pethN=cat(3,pethN,reshape(cat(2,pethNCell{g,:,:}),[nBin n2 n3]));
    % recSite
    recSite=cat(3,recSite,repmat(wp_s.uRecSite(g),[1 n2 n3]));
  end
  expChanName=repmat(expChanName,[1 1 n1]);
end
   
  
% index to original bins
binIx=pethBin>=ds.intv(1) & pethBin<=ds.intv(2);
pethBin=pethBin(binIx);
% MINIMAL original bin width
binW=min(diff(pethBin));

% normalize pethN so that it reflects the proportion of bursts which
% are 'on' at each bin
pethN=pethN./repmat(max(pethN),[size(pethN,1) 1]);

% restrict PETH variables to set of bins under consideration
pethMn=pethMn(binIx,:,:);
pethN=pethN(binIx,:,:);


function [recSite,expChanName,pethBin,pethMn,pethN,binW,columnIx]=loadTslbatchData(ds)
% load parameters needed independent of ds.typeEvent  
% ** note: this part of the function tries to accomodate data generated by
% really old variants of tslbatch
load([ds.dDir ds.dFn],'indepPa*','recSite','undef*','expChanName');
switch lower(ds.typeEvent)
  case 'start'
    load([ds.dDir ds.dFn],'pethMn','pethBin','pethN');
    % tslbatch version issue: pethBin was formerly called 'binArr'
    if ~exist('pethBin','var')
      load([ds.dDir ds.dFn],'binArr');
      if exist('binArr','var')
        pethBin=binArr;
      else
        error('the data file does not contain proper peth data');
      end
    end
  case 'end'
    load([ds.dDir ds.dFn],'spPethMn','spPethBin','spPethN');
    % check right away whether the desired peth data exist
    if ~exist('spPethBin','var')
      error('the data file does not contain peth data for burst ends');
    end
    pethBin=spPethBin;
    pethMn=spPethMn;
    pethN=spPethN;    
  otherwise
    error('illegal value for ds.typeEvent (choose ''start'' or ''end''');
end

% §§ here's to more version issues, darn it...
if ~exist('recSite','var')
  if ~exist('undef1','var')
    warning('data set does not contain information on recording site or values of column ''undef1'' - this is ONLY of concern if you want the data set restricted according to this criterion');
    recSite=repmat(ds.curRecSite,[1 size(pethMn,2) size(pethMn,3)]);
  else
    recSite=undef1;
  end
end
if ~exist('pethN','var')
  pethN=ones(size(pethMn));
end
  
% index to columns (=independent variable)
[nix,columnIx,tmpIx]=intersect(indepParLevel,ds.curIpVal(:));
if ~isequal(nix,unique(ds.curIpVal(:)))
  errordlg({'variable ''ds.curIpVal'' contains values that do not exist in the data.',...
    'Available values: ', num2str(indepParLevel)})
  return
end
% revert sorting by intersect, if any
columnIx=columnIx(tmpIx);
% index to original bins
binIx=pethBin>=ds.intv(1) & pethBin<=ds.intv(2);
pethBin=pethBin(binIx);
% MINIMAL original bin width
binW=min(diff(pethBin));

% normalize pethN so that it reflects the proportion of bursts which
% are 'on' at each bin
pethN=pethN./repmat(max(pethN),[size(pethN,1) 1]);

% restrict PETH variables to set of bins under consideration
pethMn=pethMn(binIx,:,:);
pethN=pethN(binIx,:,:);


% average electrodes from same culture?
if ds.doAvEl
  nECh=numel(expChanName);
  expName={};
  pethMn2=nan(size(pethMn));
  pethN2=nan(size(pethN));  
  ct=1;
  list=1:nECh;
  while ~isempty(list)
    curExp=expChanName{list(1)};
    commaIx=strfind(curExp,',');
    expName{ct}=curExp(1:commaIx-1);
    % index to entries of current exp
    ix=find(strncmp(curExp(1:commaIx-1),expChanName,commaIx-1));
    pethMn2(:,:,ct)=nanmean(pethMn(:,:,ix),3);
    pethN2(:,:,ct)=nanmean(pethN(:,:,ix),3);    
    list=setdiff(list,ix);
    ct=ct+1;
  end
  % overwrite, killing the superfluous preallocated rows
  pethMn=pethMn2(:,:,1:ct-1);
  pethN=pethN2(:,:,1:ct-1);
  expChanName=expName;
  % replace values by internal code for averages
  recSite=repmat(ds.curRecSite,size(pethMn,2),ct-1);
end
