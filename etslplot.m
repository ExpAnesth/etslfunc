function ph=etslplot(etslCell,varargin)
% ** function function ph=etslplot(etslCell,varargin)
% creates plot of extended time stamp lists (as many as there are in
% etslCell).
% All input parameters except etslCell are optional and must be specified
% as parameter/value pairs, e.g. as in
%          etslplot(etslCell,'timeUnit','s')
%
%                         >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT              DESCRIPTION
% etslCell                                 cell array of etsls
% style          char, 'horizontalLines'   plot style, 'hmm' is the other
%                                          option
% timeUnit       char, 's'
% hOffset        double, 0                 horizontal offset of plotted
%                                          data
%
%                         <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT           DESCRIPTION
% ph                                    array or struct of handles to plot
%                                       objects



style='horizontalLines';
timeUnit='ms';
hOffset=0;
pvpmod(varargin)

% if etslCell is not a cell, assume that user mistakenly entered an etsl
% (=array) instead of a cell array
if ~iscell(etslCell)
  etslCell={etslCell};
  nEtsl=1;
else
  nEtsl=numel(etslCell);
end
etslconst

switch style
  case 'horizontalLines'
    ph.lineHandles=[];
  case 'hmm'
    ph=gobjects(nEtsl,1);
  otherwise
    error('bad style')
end


hold on
count=0;
for g=1:nEtsl
  etsl=etslCell{g};
  switch timeUnit
    case 'ms'
      % everything is fine
    case 's'
      etsl(:,[etslc.tsCol etslc.durCol])=etsl(:,[etslc.tsCol etslc.durCol])/1000;
    otherwise
      error('check timeUnit')
  end
  [n1 n2]=size(etsl);
  if n1>0 && n2>=max([etslc.tsCol etslc.durCol])
    count=count+1;
    % start and stop times of active periods in ms
    aTimes=cumsum(etsl(:,[etslc.tsCol etslc.durCol]),2);
    
    switch style
      case 'horizontalLines'
        offs=1-count+hOffset;
        ph(g).lineHandles=line(aTimes',repmat(offs,[2 size(aTimes,1)]),'color','k','linewidth',2);
        
      case 'hmm'
        avDur=median(etsl(:,etslc.durCol));
        xArr=aTimes(:,[1 1 2 2])';
        xArr=xArr(:);
        yArr=repmat([0 1 1 0]',n1,1);
        % add a little jitter in case of multiple etsl
        if count>1
          yArr=yArr+rand(size(yArr))*.02;
          xArr=xArr+rand(size(xArr))*avDur*.02;
        end
        ph(g)=plot(xArr,yArr+hOffset,'-o');
    end
  else
    disp('etsl empty or missing column')
  end
end
hold off
niceyax;