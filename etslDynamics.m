function stats=etslDynamics(aEtsl,iEtsl,order,varargin)
% ** function stats=etslDynamics(aEtsl,iEtsl,order,varargin)
% analyzes dynamics of active and inactive events 
%
%                    >>> INPUT VARIABLES >>>
% NAME         TYPE/DEFAULT        DESCRIPTION
% aEtsl        2D arr              extended time stamp list of 'active'
%                                   events
% iEtsl        2D arr              extended time stamp list of 'inactive'
%                                   events
% order        char                any of 'aa','ai','ia','ii'
% doPlot       logical, false      if true, plots return map
% 
%                    <<< OUTPUT VARIABLES <<<
% NAME         TYPE/DEFAULT           DESCRIPTION
% d            2-column array                 

doPlot=false;
pvpmod(varargin);

% check input args
if numel(order)~=2 || any(~ismember(order,{'aa','ai','ia','ii'}))
  error('input argument ''order'' must be any combination of two letters ''a'' and ''i''');
end

etslconst;

nA=size(aEtsl,1);
nI=size(iEtsl,1);


% a 2-column array with durations of the first type of event in the first
% column and durations of the FOLLOWING second type of event in the second
% column
stats.d=zeros(0,2); 
stats.r=nan;

% tmp=diff(aEtsl
switch order

  case 'aa'
    if nA
      stats.d=[aEtsl(1:end-1,etslc.durCol) aEtsl(2:end,etslc.durCol)];
    end
      
  case 'ii'
    if nI
      stats.d=[iEtsl(1:end-1,etslc.durCol) iEtsl(2:end,etslc.durCol)];
    end
    
  case {'ai','ia'}
    if strcmp(order,'ai')
      % probably the most interesting case: duration of inactive event
      % following an active event
      e1=aEtsl;
      e2=iEtsl;
      n1=nA;
      n2=nI;
    else
      e1=iEtsl;
      e2=aEtsl;
      n1=nI;
      n2=nA;
    end
    if n1 && n2
      if e1(1,etslc.tsCol)<e2(1,etslc.tsCol)
        % first full event is a
        nEv=min(n1,n2);
        stats.d=[e1(1:nEv,etslc.durCol) e2(1:nEv,etslc.durCol)];
      else
        % first full event is i
        nEv=min(n1,n2-1);
        stats.d=[e1(1:nEv,etslc.durCol) e2(2:nEv+1,etslc.durCol)];
      end        
    end
end

if size(stats.d,1)>1
  stats.r=corrcoef(stats.d);
  stats.r=stats.r(2);
end

if doPlot
  if size(stats.d,1)>1
    plot(stats.d(:,1),stats.d(:,2),'ko')
    nicexyax;
    xlabel(order(1));
    ylabel(order(2));
    smarttext(num2str(stats.r));
  else
    smarttext('too few events');
  end
end