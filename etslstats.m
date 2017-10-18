function stats=etslstats(aEtsl,iEtsl)
% ** function stats=etslstats(aEtsl,iEtsl)
% computes a handful of pretty basic statistics from etsl, like median
% burst length, relative time spent in active or inactive periods, etc.
%
%                    >>> INPUT VARIABLES >>>
% NAME         TYPE/DEFAULT        DESCRIPTION
% aEtsl        2D arr              extended time stamp list of 'active'
%                                   events
% iEtsl        2D arr              extended time stamp list of 'inactive'
%                                   events
%
%                    <<< OUTPUT VARIABLES <<<
% NAME         TYPE/DEFAULT           DESCRIPTION
% stats        struct                 a few basic statistics


etslconst;

nA=size(aEtsl,1);
nI=size(iEtsl,1);

stats.burstRate=nan;
stats.relTimeInBurst=nan;

stats.mnBurstPeriod=nan;
stats.stdBurstPeriod=nan;
stats.geomnBurstPeriod=nan;
stats.mdBurstPeriod=nan;

stats.minBurstLen=nan;
stats.lqBurstLen=nan; % lower quartile
stats.mdBurstLen=nan;
stats.uqBurstLen=nan; % upper quartile
stats.maxBurstLen=nan;
stats.asBurstLen=nan; % asymmetry

stats.mnBurstLen=nan;
stats.stdBurstLen=nan;
stats.geomnBurstLen=nan;

stats.minSilentPerLen=nan;
stats.lqSilentPerLen=nan; % lower quartile
stats.mdSilentPerLen=nan;
stats.uqSilentPerLen=nan; % upper quartile
stats.maxSilentPerLen=nan;
stats.asSilentPerLen=nan; % asymmetry

stats.mnSilentPerLen=nan;
stats.stdSilentPerLen=nan;
stats.geomnSilentPerLen=nan;

% catch errors in some aEtsl and iEtsl spotted May 2015
if nI
  allBadIx=find((iEtsl(:,etslc.durCol))<=0);
  allBadIx2=allBadIx;
  if any(allBadIx)
    warning('faulty iEtsl - correcting')
    for g=1:numel(allBadIx)
      badIx=allBadIx(g);
      if nA>=badIx
        if badIx==1
          allBadIx2(g)=1;
        else
          [m,localBadIx]= min(abs(aEtsl(badIx-1:badIx,etslc.durCol) + iEtsl(badIx,etslc.durCol)));
          allBadIx2(g)=badIx+localBadIx-2;
        end
      else
        error('giving up')
      end
    end
    % now delete entries and the matching bursts in aEtsl
    iEtsl(allBadIx,:)=[];
    nI=nI-g;
    aEtsl(allBadIx2,:)=[];
    nA=nA-g;
  end
end

% compute time interval containing complete periods (i.e. start of
% first burst to start of last burst, or start of first silent
% periods to start of last silent period).
if nA && nI
  if aEtsl(1,etslc.tsCol)<iEtsl(1,etslc.tsCol)
    % if first transition is i->a, one way of defining cIntv is to
    % count from beginning of first burst to end of last silent
    % period (=beginning of last, incompletely recorded and therefore
    % excluded burst)
    stats.cIntv=[aEtsl(1,etslc.tsCol) sum(iEtsl(end,[etslc.tsCol etslc.durCol]))];
  else
    % equivalent inverse
    stats.cIntv=[iEtsl(1,etslc.tsCol) sum(aEtsl(end,[etslc.tsCol etslc.durCol]))];
  end
else
  stats.cIntv=[nan nan];
end

if nA
  % bursts within stats.cIntv:
  bIx=find(aEtsl(:,etslc.tsCol)>=stats.cIntv(1) & ...
    aEtsl(:,etslc.tsCol)<stats.cIntv(2));
  % - burst rate
  stats.burstRate=numel(bIx)/(diff(stats.cIntv)/1000);
  % - the inverse, burst period (the time interval from one burst
  % beginning to the next burst beginning):
  if nA>1
    df=diff(aEtsl(:,etslc.tsCol));
    stats.mnBurstPeriod=mean(df);
    stats.stdBurstPeriod=std(df);
    stats.geomnBurstPeriod=geomean(df);
    stats.mdBurstPeriod=median(df);
  end
  % - burst length median etc.
  tmp=prctile(aEtsl(:,etslc.durCol),[0 5 25 50 75 95 100]);
  stats.minBurstLen=tmp(1);
  stats.lqBurstLen=tmp(3);
  stats.mdBurstLen=tmp(4);
  stats.uqBurstLen=tmp(5);  
  stats.maxBurstLen=tmp(7);
  % asymmetry of burst length:
  % (95th percentile - 2*median + 5th percentile)/(95th percentile - 5th percentile)
  stats.asBurstLen=(tmp(6)-2*tmp(4)+tmp(2))/(tmp(6)-tmp(2));
  % - burst length (arithmetic mean & std)  
  stats.mnBurstLen=mean(aEtsl(:,etslc.durCol));
  stats.stdBurstLen=std(aEtsl(:,etslc.durCol));
  % - burst length (geometric mean)  
  stats.geomnBurstLen=geomean(aEtsl(:,etslc.durCol));
  % - relative time spent in active state
  stats.relTimeInBurst=sum(aEtsl(bIx,etslc.durCol))/diff(stats.cIntv);
end
if nI
  % if list of inactive events contains at least two entries, compute burst
  % period length based on this list as well and average with values above
  if nI>1
    df=diff(iEtsl(:,etslc.tsCol));
    stats.mnBurstPeriod=nanmean([stats.mnBurstPeriod  mean(df)]);
    stats.stdBurstPeriod=nanmean([stats.stdBurstPeriod  std(df)]);
    stats.geomnBurstPeriod=nanmean([stats.geomnBurstPeriod  geomean(df)]);
    stats.mdBurstPeriod=nanmean([stats.mdBurstPeriod  median(df)]);
  end
  % - silent period length median etc.
  tmp=prctile(iEtsl(:,etslc.durCol),[0 5 25 50 75 95 100]);
  stats.minSilentPerLen=tmp(1);
  stats.lqSilentPerLen=tmp(3);
  stats.mdSilentPerLen=tmp(4);
  stats.uqSilentPerLen=tmp(5);  
  stats.maxSilentPerLen=tmp(7);  
  stats.asSilentPerLen=(tmp(6)-2*tmp(4)+tmp(2))/(tmp(6)-tmp(2));
  stats.mnSilentPerLen=mean(iEtsl(:,etslc.durCol));
  stats.stdSilentPerLen=std(iEtsl(:,etslc.durCol));
  stats.geomnSilentPerLen=geomean(iEtsl(:,etslc.durCol));  
end


