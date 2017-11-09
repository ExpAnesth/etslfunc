function [h,binArr,lag,cbi]=tslcc(xtsl,ytsl,varargin)
% ** function [h,binArr,lag,cbi]=tslcc(xtsl,ytsl,varargin)
%    computes cross-correlations between two time stamp lists (=peri-event
%    time histograms). Most of the nomenclature used here is lent from
%    xcorr.m, a function provided by The Mathworks producing estimates of
%    cross correlations for sampled data series such as field potential
%    traces. The analogy between the two functions is best understood if
%    the time stamp lists (the input variables here) are imagined as
%    sampled data series with needles of area 1 where a spike occurs and 0
%    where this is not the case. The bin width 'binw' then corresponds to
%    the sampling interval in sampled series, and the number of bins of
%    this length that would be needed to cover the longer of the two tsl
%    corresponds to the number of data points of a sampled series. Note
%    that zero-padding of the shorter tsl is implied.
%
%          ** time unit is ms for all variables **
%
%                    >>> INPUT VARIABLES >>>
%
% NAME          TYPE/DEFAULT         DESCRIPTION
%
% xtsl,ytsl     column arrays        time stamp lists
% lag           scalar, 199.5        lags
% binw          scalar, 1            the bin width in ms
% norm          char array, 'none'   normalization of histogram:
%                           'biased' - division by nMaxLagBins, the number 
%                             of bins needed to cover the longer of the two
%                             tsl (see above)
%                           'unbiased' - each bin is divided by 
%                             nMaxLagBins-abs(lag)
%                           'coeff' - division by the average of the number
%                             of spikes on both channels
% cb            char array, 'leave'   determines value to which central bin will be set
%                                      'leave' or 'none' = leave as it is
%                                      'del' or 'nan' or 'kill' = nan
%                                      'zero' or '0' = 0
%                                      this feature is useful for autocorr plots
%                                      where central bin is mostly not wanted
% mr            scalar, 0            determines whether mean of all bins shall be
%                                    removed; 0=no, 1=yes; useful for computation of spectra
% maxArrBytes   scalar, 1e9          max size of largest temporary matrix 
%                                    tslcc.m will set up for internal
%                                    computations in bytes; set maximally
%                                    to somewhat less than the RAM
%                                    available to Matlab (GPU RAM in case
%                                    doComputeGP==true!)
% doGPUCompute  logical, true       if true, computations will run on GPU
%                                   (may be faster)
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT          DESCRIPTION
%
% h             column array          the histogram
% binArr        column array          the CENTERS of the bins of the histogram
%                                     ** NOTE: center of center bin = 0
% lag           scalar                the lag realizable given bin width
%                                     ** realize that with e.g. binw=1, the realizable lags are
%                                     for example -10.5 and 10.5, but never integer values like -10 and 10
%                                     because the central bin is centered at 0
% cbi           scalar                index to central bin

% hh 04/03
% changes 07/03:
% - switched from hist to histc (which is faster)
% - autocorrs treated differently internally (no more symmetry problems)
% changes Feb 2015:
% - now computing xtls-ytls
% - numerous rearrangements of code and bug fixes
% changes Nov 2017:
% - replaced Histc by histcounts
% - computation of difference matrix via automatic expansion of arrays
% - implemented possibility to compute on GPU

% defaults
lag=199.5;
binw=1;
cb='leave';
norm='none';
mr=0;
doGPUCompute=true;
maxArrBytes=1e9;
% modify according to input
pvpmod(varargin);

% set up bin array
finalBinArr=[fliplr(-1*(binw:binw:lag)) 0:binw:lag]';
nFinalBins=length(finalBinArr);
finalCbi=(nFinalBins+1)/2;

% check tsl
nx=length(xtsl);
ny=length(ytsl);

if ~nx || ~ny
  warning('at least one tsl is empty - output will consist of zeroes');
  h=zeros(size(finalBinArr));
  binArr=finalBinArr;
  lag=binArr(end);
  cbi=finalCbi;
else
  if isequal(xtsl,ytsl)
    auto=1;
  else
    auto=0;
  end
  
  % determine start and stop points in time common to both tsl
  startT=min([xtsl(1) ytsl(1)]);
  stopT=max([xtsl(end) ytsl(end)]);
  maxlag=stopT-startT;
  
  % cut down lag if reasonable
  if lag>maxlag
    lag=floor((maxlag-.5*binw)/binw)*binw+1.5*binw;
    lagBinBorders=[-.5*binw:binw:lag];
    nMaxLagBins=length(lagBinBorders);
  else
    % this is the number of bins of width 'binw' needed to cover
    % the time interval which includes both tsl completely (1st spike
    % at left border of 1st interval)
    nMaxLagBins=cont2discrete(maxlag,binw);
    lagBinBorders=[-.5*binw:binw:lag+.5*binw];
  end
  
  if auto
    binArr=lagBinBorders';
    cbi=length(binArr)-1;
  else
    binArr=[-fliplr(lagBinBorders(3:end)) lagBinBorders]';
    % binArr=[fliplr(-.5*binw:-binw:-lag) .5*binw:binw:lag]';
    cbi=length(binArr)/2;
  end
  
  % 'true' lag
  if auto
    lag=binArr(end);
  else
    lag=unique(abs(binArr([1 end])));
  end
  if length(lag)>1
    error('computation of bins & lags messed up');
  end
  % initialize histogram container
  h=zeros(numel(binArr)-1,1);
  
  if doGPUCompute
    xtsl=gpuArray(xtsl);
    ytsl=gpuArray(ytsl);
    binArr=gpuArray(binArr);
  end
  
  % if two double matrices of size ny by nx are too large for the assumed
  % memory accessible to Matlab divide xtsl into chunks
  xilen=min(nx,round(maxArrBytes/(2*ny*8)));
  [intrvls,intrvls_pts]=mkintrvls([0 nx],'ilen',xilen,'olap',0,'border','include');
  if size(intrvls,1)>1
    disp(['dividing calculation of difference matrix in ' int2str(size(intrvls,1)) ' chunks']);
  end
  
  for ccnt=1:size(intrvls_pts,1)
    xidx=intrvls_pts(ccnt,1):intrvls_pts(ccnt,2);
    subnx=length(xidx);
    % a matrix sized ny * nx holding all differences xtsl-ytsl (** note new
    % R2016x automatic expansion arithmetic!)
    tsDiff=xtsl(xidx)'-ytsl;
    % all values ytsl-xtsl not within interval won't appear in the
    % histograms because histcounts ignores values outside the specified
    % bins 
    curH=histcounts(tsDiff(:),binArr);
    h=h+curH';
  end
    
  % get rid of last bin
  binArr(end)=[];
  
  % shift bins to the right by binw/2 so binArr represents the center of intervals
  binArr=binArr+.5*binw;
  
  if auto
    binArr=[-1.0*flipud(binArr(2:end)); binArr];
    h=[flipud(h(2:end)); h];
  end
  if binArr(cbi)
    error('internal: identification of central bin');
  end
  % if lag had been set to maxlag the following must be true here:
  % length(binArr)==2*length(nMaxlagBinBorders)-1==2*cbi-1
  % let's indirectly check this by precomputing the unbiased normalization factor
  ubFac=([1:cbi [cbi-1:-1:1]]+nMaxLagBins-cbi)';
  if any(ubFac<1)
    error('internal: computation of bins');
  end
  
  normFac=[];
  switch norm
    case 'biased'
      normFac=nMaxLagBins;
    case 'unbiased'
      % attention: if lag  > maxlag normFac will contain zeroes & negative factors!
      normFac=ubFac;
    case 'coeff'
      normFac=(nx+ny)/2;
    case 'none'
    otherwise
      warning('normalization option not valid - no normalization applied');
  end
  
  if isempty(normFac)
  elseif length(normFac)==1
    h=h/normFac;
  else
    h=h./normFac;
  end
  
  % taking care of central bin
  switch cb
    case {'del','nan','kill'}
      h(cbi)=nan;
    case {'zero','0'}
      h(cbi)=0;
    case {'leave','none'}
      % do nothing
    otherwise
      warning('option for treatment of central bin not valid - no changes applied');
  end
  
  if mr
    h=h-mean(h([1:cbi-1 cbi+1:end]));
  end
  
  % pad with zeros if necessary
  if length(h)<nFinalBins
    binArr=finalBinArr;
    h=[zeros(finalCbi-cbi,1); h; zeros(finalCbi-cbi,1)];
  end
  
  h=gather(h);
  binArr=gather(binArr);
end
