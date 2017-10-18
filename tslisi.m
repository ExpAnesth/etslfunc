function [h,x]=tslisi(tsl,varargin)
% ** function [h,x]=tslisi(tsl,varargin)
% computes an inter-spike-interval distribution
%                         >>> INPUT VARIABLES >>>
%
% NAME             TYPE/DEFAULT            DESCRIPTION
% tsl              1d-array                time stamp list 
% bins             scalar or 1d-array,100  number of bins or bin CENTERS of the histogram
% ord              scalar, 1               order of the isi (2=time to second next event, etc.)
%                     
%                         <<< OUTPUT VARIABLES <<<
%
% NAME             TYPE/DEFAULT            DESCRIPTION
% x                column array            the histogram bins
% h                column array            the histogram (bin counts)
%
% 
%  ------- non-standard matlab functions used:
%  makerow

% default values
verbose=0;
bins=100;
ord=1;

pvpmod(varargin);
if verbose, disp(['**** ' mfilename ':']); end;


if ord>1
  % make tsl a row vector
  tsl=makerow(tsl);
  % cut down length to integer multiple of ord
  tsl=tsl(1:floor(length(tsl)/ord)*ord);
  % reshape such that rows contain every ord-th ts 
  tslLen=length(tsl);
  tsl=reshape(tsl,ord,tslLen/ord);
  % do it
  isin=diff(tsl,1,2);
  % now re-reshape
  isin=reshape(isin,1,tslLen-ord);
else
  isin=diff(tsl);
end
% make histogram
if length(bins)>1
  % add bins at either end and delete them afterwards
  binw=min(diff(bins));
  bins=[bins(1)-binw; bins(:); bins(end)+binw];
  [h,x]=hist(isin,bins);
  h=h(2:end-1);
  x=x(2:end-1);
else
  [h,x]=hist(isin,bins);
end
h=h(:);
x=x(:);
