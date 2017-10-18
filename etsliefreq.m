function [t,ief,varargout]=etsliefreq(cetsl,si,varargin)
% ** function [t,ief,varargout]=etsliefreq(cetsl,si,varargin)
%    computes the 'instantaneous event frequency' (ief) from a list of
%    events with nonzero duration. cetsl is a list of events (possibly of
%    different types) which have a start time and duration and which may or
%    may not overlap with each other. The ief specifies the number of
%    events occurring (or being 'on') at each point in time (think of the
%    events as conversations between pairs of people in a crowd: you want
%    to know, at each point in time, how many conversations are going on).
%    Be aware that the time resolution must be finer than the shortest
%    event, otherwise these events will be either over- or underrepresented
%    (in time). 
% 
%                    >>> INPUT VARIABLES >>>
%
% NAME        TYPE/DEFAULT        DESCRIPTION
% cetsl       2d array            the list of events. It must adhere to the
%                                 conventions for etsl with the exception that
%                                 events may be overlapping in time
% ****  see etslconst.m for information on variables of type 'etsl' **** 
% si          scalar              desired time resolution of the result 
%                                 (same unit as in cetsl must be specified)
% syncT       scalar,[]           synchronizing time point: if specified, the 
%                                 time raster will be aligned to this point
% emptyV      scalar, 0           the default value to insert in bins of the 
%                                 NORMALIZED ief which are not covered by any
%                                 event. Nan would be another ovious choice
%                                 which, however, is no good with area plots
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME          TYPE/DEFAULT    DESCRIPTION
% ief           2d array        'instantaneous event frequency': each event type 
%                               will be represented in one column
%                               ** columns are sorted according to tag value (lower
%                                  ones first) **
% varargout{1}  2d array        normalized ief: all components add up to one;
%                               time bins with no event are NaN
% 
%    For plots of results, use
%            area(t,ief);
%    or
%            plot(t,cumsum(ief,2));
%

% improvements:
% This function is very inefficient with sparse events spread out in time,
% since it expands a time stamp list into a regularly spaced array covering
% the whole time interval. As its only use is for plots, an improved
% version would (alternatively) spit out the results in a format suitable
% for the line command (just x and y values of the corners)

verbose=0;
syncT=[];
emptyV=0;
pvpmod(varargin);

% ------ preparatory things
etslconst;
if verbose, disp(['** ' mfilename]); end

esz=size(cetsl,1);
if ~esz, error('input variable ''cetsl'' is empty'); end
% the column holding event START times
startcol=etslc.tsCol;
% temporary new column: event STOP times
stopcol=etslc.tagCol+1;
durcol=etslc.durCol;
tagcol=etslc.tagCol;
% --- prepare cetsl, make new variables
% sort according to start times
cetsl=sortrows(cetsl,startcol);
% create temporary new column  (start + duration)
cetsl(:,stopcol)=cetsl(:,startcol)+cetsl(:,durcol);
% time raster: these are the bin CENTERS
% @to do: end points
if isempty(syncT)
  t=(cetsl(1,startcol)+si/2:si:cetsl(end,stopcol)+si/2)';
else
  t=[fliplr(syncT-si:-si:cetsl(1,startcol)-si/2)  syncT:si:cetsl(end,stopcol)+si/2]';
end
nt=length(t);
% unique event types
uevType=unique(cetsl(:,tagcol));
nuevType=length(uevType);
% preallocate ief
ief=zeros(nt,nuevType);

for eix=1:size(cetsl,1)
  % column index into ief
  cix=find(uevType==cetsl(eix,tagcol));
  % 'time' index (into bins of t): 
  tix=[t>=cetsl(eix,startcol) & t<cetsl(eix,stopcol)];
  ief(tix,cix)=ief(tix,cix)+1;
end

if nargout>2
  emptyix=~any(ief,2);
  % fill empty bins with ones to prevent division by zero
  ief(emptyix,:)=1;
  varargout{1}=ief./repmat(sum(ief,2),1,nuevType);
  % set to zero again..
  ief(emptyix,:)=0;  
  % and set bins in normalized array to desired value
  varargout{1}(emptyix,:)=emptyV;  
end