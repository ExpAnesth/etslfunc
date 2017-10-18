function [cutout,isCutout,netsl]=etsl2exc(d,si,etsl,varargin)
% ** function [cutout,isCutout,netsl]=etsl2exc(d,si,etsl,varargin)
%    extracts variable-length cutouts from continuous raw data d.
%    Cutouts are listed in extended time stamp list etsl.
%
%                    >>> INPUT VARIABLES >>>
%
% NAME           TYPE/DEFAULT              DESCRIPTION
% d              1d                        sampled data, time runs along column
% etsl           2d arr                    extended time stamp list  
% si             a) scalar                 sampling interval in us
%                b) string 'idx'           means that ALL timing information 
%                                            is specified in points (as opposed to ms)
%                                            and will also be put out as such
% win            2-element arr, [-1 1]     cutout window relative to the start and stop
%                                            time stamp defining an event (ms or pts)
%
%                    <<< OUTPUT VARIABLES <<<
%
% NAME           TYPE/DEFAULT               DESCRIPTION
% cutout         cell array                 cutouts
% isCutout       1d arr                     an index into extended time stamp list:
%                                             1=cutout could be produced for time stamp
%                                             0=cutout could NOT be produced (was 
%                                             on the border of or outside of sweep)
% netsl          1d arr                     extended time stamp list devoid of events not 
%                                             respecting borders


% improvements:
% - have option to place data in array with fixed length because in many
% applications that's what will have to happen anyways

win=[-1 1];
verbose=0;
pvpmod(varargin);

if verbose, disp(['**** ' mfilename ':']); end
etslconst;

[n1 n2]=size(d);
if n2>1
  error('d must be a single-column array')
end
nTs=size(etsl,1);

if ischar(si)
  if ~strcmpi(si,'idx')
    error('input variable si must be either a scalar or string ''idx''');
  end
else
  % ** convert to ticks (easier to calculate with) **
  win=cont2discrete(win,si*.001,'intv',1);
  etsl(:,[etslc.tsCol etslc.durCol])=cont2discrete(etsl(:,[etslc.tsCol etslc.durCol]),si*.001);
end
preTrigPts=win(1);
postTrigPts=win(2);

cutout={};
netsl=[];
isCutout=[];

if nTs
  etslconst;
  % first locate events fully contained in sweep by inspecting border points
  ix=repmat([preTrigPts postTrigPts],nTs,1)+cumsum(etsl(:,[etslc.tsCol etslc.durCol]),2);
  isCutout=((ix(:,1)>=1) & (ix(:,2)<=n1));
  % delete others from list
  etsl(~isCutout,:)=[];
  ix(~isCutout,:)=[];
  %  NEW number of time stamps..
  nTs=size(etsl,1);
  if nTs
    %  ..and tsl
    netsl=etsl;
    % initialize cutout
    cutout=cell(1,nTs);
    for g=1:nTs
      cutout{g}=d(ix(g,1):ix(g,2),1);
    end
  else
    warning('no single event in etsl suitable for cutout - putting out empty results variables');
  end    
else
  warning('etsl is empty - putting out empty results variables');
end


if ~ischar(si)
  % no need to check whether si==idx
  etsl(:,[etslc.tsCol etslc.durCol])=discrete2cont(etsl(:,[etslc.tsCol etslc.durCol]),si*.001);
end

