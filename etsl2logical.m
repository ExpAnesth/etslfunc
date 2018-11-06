function logicalIx=etsl2logical(numPts,etsl)
% ** logicalIx=etsl2logical(numPts,etsl)
% generates from extended time stamp list etsl logical array logicalIx of
% length numPts (sampling) points in which true values represent specific
% periods (e.g. bursts) in the underlying 1D time series (which has length
% numPts).
% ** Note: the time basis of both input arguments must be SAMPLING POINTS!

etslconst;
assert(~any(rem(etsl(:,etslc.tsCol),1)),...
  'time basis of extended time stamp list must be sampling points');
logicalIx=false(numPts,1);
% start and stop times of active periods in points
cumTimes=cumsum(etsl(:,[etslc.tsCol etslc.durCol]),2)-[0 1];
% indexes to data points belonging to periods, omitting incompletely
% detected periods at beginning and end
for k=1:size(cumTimes,1)
  logicalIx(cumTimes(k,1):cumTimes(k,2))=true;
end
