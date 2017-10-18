% this is a (real-life) example of the input structure needed for running
% utility_effectSize_peth

% mat file name (without extension) containing PETH data
% - perievdeal example
ds.dFn='perievdeal_Blocker_Dia';
% % - tslbatch example
% ds.dFn='Blocker_Dia';
% directory in which data reside:
% - perievdeal example
ds.dDir='e:\_data\otc_ctx\ACh\Diazepam\Blocker_Dia\periev_figs\'; 
% % - tslbatch example
% ds.dDir='e:\_data\otc_ctx\ACh\Diazepam\Blocker_Dia\Figs\';
% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
% NOTE: for a comparison of two different data sets, it is possible to
% specify two different files and directories. In this case, the file
% names/dirs must be given in a cell array, and a number of assumptions are
% made and/or rules apply, among them the following (be aware that rigorous
% checks are sparse):
% - bin sizes of the PETHs must be identical
% - ds.curIpVal below should contain two identical values (these are the
% independent levels (in most cases drug concentrations) to be compared
% among the data sets)
% Here is an example:
% ds.dFn={'perievdeal_Blocker_Dia','perievdeal_ACh_Dia'};
% ds.dDir={'d:\_data\otc_ctx\ACh\Diazepam\Blocker_Dia\periev_figs\',...
%   'd:\_data\otc_ctx\ACh\Diazepam\ACh_Dia\periev_figs\'};
% *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

% directory into which to save results
ds.exportDir='d:\hh\projects\ctx_ACh\rawFig\';
% name of exported figure
ds.exportFn=ds.dFn;
% the 'events' defining time 0 of the histograms ('start' or 'end')
ds.typeEvent='start';
% interval for peth
ds.intv=[-50 600];
% the two independent levels to be compared (in most cases, these are
% concentrations)
ds.curIpVal=[1 2];
% index to experimental IDs to omit:
% first group
ds.killer(1).ix=[];
% second group
ds.killer(2).ix=[];
% recording sites to use (set to 'av' to average rec sites within
% experiment irrespective of which and how many there are per experiment)
% ** note: if processing older data sets the code will regard the values of
% column 'undef1' as the recording site
ds.curRecSite=[5];
% ds.curRecSite='av';
% upper y limit of peth (nan for automatic scaling
ds.pethYLim=nan;
% type of effect size measure to be computed (see help of effectsz.m)
ds.esType='auroc';
% the value corresponding to a null effect of the effect size measure in
% question
ds.NullEffectVal=0.5;
% IMPORTANTISSIMO: the function to be used for averaging peths from
% different experiments (specify a function handle, @nanmedian or @nanmean)
ds.avFuncH=@nanmedian;
% y limits for effects size plot ([nan nan]) for automatic scaling
ds.esYLim=[0 1];
% colors for the two different data sets
ds.pCol=[0 0 0;1 0 0];
% print?
ds.printas=[]; % no
% ds.printas='-djpeg95'; % q&d
% ds.printas='-dpsc2'; % export into corel etc
% *** are these paired data or not?
ds.pair=1;
% shall results be plotted?
ds.doPlot=true;
% shall individual experiments be plotted?
ds.indvPlot=false;
% style in which to plot effect size measure & error bars (=confidence intervals) 
% (any of 'errorbar','boundedline','area')
ds.mesPlotStyle='boundedline';

% *** a few graphics settings
labelscale('scaleFac',.8,'fontSz',7,'lineW',.5,'markSz',2);

% *** call utility 
utility_effectSize_peth(ds);

% (alternatively, with output argument(s) - if you don't need them, don't
% specify them)
% [stats]=utility_effectSize_peth(ds);