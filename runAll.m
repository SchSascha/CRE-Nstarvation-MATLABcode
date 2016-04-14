% runAll
%
% script that executes runAnalysis and all different sorts of different 
% experiments via analysePapinModel
% things
%

% initialize COBRA toolbox
% initCobraToolbox;
%

% KO sims of single KO and module KO
% step4 = runAnalysis(4, [0.1 0.3 0.5]);
step4mod = runAnalysis(4, [0.1 0.3 0.5], [], './data/ALGdL_modules.mat', iRC1080mod);

% KD sims of single KD and module KD
% step5 = runAnalysis(5, [0.1 0.3 0.5], [0.0625 0.125 0.25 0.5]);

%
% compute ratios of mod TAG & BM / wt TAG & BM
% step4_ratio = visAnalysis(step4, 1,'step4_140609');
% step4mod_ratio = visAnalysis(step4mod, 1,'step4ppimod_140609');

% schr5_ratio = visAnalysis(step5, 2,'schr5_140609');
%

% filter results for at least 5% change in TAG and/or BM production
% filtSingles_step4     = filterRes(iRC1080, step4, 4, 1, 'singleKO');
filtSingles_step4mod     = filterRes(iRC1080, step4mod, 4, 1, 'singleKO_ppiMod');
% filtSingles_step5bm10 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM10', 'lb10');
% filtSingles_step5bm30 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM30', 'lb30');
% filtSingles_step5bm50 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM50', 'lb50');
% %
% 
% filtModules_step4     = filterRes(iRC1080, step4, 4, 2, 'moduleKO'     ,     [], './data/ALGdL_modules.mat');
% filtModules_step5bm10 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM10', 'lb10', './data/ALGdL_modules.mat');
% filtModules_step5bm30 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM30', 'lb30', './data/ALGdL_modules.mat');
% filtModules_step5bm50 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM50', 'lb50', './data/ALGdL_modules.mat');
%}
