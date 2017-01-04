% runAll
%
% header script 
% executes runAnalysis and all different types of simulations 
% based on iRC1080 
%
%

%
% initialize COBRA toolbox if not done so already by calling
% initCobraToolbox;
%

modulePath = './data/ALGdL_modules.mat'; % change to whatever your path is

% KO sims of single KO and module KO
step4 = runAnalysis(4, [0.1 0.3 0.5]);

% KD sims of single KD and module KD
% step5 = runAnalysis(5, [0.1 0.3 0.5], [0.0625 0.125 0.25 0.5]);

%
% compute and visualize ratios of mod TAG & BM / wt TAG & BM
% step4_ratio = visAnalysis(step4, 1,'step4');
% step4mod_ratio = visAnalysis(step4mod, 1,'step4ppimod');

% step5_ratio = visAnalysis(step5, 2,'step5');
%

% filter results for at least 5% change in TAG and/or BM production
% filtSingles_step4     = filterRes(iRC1080, step4, 4, 1, 'singleKO');
filtSingles_step4mod     = filterRes(iRC1080, step4mod, 4, 1, 'singleKO_ppiMod');
% filtSingles_step5bm10 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM10', 'lb10');
% filtSingles_step5bm30 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM30', 'lb30');
% filtSingles_step5bm50 = filterRes(iRC1080, step5, 5, 1, 'singleKD_BM50', 'lb50');
% %
% 
% filtModules_step4     = filterRes(iRC1080, step4, 4, 2, 'moduleKO'     ,     [], modulePath);
% filtModules_step5bm10 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM10', 'lb10', modulePath);
% filtModules_step5bm30 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM30', 'lb30', modulePath);
% filtModules_step5bm50 = filterRes(iRC1080, step5, 5, 2, 'moduleKD_BM50', 'lb50', modulePath);
%}