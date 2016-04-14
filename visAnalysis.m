%% visAnalysis
%
% visualize runAnalysis results
%
% IN: rA_res  - runAnalysis result
%     stType  - 0 - only fixed BMmod and FLUXmod (both not given)
%             - 1 - varied BMmod, no FLUXmod
%             - 2 - varied/fixed BMmod, varied FLUXmod
%   OPTIONAL
%     figName - if given save pics as pdf, otherwise no pics are presented
%
% OUT:
%     res     - gather single and module KO analysis 
%
% e.g.: step4vis = visAnalysis(step4_lb001_035,1) - step4, no FLUXmod, but 
% BMmod and no figure print as filename not provided
%

function [res] = visAnalysis(rA_res, stType, figName)

% check if pics should be save or presented as matlab panels
if exist('figName','var')
    switch stType
        case 0
            res.singleKO = ...
                visResults(1, rA_res, ...
                [figName, '_singleMOD']);
            res.singleKO = ...
                visResults(2, rA_res, ...
                [figName, '_moduleMOD']);
        case 1
            BMmods = fieldnames(rA_res);
            for i = 1:length(BMmods)
                res.(char(BMmods(i))).singleKO = ...
                    visResults(1, rA_res.(char(BMmods(i))), ...
                    [figName, '_singleMOD_', BMmods{i}]);
                res.(char(BMmods(i))).moduleKO = ...
                    visResults(2, rA_res.(char(BMmods(i))), ...
                    [figName, '_moduleMOD_', BMmods{i}]);
            end
        case 2
            BMmods = fieldnames(rA_res);
            for i = 1:length(BMmods)
                FLUXmods = fieldnames(rA_res.(char(BMmods(i))));
                for j = 1:length(FLUXmods)   
                    res.(char(BMmods(i))).(char(FLUXmods(j))).singleKO = ...
                        visResults(1, rA_res.(char(BMmods(i))).(char(FLUXmods(j))), ...
                        [figName, '_singleMOD_', BMmods{i}, '_', FLUXmods{j}]);
                    res.(char(BMmods(i))).(char(FLUXmods(j))).moduleKO = ...
                        visResults(2, rA_res.(char(BMmods(i))).(char(FLUXmods(j))), ...
                        [figName, '_moduleMOD_', BMmods{i}, '_', FLUXmods{j}]);
                end
            end
    end
else % no pics
    switch stType
        case 0
            res.singleKO = ...
                visResults(1, rA_res);
            res.singleKO = ...
                visResults(2, rA_res);
        case 1
            BMmods = fieldnames(rA_res);
            for i = 1:length(BMmods)
                res.(char(BMmods(i))).singleKO = ...
                    visResults(1, rA_res.(char(BMmods(i))));
                res.(char(BMmods(i))).moduleKO = ...
                    visResults(2, rA_res.(char(BMmods(i))));
            end
        case 2
            BMmods = fieldnames(rA_res);
            for i = 1:length(BMmods)
                FLUXmods = fieldnames(rA_res.(char(BMmods(i))));
                for j = 1:length(FLUXmods)   
                    res.(char(BMmods(i))).(char(FLUXmods(j))).singleKO = ...
                        visResults(1, rA_res.(char(BMmods(i))).(char(FLUXmods(j))));
                    res.(char(BMmods(i))).(char(FLUXmods(j))).moduleKO = ...
                        visResults(2, rA_res.(char(BMmods(i))).(char(FLUXmods(j))));
                end
            end
    end
end