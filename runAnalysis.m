%% runAnalysis
%
% runs analysePapinModel with different setups
%
% IN
%       step <int>    - which step of analysePapinModel (mostly 4 or 5)
%       BMmod <vector>   - What minimal biomass flux is required?
%                        - [0 1]
%       fluxMod <vector> - How should affected fluxes be modified?
%	modulesPath	 - path to compiled modules as input for analysis
%
% OPTIONAL 
%       model            - COBRA model
%
% OUT
%       res <struct>     - structurally organized FBA analysis performed by
%       analysePapinModel function according to manipulation resulting from 
%       input data 
%
% e.g.
% step4 = runAnalysis(4, [0.1 0.3 0.5]);
% step5 = runAnalysis(5, 0.3, [0.0625 0.125 0.25 0.5]);

function [res] = runAnalysis(step, BMmod, FLUXmod, modulesPath, model)

if ~exist('model','var')
    model = [];
end

% set PRISM reactions to use:
prism = [4 5];

switch step
    case 2 % single KO is requested
        res = analysePapinModel(2, prism);
    case {4, 5, 6, 7} % single, double or modules KO/KD is requested
        
        % create string structure names according to input data
        if exist('BMmod', 'var') && ~isempty(BMmod)
            bmMods = cell(length(BMmod), 1);
            for i = 1:length(BMmod)
                bmMods{i} = ['lb' num2str(floor(BMmod(i)*100))];
            end
        end

        if exist('FLUXmod', 'var')
            fluxMods = cell(length(FLUXmod), 1);
            for i = 1:length(FLUXmod)
                fluxMods{i} = strcat('fluxMod', num2str(floor(FLUXmod(i)*100)));
            end
        end


        % produce results
        for i = 1:length(BMmod)
            if exist('FLUXmod', 'var') && ~isempty(FLUXmod)
                for j = 1:length(FLUXmod)
                    % analysis
                    res.(char(bmMods(i))).(char(fluxMods(j))) = ...
                        analysePapinModel(step, prism, BMmod(i), FLUXmod(j), model, modulesPath);
                end
            else
                % analysis
                res.(char(bmMods(i))) = ...
                        analysePapinModel(step, prism, BMmod(i), [], model, modulesPath);   
            end
        end        
end