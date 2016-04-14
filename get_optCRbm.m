%% get_optCRbm(bm, model)
%
% compute biomass according to input and be sensitive to growth condition
% (auto - no organic carbon source (no acetate), hetero - no anorganic 
% carbon source (no photon uptake))
% 
% IN: 
%       bm    - which biomass function? 
%           1 - auto
%           2 - hetero
%           3 - mixo
%       model - a COBRA model structure
%
% OUT:
%       obj  - optimized FBA value
%       xVec - opt. flux distribution
% 
% @Sascha Schäuble

function [obj, xVec] = get_optCRbm(bm, model)

switch bm
    case 1 % autotrophic
        % no acetate uptake, since autotrophic growth
        model.lb(28) = 0;
        bmID = 62;
    case 2 % mixotrophic
        bmID = 63;
    case 3 % heterotrophic, no light uptake, since heterotrophic growth
        model.lb(1:12) = 0;
        model.ub(1:12) = 0;
        bmID = 64;
end

model = changeCbOpt(model, bmID);

res   = optimizeCbModel(model);
obj   = res.f;
xVec  = res.x;