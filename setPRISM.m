%% select PRISM reaction
%
% modifies given model to selected PRISM reactions as given in prismIDs
% assumes PRISM reaction IDs are [1:12] 
%
% take care whether only ub or also lb should be updated 
%
function res = setPRISM(prismIDs, model)
tmpLB = model.lb(prismIDs);
tmpUB = model.ub(prismIDs);

model.lb(1:12) = 0;
model.ub(1:12) = 0;

% model.lb(prismIDs) = tmpLB; % if executed, PRISM flux is fixed
model.ub(prismIDs) = tmpUB;

res = model;