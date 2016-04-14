%% protocol for analysis of Papin network iRC1080, Chang et al 2011
%
% IN:
%       step  - choose which analysis to take (see description in switch/case below)
%       affPrism - <vector> Which prism reactions should be used? 
%   OPTIONAL
%       bmdrop   - minimal lower bound for biomass (bm) reaction
%                - default: 0.3
%       kdown    - if knock down instead of knock out is investigated,
%       provide degree of knock down
%       model    - alternative (modified) Papin model
%	modulesPath - path to compiled modules as input for analysis

function res = analysePapinModel(step, affPrism, bmdrop, kdown, model, modulesPath)
global cr

% assign bmdrop if not given
if ~exist('bmdrop', 'var')
    bmdrop = 0.3;
end

% load WT model anew:
if exist('model','var') & ~isempty(model)
    cr = model;
else
    cr = get_model();
end

    
switch step
    % Step 1
    % =========
    % compute WT biomass for all three types of biomass for all light sources
    %
    case  1
        res = get_BMopt();
    
    % Step 2
    % =========
    % simulate single KO of all genes and provide BM
    case 2
        tic;
        res = singleKO(affPrism);
        toc;

    % Step 3
    % =========
    % simulate KO of all modules and provide BM
    case 3
        tic;
        res = moduleKO(affPrism, modulesPath);
        toc;
        
    % Step 4
    % =========
    % extract TAG production and simulate gene and module KO
    % 'bmdrop' enables maximisation of TAG production
    case 4
        tic;
        res.a = maxTAGprod_WT(affPrism, bmdrop);
        res.b = maxTAGprod_singleMOD(affPrism, bmdrop);
        res.c = maxTAGprod_moduleMOD(affPrism, bmdrop);
        toc;
    
    % Step 5
    % =========
    % extract TAG production and simulate gene and module KD
    % same as step four, except for 'kdown' - knock down instead of knock out
    %
    % IN: bmdrop - min. lb biomass flux
    %     kdown  - what fraction of wt flux for knock down?
    case 5
        tic;
        res.a = maxTAGprod_WT(affPrism, bmdrop);
        res.b = maxTAGprod_singleMOD(affPrism, bmdrop, kdown);
        res.c = maxTAGprod_moduleMOD(affPrism, bmdrop, kdown);
        toc;           
    
        
end



% =========================================================================
% =========================================================================


% max TAGprod will be comprised by three functions:
% =================================================
%
% maxTAGprod_WT gives back WT biomass and TAG production FBA
% maxTAGprod_singleMOD
% maxTAGprod_moduleMOD


%% maxTAGprod_WT
%
% extract TAG production and compute optBM and optTAG independently
% basis: WT-model (no KO or KD)
% bmDrop * lb(bm) necessary to calculate opTAG while still having bm feasible
%
% IN: prism  - which prism reactions should be used as light source?
%     bmDrop - [0 1] what fraction of WT biomass should be kept?
%
% OUT: res.
%       BM.optBM  - optimized WT bm value
%       BM.TAG    - associated Tag flux, when BM is maximized
%       BM.x      - associated flux distribution
%       TAG.optBM - optimized WT TAG value
%       TAG.TAG   - associated BM flux, when TAG is maximized
%       TAG.x     - associated flux distribution
%       model     - used model
%
%
function res = maxTAGprod_WT(prism, bmDrop)
global cr

% 0.
%
% modify base network
% extract TAG production and formulate as extra function
tmp = modModel_extraTAG(cr);

% set PRISM 
tmp = setPRISM(prism, tmp);

% restrict uptake:
% ----------------
% no starch
tmp.lb(27) = 0;
% no no3
tmp.lb(17) = 0;
% no nh4
% tmp.lb(16) = 0;

res.model = tmp;

% 1. get opt. WT growth yield & TAG flux alongside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:3
    [res.BM.optBM(i), foo] = get_optCRbm(i, tmp);
    res.BM.TAG(i)   = foo(end);
    res.BM.x(:,i)       = foo;
end
    
% 1.1 restrict minBM flux to bmDrop of WT_bm flux & maximize TAG prod.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:3
    tmp2 = tmp;
    tmp2.lb(62+i-1) = bmDrop*res.BM.optBM(i);
    if i == 1
        % no acetate uptake for auto
        tmp2.lb(28)   = 0;
    elseif i == 3
        tmp2.lb(1:12) = 0;
        tmp2.ub(1:12) = 0;
    end
    % select TAG production as obj function
    tmp2 = changeCbOpt(tmp2, length(tmp2.rxns)); 
    poo = optimizeCbModel(tmp2);
    res.TAG.optTAG(i) = poo.f;
    res.TAG.BM(i) = poo.x(62+i-1);
    res.TAG.x(:,i) = poo.x;
end





%% maxTAGprod_singleMOD
%
%
% IN: prism  - which prism reactions should be used as light source?
%     bmDrop - [0,1] what fraction of WT biomass should be kept?
%   OPTIONAL
%     kdown  - indicate extend of knock down instead of KO
%
% OUT: res <struct>
%       .bm     - opt BM value after single gene mod (KO or KD, if 0 < kdown < 1)
%       .TAG    - opt. FBA value for TAG production after single gene related 
%       rxn modification to kdown*ub(WT-affected-flux-value)
%       .TAGx   - complete flux distribution according to .TAG-FBA
%       
%
%
function res = maxTAGprod_singleMOD(prism, bmDrop, kdown)
global cr

% handle kdown
if ~exist('kdown','var')
    kdown = 0;
end

% 0.
%
% modify base network
% extract TAG production and formulate as extra function
tmp = modModel_extraTAG(cr);

% set PRISM 
tmp = setPRISM(prism, tmp);

% restrict uptake:
% ----------------
tmp.lb(27) = 0; % no starch
tmp.lb(17) = 0; % no no3
% tmp.lb(16) = 0; % no nh4


% 2. BM & TAG production after SingleKO
%
% biomass ...
% res.singleKO.WT_bm  = singleKO([4 5], tmp);

% TAG ... 
%   restrict lb-biomass to bmDrop of WT-biomass
for i = 1:3 % all biomass functions
    booBM = get_optCRbm(i,tmp); 
    tmp2 = tmp;
    tmp2.lb(62+i-1) = bmDrop*booBM; % modify current bm to carry at least bmDrop fraction
    if i == 1
        tmp2.lb(28) = 0; % no ac for autotrophic
    elseif i == 3
        tmp2.lb(1:12) = 0;
        tmp2.ub(1:12) = 0; % no light for heterotrophic growth
    end
    
    % TAG production after single KOs
    numGenes = length(tmp.genes);
    TAGid = length(tmp2.rxns);
    for j = 1:numGenes        
        [bohoo, foo, flag] = get_optFBA([62+i-1 TAGid], tmp2, [], j, kdown);
        if flag(1) % BM solution has been found
            res.bm(j,i)         = bohoo(1); % assign respective biomass
        else
            res.bm(j,i)         = NaN;
        end
        if flag(2) % TAG solution has been found ...
            res.TAG(j,i) = bohoo(2);
            res.TAGx(i, j, :)  = foo(:,2);
        else % ... else assign NaN as solution
            res.TAG(j,i) = NaN;            
            res.TAGx(i, j, :)  = NaN;
        end
        
    end
end





%% maxTAGprod_moduleMOD
%
%
% IN: prism  - which prism reactions should be used as light source?
%     bmDrop - [0,1] what fraction of WT biomass should be kept?
%     modulesPath - 
%   OPTIONAL
%     kdown  - [0,1] if given simulate knock down instead of KO (default = 0)
%              will be processed in rxn acitivity (lb/ub compared to wt)
%
% OUT:
%   res<struct>
%       .wtBM           - 
%       .wtTAG            -
%       .moduleKO.TAG   -
%       .moduleKO.bm    -
%
function res = maxTAGprod_moduleMOD(prism, bmDrop, modulesPath, kdown)
global cr

% handle kdown
if ~exist('kdown','var')
    kdown = 0;
end

% 0.
%
% modify base network
% extract TAG production and formulate as extra function
tmp      = modModel_extraTAG(cr);
TAGrxnID = length(tmp.rxns); % used for optimizations later on

% set PRISM 
tmp      = setPRISM(prism, tmp);


% restrict uptake:
% ----------------
tmp.lb(27) = 0; % no starch
tmp.lb(17) = 0; % no no3
% tmp.lb(16) = 0; % no nh4

% 0.1. compute biomass & TAG of WT after bmDrop to enable computation of 
% fractions later on
% --------------------------------------------------------------------------
for i = 1:3
    booBM = get_optCRbm(i,tmp); % handles specific BM function properties on its own (no ac for auto, no light for hetero)
    tmp2 = tmp;
    tmp2.lb(62+i-1) = bmDrop*booBM;     % restrict lb of WT_bm to bmDrop
    tmp2.ub(62+i-1) = booBM;            % for auto some constellations allow for huge bm increase
    % auto growth -> no acetate uptake
    if i == 1
        tmp2.lb(28) = 0;
    elseif i == 3
        tmp2.lb(1:12) = 0; % no light for hetero
        tmp2.ub(1:12) = 0; % no light for hetero
    end
    % biomass & TAG given the imposed restrictions to fluxes
    res.wtBM(i)   = booBM; % unnecessary: res.WTbm(i) = get_optCRbm(i,tmp2);
    res.wtTAG(i)  = get_optFBA(TAGrxnID, tmp2);
end
%}


% 1. perform KO of modules
% ------------------------
% identify KO genes along

% load Adrián's module data:
modules = load(modulesPath);
res.moduleKO.TAG  = zeros(size(modules,1), 3);
res.moduleKO.bm = zeros(size(modules,1), 3);


% maxJ in order to not compute every time dim length (see while loop below)
maxJ = size(modules,2);
for i = 1:size(modules, 1)
    
    % create KOidList related to the given modules
    % KOidList contains list of gene-IDs as given in the model
    KOidList = []; j = 2;
    while ~isempty(modules{i,j})
        % select only gene tag from Adrians data, like 'g1224'
        % add '.' to the end to make sure the right gene is identified
        KOcand = modules{i,j} ...
            (strfind(modules{i,j},'g'): ...
            strfind(modules{i,j},'_')-1);
        KOcand = strcat(KOcand,'.');
        
        % build up list of KO candidates
        KOidList = union(KOidList, find(...
            ~cellfun(@isempty,strfind(tmp.genes,KOcand))));
        j = j+1;
        if j > maxJ 
            break
        end
    end
    
    % if j == 2 the module was empty and the KO/KD has no effect on metabolism 
    if j == 2
        res.moduleKO.bm(i,:)  = res.wtBM;
        res.moduleKO.TAG(i,:) = res.wtTAG;
    else
        for m = 1:3        
            tmp2 = tmp;
            
            % restrict lb of WT_bm to bmDrop
            tmp2.lb(62+m-1) = bmDrop*res.wtBM(m);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % since res.wtBM = bmDrop * booBM (Erg von get_optCRbm(m,tmp))
            % booBM = get_optCRbm(m,tmp);   % DEPRECATED
            % see above (tmp2.lb(62+i-1) = bmDrop*booBM;):
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            if m == 1 % auto growth -> no acetate uptake
                tmp2.lb(28) = 0;
                % furthermore, restrict ub to WT value, since for auto, 
                % some constellations allow for huge bm increase...
                tmp2.ub(62+m-1) = res.wtBM(m); 
            elseif m == 3 % heterotrophic growth -> no light uptake
                tmp2.lb(1:12) = 0;
                tmp2.ub(1:12) = 0;
            end

            % compute FBA for bm and TAG 
            [poohoo, ~, flag] = ...
                get_optFBA([62+m-1, TAGrxnID], tmp2, [], KOidList, kdown);
            if flag(1) == 0
                % wtBM ist irreführend - bitte umbenennen
                res.moduleKO.bm(i,m) = NaN;
            else
                res.moduleKO.bm(i,m) = poohoo(1);
            end
            if flag(2) == 0
                res.moduleKO.TAG(i,m) = NaN;
            else
                res.moduleKO.TAG(i,m) = poohoo(2);
            end
        end
    end
end





%% modModel_extraTAG
%
% pull complete TAG production out of biomass function and set as extra
% function
%
% IN:  CR model, where TAG should be an extra reaction
%
% OUT: modified model
%
function res = modModel_extraTAG(model)

% biomass entries 128-169 for auto and mixo
% biomass entries 127-168 for hetero
% ausreichend das Ganze für auto zu machen, da die IDs für TAGs unabh./ von
% der BM rxn gleich sind
tmpIDs = find(model.S(:,62));
tmpIDs = tmpIDs(128:169);

% create new reaction for TAG production with S coefficients according to 
% biomass function
model = addReaction(model,'TAGprod', ...
    model.mets(tmpIDs), model.S(tmpIDs, 62), false);

% remove TAG from Biomass functions
model.S(tmpIDs, 62:64) = 0;

res = model;


%% moduleKO
%
% simuliere KO aller Module von Adrian und gib jeweilige BM zurück
%
% IN: prism - which prism
%   OPTIONAL
%     model 
%
% OUT: 
%     res.WT       - WT biomass
%     res.moduleKO - bm after module KO
%
function res = moduleKO(prism, modulesPath, model)
global cr

% check whether specific model is given, otherwise choose unmodified cr
% model
if ~exist('model', 'var')
    tmp = cr;
else
    tmp = model;
end

tmp = setPRISM(prism, tmp);

% load module data:
modules = load(modulesPath);

% restrict uptake:
% ----------------
tmp.lb(27) = 0; % no starch
tmp.lb(17) = 0; % no no3
% tmp.lb(16) = 0; % no nh4
% tmp.lb(28) = 0; % no ac


% 0. berechne biomass Opt von WT
res.WT(1) = get_optCRbm(1,tmp);
res.WT(2) = get_optCRbm(2,tmp);
res.WT(3) = get_optCRbm(3,tmp);

% 1. perform KO of modules
%   identify KO genes along
res.moduleKO = zeros(size(modules,1), 3);

% maxJ in order to not compute every time dim length (see while loop below)
maxJ = size(modules,2); 

for i = 1:size(modules, 1)
    tmpKO = tmp;
    j = 2;
    KOidList = [];
    while ~isempty(modules{i,j})
        % select only gene tag from Adrians data, like 'g1224'
        % add '.' to the end to make sure the right gene is identified
        KOcand = modules{i,j} ...
            (strfind(modules{i,j},'g'): ...
            strfind(modules{i,j},'_')-1);
        KOcand = strcat(KOcand,'.');
        % build up list of gene KO candidates
        KOidList = union(KOidList, find(...
            ~cellfun(@isempty,strfind(tmpKO.genes,KOcand))));
        j = j+1;
        if j > maxJ 
            break
        end
    end
    if j == 2 % j wie initialisiert -> daher kein Einfluss des KO auf NW
        res.moduleKO(i,:) = res.WT;
    else
        res.moduleKO(i,:) = get_BMgeneKO(KOidList, tmp);
        if res.moduleKO(i,1) == -1 
            res.moduleKO(i,:) = res.WT;
        end
    end
end





%% singleKO
%
% simulate all single KOs and influence on biomass
%
% IN:
%       prism - which prism reaction to use?
%   OPTIONAL
%       model - COBRA iRC1080 model
%
% OUT:
%       res.WT       - (1,1:3) vec, biomass w/o KO
%       res.singleKO - (length(genes),3) matrix, biomass after KO of one
%       gene at a time
%       res.nwEffect - <bool>(length(genes),1) vector,
%                    - 1 - gene KO affects nw, 0 - it does not
function res = singleKO(prism, model)
global cr

% check whether specific model is given, otherwise choose unmodified cr
% model
if ~exist('model','var')
    tmp = cr;
else
    tmp = model;
end
tmp = setPRISM(prism, tmp);

% restrict uptake: 
tmp.lb(27) = 0; % no starch
tmp.lb(17) = 0; % no no3
% tmp.lb(16) = 0; % no nh4
% tmp.lb(28) = 0; % no acetate ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. compute biomass from WT
res.WT(1) = get_optCRbm(1,tmp);
res.WT(2) = get_optCRbm(2,tmp);
res.WT(3) = get_optCRbm(3,tmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. get influence of singleKOs

% ----------------------------------
% res.singleKO pro biomass Def (auto, mixo, hetero)
res.singleKO = zeros(length(cr.genes), 3);

% -----------------------------------------------------------------
% res.nwEffect - bool: gene influences reaction -> 1, 0 otherwise
res.nwEffect = zeros(length(cr.genes),1);

for i = 1:length(res.singleKO)
%     tmpMod = tmp;

    tmpBM  = get_BMgeneKO(i, tmp);
    if tmpBM(1) == -1
        res.singleKO(i,:) = res.WT;
    else
        res.nwEffect(i)   = 1;
        res.singleKO(i,:) = tmpBM;
    end
end





%% get_BMgeneKO
%
% assuming complete KO
% translate gene KO to reaction KO
% give biomass values back
%
% IN:
%       KOs     = list of gene IDs that are going to be knocked out
%       model   = COBRA model structure
% OUT:
%       res     = (1,3) vector for three biomass values
%
function res = get_BMgeneKO(KOs, model)
res = zeros(1,3);

KOrxnsVec = logical(mapGeneKOtoRxns(model, KOs));

if sum(KOrxnsVec) == 0 % WT situation
    res = [-1 -1 -1];
else
    % gene KO causes reaction KO
    model.lb(KOrxnsVec) = 0;
    model.ub(KOrxnsVec) = 0;
    for i = 1:3
        res(i) = get_optCRbm(i,model);        
    end
end





%% get_BMopt
%
% calculates biomass WT value for all light resources and all biomass
% functions
%
% OUT: res<struct>
%   .opt:
%       1. Dim: active PRISM-reaction (length = 12) * 3 (one per biomass version)
%       2. Dim: biomass reaction (length = 3)
%               PRISM setting (12)
%               = 15
%   .xVec:
%       optimized flux vector
%
function res = get_BMopt()
global cr

% 12 PRISM reactions* #BM functions & 3 biomass functions + 12 PRISM ub status
res.opt  = zeros(12*3, 3+12);
% 3*12 xVec distribution
res.xVec = zeros(length(cr.c), 3*12);

tmpLB = cr.lb(1:12);
tmpUB = cr.ub(1:12);

for i = 3:3:36
    tmpCR = cr;
    % restrict uptake: 
    tmpCR.lb(27) = 0; % no starch
    tmpCR.lb(17) = 0; % no no3
%     tmpCR.lb(28) = 0; % no acetate ...
%     tmpCR.lb(24) = 0; % no O2 -> lethal
%     tmpCR.ub(13) = 0; % no h production
    
    tmpCR.lb(1:12) = 0;
    tmpCR.ub(1:12) = 0;
    % reassign light source flux
    tmpCR.ub(i/3)    = tmpUB(i/3);
    
    for j = 0:2
        [tmp.o, tmp.x]       = get_optCRbm(j+1, tmpCR);
        res.opt(i-2+j,j+1)   = tmp.o;
        res.opt(i-2+j,4:end) = tmpCR.ub(1:12)';
        res.xVec(:,(i-2)+j)  = tmp.x;
    end
        
end





%% get_model
%
% loads Papin model and hands it back
%
function res = get_model()

load ./model/iRC1080.mat
res = iRC1080;