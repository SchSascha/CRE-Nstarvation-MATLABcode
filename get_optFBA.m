%% get_optFBA(obj, model, KOlist)
%
% compute FBA per obj function at a time and give optimized value and x vec back
%
% PRE: solver must be loaded (simplest way: initCobraToolbox, if configurated)
% 
% IN: 
%       objID  - list of objective functions IDs (IDs refer to model rxns)
%       model  - a COBRA model structure
%   OPTIONAL
%       rxn_KOlist  - list of rxn IDs that have to be KO
%       gene_KOlist - list of gene IDs that have to be KO 
%                     (if given, rxn_KOlist is NOT given!)
%       FLUXmod       - if given simulate knocFLUXmod instead of KO
%                       - infers calculation of WT model without knock out
%
% OUT:
%       obj<vector>  - list of optimized FBA value i (length = length(objID))
%       xVec<vector> - (:,i): opt. flux distribution i
%       flag<vector> - 1 - feasible FBA
%                    - 0 - no feasible solution
%                    - length = length(objID)
%
% @Sascha Schäuble
%
function [obj, xVec, flag] = get_optFBA(objID, model, rxn_KOlist, gene_KOlist, FLUXmod)

% go through obj list and assign FBA value for each
obj  = zeros(length(objID), 1);
xVec = zeros(length(model.rxns), length(objID));
flag = zeros(length(objID), 1);

if exist('gene_KOlist','var') && ~isempty(rxn_KOlist)
    warning('Contradicting input data (cannot give gene AND rxn KO list)');
else
    % identify rxns that will be modified by KO/KD
    KOrxnsVec = [];
    if exist('gene_KOlist','var')
        KOrxnsVec = find(mapGeneKOtoRxns(model, gene_KOlist));
    elseif exist('rxn_KOlist','var')
        KOrxnsVec = rxn_KOlist;
    end
    
    for i = 1:length(objID)
        tmp = model;

        % check whether rxn-flux modifications (KO/KD/OE) need to be
        % applied
        if ~isempty(KOrxnsVec)
            if ~exist('FLUXmod', 'var') % KO assumed
                tmp.lb(KOrxnsVec) = 0;
                tmp.ub(KOrxnsVec) = 0;
            else
                % get WT flux distribution in order to apply modification
                % for 'KOrxnsVec' affected rxn-fluxes
                [~, WTxvec] = get_optFBA(objID(i), tmp);
                
                % knock down:
                if FLUXmod < 1
                    % modify fluxes, respect reversibility...
                    %
                    % identify WT negative fluxes for modification
                    negs = intersect(KOrxnsVec, find(WTxvec < 0));
                    % identify WT positive fluxes
                    pos  = intersect(KOrxnsVec, find(WTxvec > 0));
                    if ~isempty(negs)
                        tmp.lb(negs) = WTxvec(negs) * FLUXmod;
                    end
                    if ~isempty(pos)
                        tmp.ub(pos)  = WTxvec(pos)  * FLUXmod;
                    end
                % overexpression (Fix WT flux to given FLUXmod.):    
                elseif FLUXmod >= 1
                    tmp.lb(KOrxnsVec) = WTxvec(KOrxnsVec) * FLUXmod;
                    tmp.ub(KOrxnsVec) = WTxvec(KOrxnsVec) * FLUXmod;
                else
                    error('Negative flux modification is not allowed!');
                end
                
            end
        end
        
        
        % actual calculation:
        % -------------------
        tmp    = changeCbOpt(tmp, objID(i));
        FBAsol = optimizeCbModel(tmp);
        obj(i) = FBAsol.f;
        if ~isempty(FBAsol.x)
            xVec(:,i) = FBAsol.x;
            flag(i)   = 1;
        end
    end
end

    
    
    
    
% DEPRECATED code since 04/23/14
% remove if the above works fine
%{    
if exist('FLUXmod', 'var')
    [tmp, WTxvec] = get_optFBA(objID, model);
end

if exist('gene_KOlist','var')
    if ~isempty(rxn_KOlist)
        warning('Contradicting input data (cannot give gene AND rxn KO list)');
    else
        KOrxnsVec = logical(mapGeneKOtoRxns(model, gene_KOlist));
        if sum(KOrxnsVec) > 0
            if ~exist('FLUXmod', 'var')
                model.lb(KOrxnsVec) = 0;
                model.ub(KOrxnsVec) = 0;
            else
                model.lb(KOrxnsVec) = WTxvec(rxn_KOlist) * FLUXmod;
                model.ub(KOrxnsVec) = WTxvec(rxn_KOlist) * FLUXmod;
            end
        end
    end
% handle KO rxn candidates, if existing
elseif exist('rxn_KOlist','var')
    if ~exist('FLUXmod', 'var')
        model.lb(rxn_KOlist) = 0;
        model.ub(rxn_KOlist) = 0;
    else
        model.lb(rxn_KOlist) = WTxvec(rxn_KOlist) * FLUXmod;
        model.ub(rxn_KOlist) = WTxvec(rxn_KOlist) * FLUXmod;
    end
end




for i = 1:length(objID)
    tmp           = model;
    tmp           = changeCbOpt(tmp, objID(i));
    poo           = optimizeCbModel(tmp);
    obj(i)    = poo.f;
    if ~isempty(poo.x)
        xVec(:,i) = poo.x;
        flag(i)   = 1;
    end
end
%}
