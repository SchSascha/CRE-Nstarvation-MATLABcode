%% filterRes
%
% analysing modules of step 4 and 5, when running analysePapinModel
%


% IN:
%   model      - idR iRC1080
%   preResults - Output of analysePapinModel (most often step 4 or 5)
%   step       - redundant, needs to be congruent with 'preResults'
%   filterType - 1: singleKO, 2: moduleKO
%   fileName   - Name of csv output file
%  OPTIONAL
%   step5_bmLB - parameter for convenience
%              - e.g. lb10, lb30, oder lb50
%   modulePath - if filterType==2 required
% OUT:
%   res <cell> - csv ready cell array with all modules that show at least
%   5% drop in bm or TAG production
%
% Bsp.:
%   bm = find(schritt4_lb001_035.thres1.c.moduleKO.WTbm(:,2) / ...
%       schritt4_lb001_035.thres1.c.WTbm(2) < 0.95);
%   tag = -"-

function res = filterRes(model, preResults, schritt, filterType, fileName, step5_bmLB, modulePath)
global cr
cr = model;

% KO

% compute BM and TAG ratios

% step4vis   = visAnalysis(step4_lb001_035, 1);

switch schritt
    case 4  % KO simulation
        visAnRes  = visAnalysis(preResults, 1);
        threshold = fieldnames(preResults);
        fieldLength = length(threshold);       
    case 5  % KD simulation
        visAnRes  = visAnalysis(preResults, 2);
        visAnRes  = visAnRes.(char(step5_bmLB)); 
        threshold = fieldnames(visAnRes);
        fieldLength = 4;
    otherwise
        warning('Step not recognized. Options are 4 or 5.');
end

% build matrix across all thresholds for TAG and BM
if filterType == 1
    type = 'singleKO';
elseif filterType == 2
    type = 'moduleKO';
else
    error('Filtertype not recognized! 1-singleMOD, 2-moduleMOD');
end

for i = 1:length(threshold)
    TAG.auto(:,i)   = visAnRes.(char(threshold(i))).(char(type)).TAGratio(:,1);
    BM.auto(:,i)    = visAnRes.(char(threshold(i))).(char(type)).BMratio(:,1);
    TAG.mixo(:,i)   = visAnRes.(char(threshold(i))).(char(type)).TAGratio(:,2);
    BM.mixo(:,i)    = visAnRes.(char(threshold(i))).(char(type)).BMratio(:,2);
    TAG.hetero(:,i) = visAnRes.(char(threshold(i))).(char(type)).TAGratio(:,3);
    BM.hetero(:,i)  = visAnRes.(char(threshold(i))).(char(type)).BMratio(:,3);
end


clear stat

% check whether ALL 
stat.TAG(:,1) = min(TAG.auto,[],2);
stat.TAG(:,2) = min(TAG.mixo,[],2);
stat.TAG(:,3) = min(TAG.hetero,[],2);
stat.BM(:,1)  = min(BM.auto,[],2);
stat.BM(:,2)  = min(BM.mixo,[],2);
stat.BM(:,3)  = min(BM.hetero,[],2);

% How relates TAG to BM?
ID.bm      = (stat.BM  < 0.96);
ID.TAG     = (stat.TAG < 0.96);
ID.bmDrop  = (ID.bm - ID.TAG) == 1;
ID.TAGDrop = (ID.TAG - ID.bm) == 1;
ID.both    = (ID.bm + ID.TAG) == 2;



IDtype = {'bmDrop' 'TAGDrop' 'both'};

IDtag = {'at least one >5% BMdrop among conditions' ...
    'at least one >5% TAGdrop among conditions' ...
    'at least one condition with both, BM & TAG drop >5%' };


% build result cell Matrix:
% -------------------------

% header
% ------
res = cell(1,7+fieldLength); %22

if filterType == 1
    res(1,2) = {'single ID KO'};
else
    res(1,2) = {'module ID KO'};
end
for i = 7:(7+fieldLength-1)
    res(1,i) = threshold(i-6);
end

res(1,1) = {'biomass type'};
res(1,3) = {'gene'};
res(1,4) = {'rxn ID'};
res(1,5) = {'rxn'};
res(1,6) = {'pathway'};
% res(1,15) = {'0.01'};


BMtype = {'autotrophic' 'mixotrophic' 'heterotrophic'};

if filterType == 2
    % load module data:
    modules = load(modulesPath);
    maxJ = size(modules,2);
end

resNum = 2; % iterator für Resultate

% handle simulations of single gene modifications
% -----------------------------------------------
if filterType == 1
    for q = 1:3 % bm type
        IDtypeVisited = 1; % which result set? BMdrop, TAGdrop, oder both?
        for i = 1:3 % BMtype    
            tmp = find(ID.(char(IDtype(q)))(:,i));
            if ~isempty(tmp)
                % check for new line
                if IDtypeVisited == 1 % ID type changed -> new line
                    res(resNum,1) = IDtag(q);
                    res(resNum,2:end) = {''};
                    resNum = resNum+1;
                end
                IDtypeVisited = 0;
                % iterate all hits of current experiment:
                for n = 1:length(tmp)
                    % make sure TAG AND BM opt. produced feasible solutions
                    if ~isnan(stat.TAG(tmp(n),i)/stat.BM(tmp(n),i))
                        % tmp(n) contains the single gene KO/KD
                        KOrxnList = find(mapGeneKOtoRxns(model, tmp(n)));
                                
                        % fill result matrix
                        % ------------------
                           
                        res(resNum,1) = BMtype(i);
                        res(resNum,2) = {num2str(tmp(n))};
                        
                        % setup gene string
                        res(resNum,3) = model.genes(tmp(n));
                        
                        % setup rxn ID string
                        str = '';
                        for boo = 1:length(KOrxnList)
                            str = [str, num2str(KOrxnList(boo)), '; '];
                        end
                        res(resNum,4) = {str(1:end-2)};
                        % setup rxn string
                        str = '';
                        for boo = 1:length(KOrxnList)
                            str = ...
                                [str, model.rxnNames{KOrxnList(boo)}, '; '];
                        end
                        res(resNum,5) = {str(1:end-2)};    
                        % setup pw string
                        pw = unique(model.subSystems(KOrxnList));
                        str = '';
                        for boo = 1:length(pw)
                            str = ...
                                [str, model.subSystems{KOrxnList(boo)}, '; '];
                        end
                        res(resNum,6) = {str(1:end-2)};
                        % TAGs & BMs
                        for boo = 1:fieldLength
                            res{resNum,boo+6} = ...
                                visAnRes.(char(threshold(boo))).(char(type)).TAGratio(tmp(n),i) / ...
                                visAnRes.(char(threshold(boo))).(char(type)).BMratio(tmp(n),i);                    
                        end
                        resNum = resNum+1;
                    end
                end
        

            end
        end
    end
    
% handle simulations of module gene modifications
% -----------------------------------------------
elseif filterType == 2
    for q = 1:3 % bm type
        IDtypeVisited = 1;
        for i = 1:3 % BMtype    
            tmp = find(ID.(char(IDtype(q)))(:,i));
            if ~isempty(tmp)
                if IDtypeVisited == 1 % ID type changed -> new line
                    res(resNum,1) = IDtag(q);
                    res(resNum,2:end) = {''};
                    resNum = resNum+1;
                end
                IDtypeVisited = 0;
                for n = 1:length(tmp)
                    % make sure TAG AND BM opt. produced feasible solutions
                    if ~isnan(stat.TAG(tmp(n),i)/stat.BM(tmp(n),i)) 
                        % speed up, if result is already present:
                        [test.a, test.b] = ismember(num2str(tmp(n)), res(1:end-1,2));
                        % test whether result is present || biomass type
                        % was changed:
                        if ~test.a || ~strcmp(BMtype(i), res(test.b,1))
                            KOidList = []; j = 2;
                            while ~isempty(modules{tmp(n),j})
                                % select only gene tag from module data, like 'g1224'
                                % add '.' to the end to make sure the right gene is identified
                                KOcand = modules{tmp(n),j} ...
                                    (strfind(modules{tmp(n),j},'g'): ...
                                    strfind(modules{tmp(n),j},'_')-1);
                                KOcand = strcat(KOcand,'.');
                                % build up list of KO candidates
                                KOidList = union(KOidList, find(...
                                    ~cellfun(@isempty,strfind(model.genes,KOcand))));
                                j = j+1;
                                if j > maxJ 
                                    break
                                end
                            end
                            KOrxnList = find(mapGeneKOtoRxns(model, KOidList));

                            % fill result matrix
                            % ------------------

                            res(resNum,1) = BMtype(i);
                            res(resNum,2) = {num2str(tmp(n))};
                            % setup gene string
                            moduleStr = '';
                            for boo = 2:j-1
                                moduleStr = [moduleStr, modules{tmp(n),boo}, '; '];
                            end
                            res(resNum,3) = {moduleStr(1:end-2)};
                            % setup rxn ID string
                            moduleStr = '';
                            for boo = 1:length(KOrxnList)
                                moduleStr = [moduleStr, num2str(KOrxnList(boo)), '; '];
                            end
                            res(resNum,4) = {moduleStr(1:end-2)};
                            % setup rxn string
                            moduleStr = '';
                            for boo = 1:length(KOrxnList)
                                moduleStr = ...
                                    [moduleStr, model.rxnNames{KOrxnList(boo)}, '; '];
                            end
                            res(resNum,5) = {moduleStr(1:end-2)};
                            % setup pw string
                            pw = unique(model.subSystems(KOrxnList));
                            moduleStr = '';
                            for boo = 1:length(pw)
                                moduleStr = ...
                                    [moduleStr, model.subSystems{KOrxnList(boo)}, '; '];
                            end
                            res(resNum,6) = {moduleStr(1:end-2)};
                            % TAGs & BMs
                            for boo = 1:fieldLength
                                res{resNum,boo+6} = ...
                                    visAnRes.(char(threshold(boo))).(char(type)).TAGratio(tmp(n),i) / ...
                                    visAnRes.(char(threshold(boo))).(char(type)).BMratio(tmp(n),i);
                            end

                        else % module is already listed
                            priorID = find(strcmp(num2str(tmp(n)), res(:,2)));
                            priorID = priorID(1);
                            res(resNum,1) = BMtype(i);

                            for boo = 1:size(res,2)
                                res(resNum,boo) = res(priorID,boo);
                            end
                        end

                        resNum = resNum+1;
                    end
                end
            end 
        end
    end
end

if exist('fileName','var')
    cell2csv(strcat('../results/', fileName, '.csv'),res,'$');
end
