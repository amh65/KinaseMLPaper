% saves only the necessary information into two files
% #1 for general input data in processed form
% #2 for the trained neural net

% Change log
% 2024-09-28 Add linear regression model info to optional save items
% 2024-10-17 TODO: make the various switches optional

SAVE_GEN = false;%true;%

if SAVE_GEN
    GenFileName = 'GeneralProcessedData.mat';



    % collect and save generic info

    % important but not large

    % NT - number of all kinases
    % LigandInfoKinases - list of kinases with ligand binding gata
    % NLK -- length of the above (number of kinases with ligand info)
    % NLC - number of ligand groups

    % large structures
    % KGOT(NT,..) -- table with kinase identifiers and gene ontology (?) info

    % Tpca -- PCA info for the kinase input side


    % LCAMTred(NLK,NLC) -- kinase x ligand group activities used in most applications
    % LCAMTs(NLK,NLC) -- same as LCAMTred but normalized COLUMN-WISE to [0,1]

    save(GenFileName,...
        "NT","NLK","NLC",...
        "LigandInfoKinases",...
        "LCAMTred","LCAMTs",...
        "Tpca","Lpca",...
        "KGOT","-v7.3");

    return

end

%% net specific items

% name of the old style net file
OldNetFileName = mfn;
%OldNetFileFolder = NN_TRAIN_SAVE_FOLDER;

NewNetFileName = ['ShortPlus_',OldNetFileName];

% check for a trained nn + tr
if ~exist('net','var')
    fprintf([ReturnMessage, 'NN is missing!\n'])
    return;
elseif ~exist('tr','var')
    fprintf([ReturnMessage, 'NN present but training info is missing!\n'])
    return;
end


%% determine the major attributes of this run instance

% also filter by these attributes (only save the type we are interested in)

% MULTI_NET
if (~exist('MULTI_NET','var') || MULTI_NET==0) && net.input.size == 0
    if exist('mynets','var')
        MULTI_NET = true;%'tr','var'))
    else
        fprintf('Mislabelled multinet? - net input size is zero and there is no mynets!\n')
        return;
    end
end

% LASSO
% check for reduced ligand list mode ("lasso" mode) -- phase this out
if exist('RedLpca','var')
    LASSO_MODE = true;
    SNLC = length(SLG); %number of ligand groups after lasso
else
    LASSO_MODE = false;
end

if ~exist('NEW_PROJECTION','var')
    % option for alternative basis (instead of SVD) as in AltPCA.m
    % added circa 2024-03-24
    NEW_PROJECTION = false;
end


if ~exist('CENTER_TRAIN','var')
    % option to map the SVD coordinate vectors to [-1,1] instead of [0,1]
    % added circa 2024-03-27
    CENTER_TRAIN = false;
end

if ~exist('REDO_PCA','var') || (exist("OVERRIDE_REDO_PCA","var") && OVERRIDE_REDO_PCA)
    % use a PCA basis derived only from the training (+ validation) set , added circa 2024-04-24
    REDO_PCA=false;
end

RunNoteString = '';

if exist('LASSO_MODE','var') && LASSO_MODE, RunNoteString = [RunNoteString, 'LASSO ']; end
if exist('MULTI_NET','var') && MULTI_NET, RunNoteString = [RunNoteString, 'MULTI ']; end
if exist('NEW_PROJECTION','var') && NEW_PROJECTION, RunNoteString = [RunNoteString, 'AltBasis ']; end
if exist('CENTER_TRAI','var') && CENTER_TRAIN, RunNoteString = [RunNoteString, 'CenteredSVD ']; end
if exist('REDO_PCA','var') && REDO_PCA, RunNoteString = [RunNoteString, 'RedoPCA ']; end

% TODO : phase this out eventually
if ~strcmp(strtrim(RunNoteString),'RedoPCA')
    fprintf('Short file conversion only works for RedoPCA type files\n');
    return
end

%% collect and save net specific items
% net , tr -- main network file and training record
% Lpca -- basis vectors in the PCA scheme (depend on train/val/test split)


save(NewNetFileName,...
    "mfn",...
    "REDO_PCA",...
    "net","tr",...
    "Xdata","nmax1","nmax2",...
    "-v7.3");

OptVarList = {...
    "LASSO_MODE","MULTI_NET","NEW_PROJECTION","CENTER_TRAIN",...
    "Lpca",...
    "LpcaS","NLCS","NTK","NOK","TrainSet","OtherSet",...
    "Lnew",...
    "RedLpca","SLG","SNLC",...
    "mynets","nnets"...
    "nmax1_lin","nmax2_lin","BClin","CClin",...
    };

for ivar=1:length(OptVarList)
    ThisVar = OptVarList{ivar};
    if exist(ThisVar,"var"), save(NewNetFileName,ThisVar,"-append"); end
end




