% replication code #2 in the sequence
% selects the Train/Validation/Test sets from the [455] documented kinases
% creates the LIK and LpcaS structures and identifies active ligands

% LIK --
% LpcaS -- PCA decomposition using only TRAIN set LG aff. vectors as basis
%          coordinate vectors are generated for all LG affinity vectors

%% option switches [some should be obsolete and unnecessary]

% switches to phase out ...
% TODO make sure this happens 

% ** ALWAYS FALSE **
NEW_PROJECTION = false; % use an alternative to the PCA basis
CENTER_TRAIN = false; % if true: the basis mapping has no offset
% ** ALWAYS FALSE **

% TRUE by default
REDO_PCA =  true;% PCA (based on Train set) / get coord vectors (for all) 

% ------------
% KEEP

% normally FALSE, keep as an option (?)
USE_TRACE_SPLIT = false; % if true, it uses the split from the previous tr object

SAVE_THIS_RUN = true;%false; % whether to save the trained net and related info
if ~exist('NN_TRAIN_SAVE_FOLDER','var')
    NN_TRAIN_SAVE_FOLDER = './';
end

%% input files

PML.TVTFile = [PML.DataPath,'TVTSets.xlsx']; % file with test/.. assignments
ReferenceFile = 'GeneralProcessedData.mat'; % file with general info
NNReferenceFile = 'samplerun.mat';  % file with a completed NN run

%% 

% if the main data is not present, load the reference file
%if ~exist('Tpca','var') || ~exist('Lpca','var') || ~exist('LCAMTs','var')
if ~exist('Tpca','var') || ~exist('LCAMTs','var')
    fprintf('Loading all data from %s ...',ReferenceFile)
    load(ReferenceFile);
    fprintf('done\n');
end

% choose the train/val/test split and run LIKsetup
if USE_TRACE_SPLIT
    % this is needed when we want to follow the same data partition as in
    % a reference file (typically a previous NN run)
    if ~exist('tr','var')
        fprintf('Loading training info only from %s ...',ReferenceFile)
        load(ReferenceFile,"tr")
        fprintf('done\n');
    end
else
    SetTable = readtable(PML.TVTFile);
    CurrentTVC = SetTable.Set52;
end


if ~exist('LIK','var') || ~isfield(LIK,'indivhits')
    LIKsetup
else
    % override any previous assignments in LIK
    LIK.TVC = CurrentTVC;
end

% perform PCA and construct coordinate vectors using the specified train
% set
RedoPCA;

