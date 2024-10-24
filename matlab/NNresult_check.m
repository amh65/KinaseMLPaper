% computes measures to evaluate a trained NN for ligand affinity predictions
% usage: expects a completed NN run, for example
% >> load Processed/set4_100_1000_500_500_200_20_72s_485epochs_perf1.039e-03.mat

% global Rez  % structure to hold NN output and reconstructed affinities

% settings 

OVERRIDE_LASSO = false;%true;%
OVERRIDE_REDO_PCA = false;%true;
OMIT_LINREG = false;%true;%

LINEAR_TEST = false;%true;%  replaces linear regression output for the NN - PHASING OUT 
if LINEAR_TEST, SaveDir = [SaveDir,'LinReg']; end

NO_FIGURES = true;%false;%  % if true, no plots are made true;%true;%true;%
SAVEFIGURES = false;% true;%

%SaveDir = './scratch/Figs75004/';

 % ** Change Log Started Jan 12 2024 **
% 2024-01-08 - 2024-01-10
%    count predicted hits using a cutoff value
%    true / false positives and negatives (fig. 765_
%    TPR, FPR and ROC  curves (fig.100 101)
%    area under ROC curve fig. 102 103
% 2024-01-16 - 2024-01-17
%    added calculation of ROC curves for the entire ensemble ("global")
%    retiring individual (one-off) ROC calculation done earlier for "kk"
%    clean up the selection of kinases of interest
% 2024-01-19
%    update so it can also process results where only a subset of all
%    ligands is included in the NN calculation:
%
%    if the right objects exist (RedLpca and the list of ligand SLG), it
%    performs the analysis using the reduced ligand set etc. ...
%
%     ... unless OVERRIDE_LASSO is set to true; this forces the normal
%     analysis, using Lpca and the full set of ligands
%
%     REMOVED fig. 765 that used LIKhitcounts and also the RefScore
%     value used only tere
%
%  2024-02-02
%     Cleaning up
%
%
% 2024-02-16
%
%     TODO: estimate ROC for the SVD-only predictions - done 2/22
%     TODO: choose a cutoff (global, maybe also by kinase)
%           -- global cutoff done 2/23
%           -- individual cutoff done 2/26 more or less
%     TODO: for a given cutoff, count TP, FP, etc.
%           output lists of hits / misses (ligand groups)
%     TODO: replace SVD basis with something better (PC Regression etc)
%
%     TODO: (minor) set a single (original) cutoff value and propagate it throughout,
%     instead of using LCAMTred<4 each time
%
%     TODO: add an initial check / warning to make sure the initial data is
%     there
%
%     2024-02-21 cleaning up: loop-free calculation of LIKPredHits.mean etc.
%
%     2024-02-23 added global cutoff and by-kinase TP, FP, etc. based on
%     that
%
%     2024-03-15 merging two new options in the NN pipeline
%                (1) multiple nets -- done circa 2024-03-22 or before
%                (2) alternative to PCA -- done circa 2024-03-28 or before
%
%     2024-03-30
%                 adding plain linear regression for comparison
%                 cleaning up figures
%
%     2024-04-01 - ...
%                 cleaned up figures
%                 added measures for hit identification efficiency
%                 also new figure 310
%
%    2024-04-24 - new develop version (previous was frozen as 2024-04-01.m
%
%                 adding REDO_PCA flag to use reduced SVD vectors from
%                 RedoPCA.m and nn_Test_develop_2024_04_22.m
%
%                 moved LIKpredhits , LIHprednonhits into sub-structs in
%                 Lpca: Lpca.PredHits, Lpca.PredNonHits
%                 for the REDO_PCA option, these are under LpcaS
%
% 2024-05-23 saved the previous develop version with this date, current one
% continues as develop version
%
% 2024-05-28 reorganize the analysis to make it more transparent
%
%            TODO: Compute all stats for the three (TVC) groups separately
%             added ThreeSets{#} #=1,2,3 with lists of kinase indices
%             added Linds3{,} 3x2 cell similar to Linds 
%
%             renamed *pca.checkP --> *pca.checkNN
%
%            TODO: Minimal structure to hold items relating to a model run
% 
%               New object: Rez ---> RSA
%               Use a function to populate it


%% NOTES NO CODE: ** Main variables etc. **
%
%  net, tr -- pre-existing trained NN and trace [net, tr] = train(neT)
%
%  LigandInfoKinases -- list of kinase indices that have ligand information
%  NLK -- size of the above
%  NT -- number of all kinases
%
%  LCAMT(NT,NLC) -- group affinities by kinase; zero for undocumented ones
%  LCAMTred(NLK,NLC) -- same as above but only the ligand info kinases
%
%  "hits" are defined as LCAMT > 4
%
%  LIK -- general info on documented kinases
%       .interactions(NLK) -- number of "hits" by kinase
%       .level(NLK) -- 1 or 2 corresp to Tclin resp. Tchem
%       .TVC(NLK) -- 1 (train), 2 (val), 3 (test) -- from NN training
%       .lind -- index in LigandInfoKinases [DELETED]
%       .trueind -- kinase index in the general list
%       NOTE: LIKTable is a table containing the above fields of the LIK struct
%       .sort(NLK) -- kinases sorted by number of interactions (descending)
%       .unsort(NLK) -- inverse of the above
%                        LIK.unsort(i) = rank of kinase i by int.count
%       .nonhits(NLK) -- number of nonhits by kinase
%                        same as NLC - LIK.interactions
%       .indivhits(NLK,NLC)   -- logical flags, 1 if LCAMTred(i,j)>4
%       .indivnonhits(NLK,NLC) -- 1 - the above
%       .AUC_NN(NLK,1) -- area under the ROC curve, using the NN pred
%       .AUC_SVD(NLK,1) -- area under the ROC curve, using SVD-only
%
% Lpca -- PCA and other info by kinase
%           this was originally set up specifically for the PCA output
%           part of it is initially calculated in the main sequence,
%           specifically in KinaseGOCheck in preparation for training a NN
%
%   .coeff -- this one and up to .mu are outputs of the pca(LCAMTs) function
%   .score -- only coeff and score are used in the NN stuff
%   .latent -- not used
%   .tsquared -- not used
%   .explained -- used for a plot in KinaseGOCheck
%   .mu -- not used
%   .mean(NLK,NLC) -- mean by columns expanded to match the size of LCAMTs;
%       Lpca.mean = repmat(mean(LCAMTs,1),size(LCAMTs,1),1);
%       used to reconstruct the original vectors:
%   .check(NLK,NLC) -- reconstructed vectors
%       Lpca.check = Lpca.score * Lpca.coeff' + Lpca.mean;
%       max(max(LCAMTs - Lpca.check))
%   .checkSVD(NLK,NLC) -- Lpca.check renamed in NNresult)check
%       NOTE: in NNresult_check only a subset of the PCA components are
%       used so Lpca.check is a lossy recontrsuction, that can be used to
%       estimate the loss due to PCA truncation
%
%   Struct similar to Lpca: Lnew (for NEW_PROJECTION) and LpcaS (REDO_PCA)



%% mode selection and checks
% expect to run this after basic data and a trained NN have been loaded
ReturnMessage = 'NNresult_check should be run after basic data and a trained NN have been loaded -- ';
if ~exist('KGOT','var')
    fprintf([ReturnMessage, 'Kinase input information is missing!\n'])
    return;
elseif ~exist('LigandInfoKinases','var')
    fprintf([ReturnMessage, 'Ligand information is missing!\n'])
    return;
elseif ~exist('net','var')
    fprintf([ReturnMessage, 'NN is missing!\n'])
    return;
elseif ~exist('tr','var')
    fprintf([ReturnMessage, 'NN present but training info is missing!\n'])
    return;
elseif (~exist('MULTI_NET','var') || MULTI_NET==0) && net.input.size == 0 
        if exist('mynets','var')
            MULTI_NET = true;%'tr','var'))
        else
            fprintf('Mislabelled multinet? - net input size is zero and there is no mynets!\n')
            return;
        end
end


% check for reduced ligand list mode ("lasso" mode) -- phase this out
if exist('RedLpca','var') & ~OVERRIDE_LASSO
    LASSO_MODE = true;
    SNLC = length(SLG); %number of ligand groups after lasso
else
    LASSO_MODE = false;
end


if ~exist('MULTI_NET','var')
    % multi-net option added 2024-03-15
    % idea is to have multiple "parallel" NN objects that predict subsets of
    % the SVD (or generally, reduced) representation of the outputs
    MULTI_NET = false;
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

if ~exist('REDO_PCA','var') || OVERRIDE_REDO_PCA
    % use a PCA basis derived only from the training (+ validation) set , added circa 2024-04-24
    % the "shortbasis" results use only the training set for the PCA basis;
    % this is standard in the paper (to be submitted circa September 2024)
    REDO_PCA=false;
end


RunNoteString = '';
if LASSO_MODE, RunNoteString = [RunNoteString, 'LASSO ']; end
if MULTI_NET, RunNoteString = [RunNoteString, 'MULTI ']; end
if NEW_PROJECTION, RunNoteString = [RunNoteString, 'AltBasis ']; end
if CENTER_TRAIN, RunNoteString = [RunNoteString, 'CenteredSVD ']; end
if REDO_PCA, RunNoteString = [RunNoteString, 'RedoPCA ']; end   


%% recover the run code from the saved file name 
% try to identify the run 'Code' - a random integer that is part of the
% file name; the file name patterns is usually
mfnr = mfn;
RunCode = '';
for itk = 1:4+size(net.layers,1)
    [RunCode, mfnr] = strtok(mfnr,'_');
end

SaveDir = ['./scratch/Figs',RunCode,'/'];
if ~exist(SaveDir,'dir') & SAVEFIGURES
    mkdir(SaveDir);
end


%% set up the LIK struct - may be unnecessary

if ~exist('LIK','var') || ~isfield(LIK,'indivhits') || ~isfield(LIK,'TVC') % || ~exist('ThreeSet','var')
    LIKsetup
end

%% use Lthis instead of Lpca to accommodate alt projection 
%  and unified processing in general

if NEW_PROJECTION % import the main structure with vectors , bases, coordinates
    Lthis = Lnew;
elseif REDO_PCA
    % check that we have a SHORT BASIS (i.e. only the train set is used for
    % the PCA basis)
    if NTK > length(tr.trainInd)
        fprintf('*** Basis kinases: %d Train Kinases: %d -- skipping this run\n',NTK,length(tr.trainInd));
        return
    end
    if ~isfield(LpcaS,'TrainLig') || ~exist('LCAMTb','var') || ~exist('NLCS','var')
        % which ligands are [not] hit by the train set ?
        LCAMTb = LCAMTs; % binary LCAMT, only 1's and 0's
        LCAMTb(LCAMTs>0) = 1;

        % identify ligands that interact with the train set
        LpcaS.TrainLig = find(sum(LCAMTb(TrainSet,:),1)>0);
        NLCS = length(LpcaS.TrainLig);
    end

    if ~isfield('LpcaS','OtherLig')
        LpcaS.OtherLig = setdiff(1:NLC, LpcaS.TrainLig);
    end
    Lthis = LpcaS;
else
    Lthis = Lpca;
end

if REDO_PCA % active ligands & kinases used vary by case
    Lthis.ActiveLigands = LpcaS.TrainLig;
    Lthis.ActiveLigandCount = NLCS;
    %Lthis.NSV = NTK - 1 ; % # of singular vectors 
    % max # of singular vectors to use
    Lthis.NSV = rank(LCAMTs(TrainSet,:)); % changed 2024-09-27
elseif LASSO_MODE

    Lthis.ActiveLigands = SLG;
    Lthis.ActiveLigandCount = SNLC;

    % some setup is done in LIKsetup
    if ~isfield(LIK,'indivhits_red')
        LIKsetup;
    end

    Lthis.NSV = NLK - 1; % number of singular vectors
else % default (all ligands are considered "active")
    Lthis.ActiveLigands = 1:NLC;
    Lthis.ActiveLigandCount = NLC;
    Lthis.NSV = NTK - 1;
end


%% NOTES NO CODE: original, SVD projected, and NN predicted ligand interaction vectors

%  LCAMT has the max activities (of any type) by cluster
%        it is constructed in Pharos_test.m
%        data is available only for LigandInfoKinases so we use
%        LCAMTred = LCAMT(LigandInfoKinases,:)

% the SVD projection is initially done in KinaseGOCheck.m
%       LCAMTred is defined there
%       then scaled so
%       each column in LCAMTs [corresp to a ligand group] is in [0,1]

% LCAMTs is used in the pca function - output is in the Lthis struct
%       Lpca.check = Lpca.score * Lpca.coeff' + Lpca.mean
%       max(max(abs(LCAMTs - Lpca.check))) should be [machine] zero
%       !! only if there is no truncation !!
%
% The PCA vectors are in Lpca.score [in the rows]
%      These are also scaled to [0,1] by column [each column is in [0,1] ]
%      Lpca.sscore is the scaled version
%
% NN training is done with Ydata = Lpca.sscore(:,1:nmax2)
%       i.e. keeping only the first nmax2 columns [components]
%
%       YY = net(Xdata') is the prediction [regression result] using the
%       trained NN
%% retrieve NN (net) parameters and generate identifier string

if ~MULTI_NET
    thisnet = net;
    thistr = tr;
    multiprefix ='';
else
    thisnet = mynets.net{1};
    thistr = mynets.tr{1};
    multiprefix = sprintf('multi_%dx_',nnets);
end

% for backward compatibility ...
if ~exist('RAW_SVD','var'), RAW_SVD=false; end

netstring = sprintf('%d',thisnet.inputs{1}.size);
for il=1:thisnet.numLayers
    netstring = [netstring sprintf('_%d',thisnet.layers{il}.size)];
end

netstring = [multiprefix netstring];
% the outputs form the final layer so that one is already included
netfileID = [netstring sprintf('_%depochs',thistr.num_epochs)];

% use the RunCode (just random digits) as the [neural net] run / file ID
if exist('RunCode','var')
    netfileID = RunCode;
end

%save(netfilename,'-v7.3');

%% construct predicted affinity vectors

% get the NN predicted Y using a pre-trained net
if ~MULTI_NET
    YY = net(Xdata');
    YY = YY';
    nmaxx2 = nmax2;
else
    YY = zeros(NLK,nnets*nmax2);
    for inet = 1:nnets
        %thisnet = mynets.net{inet};
        %rr = mynets.rr{inet};
        thisYY = mynets.net{inet}(Xdata');
        YY(:,mynets.rr{inet}) = thisYY';
    end
    nmaxx2 = nnets * nmax2;
end


% map back to the original variables ..
if ~LASSO_MODE
    if RAW_SVD
        % if these are the NN predicted UNscaled PCA vectors ("scores")
        Lthis.scoreNN = YY;
    else
        % these are the NN predicted scaled PCA vectors ("scores")
        if CENTER_TRAIN
            Lthis.sscoreNN = (YY + 1)*0.5;
        else
            Lthis.sscoreNN = YY;
        end


        if ~isfield(Lthis,'minscore')
            Lthis.minscore = repmat(min(Lthis.score,[],1),size(Lthis.score,1),1);
            Lthis.maxscore = repmat(max(Lthis.score,[],1),size(Lthis.score,1),1);
        end

        % undo the [0,1] scaling on the PCA side
        % scoreP is the predicted, scaled set of scores
        Lthis.scoreNN = Lthis.minscore(:,1:nmaxx2) + ...
            Lthis.sscoreNN .* (Lthis.maxscore(:,1:nmaxx2) -...
            Lthis.minscore(:,1:nmaxx2));


    end

    % map predicted affinities back to the original representation (undo the PCA projection)
    Lthis.checkNN = Lthis.mean + Lthis.scoreNN * Lthis.coeff(:,1:nmaxx2)';

    % for comparison, also reconstruct the pre-PCA vectors using only nmax2 principals
    Lthis.checkSVD = Lthis.mean + Lthis.score(:,1:nmaxx2) * Lthis.coeff(:,1:nmaxx2)';

    % compute linear regression prediction (can be turned off)
    if ~OMIT_LINREG || LINEAR_TEST
        %linear regression for comparison
        linreg;
        Lthis.scoreLin = Lthis.minscore + Lthis.sscoreLin .* (Lthis.maxscore - Lthis.minscore);
        Lthis.checkLin = Lthis.mean + Lthis.scoreLin * Lthis.coeff';
    end

    if LINEAR_TEST
        %replace linear regression for the NN prediction
        Lthis.checkNN = Lthis.checkLin;
        fprintf('** caution ** LINEAR REGRESSION FOR COMPARISON - %d COMPONENTS\n',nmax2_lin);
    else
        % verification output
        fprintf('** Train:%d Val:%d Test:%d -- 2 groups: Train:%d Other:%d\n',...
            length(find(LIK.TVC==1)),length(find(LIK.TVC==2)),length(find(LIK.TVC==3)),...
            length(TrainSet),length(OtherSet));
    end

    if REDO_PCA
        % extra safety, make sure the predicted affinities of
        % non-participating ligands are zero
        Lthis.checkNN(:,Lthis.OtherLig)=0;
        Lthis.checkSVD(:,Lthis.OtherLig)=0;
        if isfield('checkLin','Lthis')
            Lthis.checkLin(:,Lthis.OtherLig)=0;
        end
    end

    % NOTE:
    % Lpca.checkNN (the NN prediction based on truncated SVD) and
    % Lpca.check (based on the exact truncated SVD)
    %    comparable to the *scaled* activity vectors in LCAMTs

else % lasso mode
    % map back to the original variables
    RedLpca.sscoreNN = YY;

    % undo the [0,1] scaling on the PCA side
    RedLpca.scoreNN = RedLpca.minscore(:,1:nmaxx2) + ...
        RedLpca.sscoreNN .* (RedLpca.maxscore(:,1:nmaxx2) -...
        RedLpca.minscore(:,1:nmaxx2));

    % map back to the original representation
    RedLpca.checkNN = RedLpca.mean + RedLpca.scoreNN * RedLpca.coeff(:,1:nmaxx2)';

    % reconstruct the pre-PCA vectors using only nmax2 principals
    RedLpca.check = RedLpca.mean + RedLpca.score(:,1:nmaxx2) * RedLpca.coeff(:,1:nmaxx2)';

    % analysis to compare with the regular (full set NN model outputs)

    % put the reconstructed values into full size interaction arrays
    Lpca.checkNN = zeros(size(LCAMTred));
    Lpca.checkSVD = Lpca.checkNN;

    Lpca.checkNN(:,SLG) = RedLpca.checkNN;
    Lpca.checkSVD(:,SLG) = RedLpca.check;

    Lthis.checkNN = Lpca.checkNN;
    Lthis.checkSVD = Lpca.checkSVD;

end

%% stats of hit / nonhit predicted values by kinase


% LIKPredHits etc. -- renamed Lpca.PredHits 
%
% V2: compute the range and mean of predicted affinities for hit / nonhits
% updated to loop-free version 2024-02-23
%
% moving LIK<...> to Lpca.<...> 2024-04-26 eg. LIKPredHits is now
% Lpca.PredHits (ie. a sub-struct)
% PredHits -- predicted affinities for the hit ligands
% PredNonHits -- same but for the non-hit ligands

if ~REDO_PCA & ~LASSO_MODE
    % default case, all ligands are "active"

    if ~isfield('Lpca','interactions')
        Lpca.interactions = LIK.interactions;
    end

    % interaction count by kinase, restricted as applicable
    Lthis.FilteredInteractions = Lpca.interactions;
    Lthis.FilteredNonHits = NLC - Lpca.interactions;

    % semi-binary - predicted affinities for hits only
    PredHits.all = Lthis.checkNN .* LIK.indivhits; % 0 if nonhit, value if hit

    PredHits.mean = sum(PredHits.all,2)./Lthis.FilteredInteractions ;
    ww = PredHits.all; ww(ww==0)=nan;
    PredHits.low = min(ww,[],2);
    PredHits.high = max(ww,[],2);

    PredNonHits.all = Lthis.checkNN.*LIK.indivnonhits;
    PredNonHits.mean = sum(PredNonHits.all,2)./Lthis.FilteredNonHits;
    ww=PredNonHits.all; ww(ww==0)=nan;
    PredNonHits.low = min(ww,[],2);
    PredNonHits.high = max(ww,[],2);

    Lthis.PredHits = PredHits;
    Lthis.PredNonHits = PredNonHits;

elseif REDO_PCA
    % hits, interaction counts etc. restricted to the ligands interacting
    % with the train (+ val) kinase set
    Lthis.FilteredInteractions = sum(LCAMTb(:,LpcaS.TrainLig),2); % int. counts
    Lthis.FilteredNonHits = NLCS - Lthis.FilteredInteractions; % counts of noninteracting ligands

    % interaction count by kinase, restricted as applicable
    % Lthis.FilteredInteractions = LpcaS.interactions;
    % Lthis.FilteredNonHits = LpcaS.nonhits;

    PredHits.all = Lthis.checkNN .* LIK.indivhits; % 0 if nonhit, value if hit
    PredHits.all = PredHits.all(:,LpcaS.TrainLig); % affinities

    PredHits.mean = sum(PredHits.all,2) ./ Lthis.FilteredInteractions;
    ww = PredHits.all; ww(ww==0)=nan;
    PredHits.low = min(ww,[],2);
    PredHits.high = max(ww,[],2);

    PredNonHits.all = Lthis.checkNN.*LIK.indivnonhits;
    PredNonHits.all = PredNonHits.all(:,LpcaS.TrainLig); % affinities

    PredNonHits.mean = sum(PredNonHits.all,2) ./ Lthis.FilteredNonHits ;
    ww = PredNonHits.all; ww(ww==0)=nan;
    PredNonHits.low = min(ww,[],2);
    PredNonHits.high = max(ww,[],2);

    Lthis.PredHits = PredHits;
    Lthis.PredNonHits = PredNonHits;

elseif LASSO_MODE

    % interaction count by kinase, restricted as applicable
    Lthis.FilteredInteractions = LIK.interactions_red;
    Lthis.FilteredNonHits = LIK.nonhits_red;

    PredHits.all = RedLpca.checkNN .* LIK.indivhits_red;

    PredHits.mean = sum(PredHits.all,2) ./ Lthis.FilteredInteractions;
    ww = PredHits.all; ww(ww==0)=nan;
    PredHits.low = min(ww,[],2);
    PredHits.high = max(ww,[],2);

    PredNonHits.all = RedLpca.checkNN .* LIK.indivnonhits_red;

    PredNonHits.mean = sum(PredNonHits.all,2) ./ Lthis.FilteredNonHits ;
    ww = PredNonHits.all; ww(ww==0)=nan;
    PredNonHits.low = min(ww,[],2);
    PredNonHits.high = max(ww,[],2);

    Lthis.PredHits = PredHits;
    Lthis.PredNonHits = PredNonHits;

end

clear ww; clear PredHits; clear PredNonHits;

%% TPR, FPR etc. curves - for the entire range of possible cutoffs

% sort the kinases by the predicted affinity then plot the (normalized)
% index on y and the predicted affinity on x
%
% reminder -- how to use sort()
% [B,I] = sort(___) also returns a collection of index vectors .
% I is the same size as A and describes the arrangement
% of the elements of A into B along the sorted dimension.
% For example, if A is a vector, then B = A(I).

if LASSO_MODE
    [RedLpca.sort_checkNN, RedLpca.predsort] = sort(RedLpca.checkNN,2,'descend');
    [RedLpca.sort_checkSVD, RedLpca.predsortSVD] = sort(RedLpca.check,2,'descend');
elseif REDO_PCA
    % this is necessary because we reshuffle the ligands of interest only
    % into the same places
    [LpcaS.sort_checkNN, LpcaS.predsort] = sort(Lthis.checkNN(:,Lthis.TrainLig),2,'descend');
    [LpcaS.sort_checkSVD, LpcaS.predsortSVD] = sort(Lthis.checkSVD(:,Lthis.TrainLig),2,'descend');
else %default case
    [Lthis.sort_checkNN, Lthis.predsort] = sort(Lthis.checkNN,2,'descend');
    [Lthis.sort_checkSVD, Lthis.predsortSVD] = sort(Lthis.checkSVD,2,'descend');% same as above for SVD only
end

% RezScaledAffinities -> RSA for short
RSA = cell(4,1); % will hold pred affinities and stats for NN,SVD,Lin
% RezGlobal -> RG for short % for global stats
% (indep, of method but specific to the
% model, typically dependent on the choce of active ligands)
clear RG ;% reinitialize to avoid confusion

% multiple result structures RSA{#}
if REDO_PCA
    %RSA{1}.label = 'NN';
    RSA{1}.label = sprintf('    NN %3d %3d ',nmax1_lin,nmax2_lin);
    %RSA{2}.label = 'SVD';
    RSA{2}.label = sprintf('   SVD %3d %3d ',nmax1_lin,nmax2_lin);
    RSA{3}.label = sprintf('   Lin %3d %3d ',nmax1_lin,nmax2_lin);
    %RSA{4}.label = 'Combo';
    RSA{4}.label = sprintf(' Combo %3d %3d ',nmax1_lin,nmax2_lin);
    % all RSA ligand fields restricted to active ligands
    RSA{1}.check = Lthis.checkNN(:,Lthis.ActiveLigands);
    RSA{2}.check = Lthis.checkSVD(:,Lthis.ActiveLigands);
    RSA{3}.check = Lthis.checkLin(:,Lthis.ActiveLigands);

    RSA{4}.check = max(RSA{1}.check,RSA{3}.check);%0.5*(RSA{1}.check+RSA{3}.check);%
    % RezGlobal for common items
    % will replace LIKall.hits etc.
    RG.AllHits = reshape(LIK.indivhits(:,Lthis.ActiveLigands),1,[]); % hit  flag 1 if hit, 0 if not
    RG.AllNonHits = reshape(LIK.indivnonhits(:,Lthis.ActiveLigands),1,[]); % nonhit flag, same idea
    % NOTE: RG.AllHits/ .AllNonHits are the same as LIKall.hits / .nonhits
    for iset = 1:3
        % approach (1) index vectors of the same size NLK x Lthis.ActiveLigands
        % matching RSA{#}.check
        % obtain them from Linds3{3,2} [ iset x hit / nonhit ]
        RG.HitsbySet{iset} = reshape(Linds3{iset,1}(:,Lthis.ActiveLigands),1,[]);
        RG.NonHitsbySet{iset} = reshape(Linds3{iset,2}(:,Lthis.ActiveLigands),1,[]);
        RG.HitCountbySet(iset) = sum(RG.HitsbySet{iset});
        RG.NonHitCountbySet(iset) = sum(RG.NonHitsbySet{iset});
        % % approach (2)index vectors to fit RSA{imethod}.predsbyset{iset}
        % each with size defined by the ligand set
        % RG.HitsbySet{iset} = reshape(LIK.indivhits(ThreeSets{iset},Lthis.ActiveLigands),1,[]);
        % RG.NonHitsbySet{iset} = reshape(LIK.indivnonhits(ThreeSets{iset},Lthis.ActiveLigands),1,[]);
    end


    for imethod = 1:4
        % sort predicted affinities and note the ordering
        % sort in DESCENDING order of predicted values because the TPR / FPR
        % are based on the number of hits / nonhits ABOVE the cutoff
        [RSA{imethod}.sort_check, RSA{imethod}.predsort] = sort(RSA{imethod}.check,2,'descend');
        % reshape the predicted value arrays for use in global stats
        % will replace LIKall.preds
        RSA{imethod}.allpreds = reshape(RSA{imethod}.check,1,[]);
        [RSA{imethod}.allsort_vals, RSA{imethod}.allsort_pos] = sort(RSA{imethod}.allpreds,'descend');
        % approach (1) -- use the full set of predicted values (all
        % kinases x active ligands, flattened then sorted
        % use the index vectors based on Linds3 and reorder them

        % this would work with approach (2)
        % for iset=1:3
        %     % will sort of replace LIKtrain LIKrest
        %     % matching index vectors in RG.[Non]HitsbySet{iset}
        %     RSA{imethod}.predsbyset{iset} = reshape(RSA{imethod}.check(ThreeSets{iset},:),1,[]);
        % end
    end

end
%% generate TPR, FPR etc curves by kinase
% NO CODE explanation:
% we need:
% TP[R] = True Positive [Rate] ("rate" is relative to all hits)
% FP[R] = False Positive [Rate] ("rate" is relative to all NONhits)
% various measures derived from these:
% ROC curve is TPR(cutoff) vs. FPR(cutoff)
% Youden (Youden's J)
% YJ = TPR - FPR; its maximum location used as an optimal cutoff value
% YU = TP - FP (i.e. the difference between non-normalized TP - FP
%
% Y is for the Youden statistic, TPR - FPR used to find an optimal cutoff
%
% will TP< FP, etc. generate it by counting the hits and nonhits up to a given predicted
% affinity value (as obtained in Lpca.checkNN sorted in Lpca.sort_checkP)

%% *** OPTION: RUN FROM HERE TO END if using a saved NN + LR file ***

%% 

% initialize
% pad the cumulative sum so we start from all zeros
FPRbyKinase = zeros(NLK,Lthis.ActiveLigandCount+1);
TPRbyKinase = FPRbyKinase;

FPRbyKinaseSVD = FPRbyKinase;
TPRbyKinaseSVD = FPRbyKinase;

YbyKinase = FPRbyKinase;
YbyKinaseSVD = FPRbyKinase;

for imethod = 1:4
    RSA{imethod}.FPRbyKinase = zeros(NLK,Lthis.ActiveLigandCount+1);
    RSA{imethod}.TPRbyKinase = zeros(NLK,Lthis.ActiveLigandCount+1);
    RSA{imethod}.YbyKinase = zeros(NLK,Lthis.ActiveLigandCount+1);
    RSA{imethod}.hitsintop = zeros(NLK,1); % number of hits in top X%
    RSA{imethod}.FPfortophits = zeros(NLK,1); % FP rate for top Y of hits
end

% also do the recall and precision numbers here

% Recall / precision: 
% How many hits are in the top X% [...] of all predictions 
% (i.e. top X% of all active ligands, hit or nonhit)
Lthis.topfraction = 0.050;%0.010;%   %fraction X% we want; 
Lthis.topind = ceil(Lthis.topfraction*Lthis.ActiveLigandCount); % index of that fraction

Lthis.hitsintop = zeros(NLK,1); % will hold the answer (# hits in top x%)

% How many nonhits (false positives) to capture y% of true hits
Lthis.hitfraction = 0.5;% fraction y% we want
Lthis.hitind = ceil(Lthis.hitfraction * Lthis.FilteredInteractions); % # of hits we are looking to recover

Lthis.FPfortophits = zeros(NLK,1); % will hold the answers (# FP above y% of hits)

if LASSO_MODE % Lasso mode

    for kin=1:NLK

        NNinds  = RedLpca.predsort(kin,:);
        SVDinds = RedLpca.predsortSVD(kin,:);

        TPRbyKinase(kin,2:end)      = cumsum(LIK.indivhits_red(kin,NNinds)) / LIK.interactions_red(kin) ;
        FPRbyKinase(kin,2:end)      = cumsum(LIK.indivnonhits_red(kin,NNinds)) / LIK.nonhits_red(kin);

        TPRbyKinaseSVD(kin,2:end)      = cumsum(LIK.indivhits_red(kin,SVDinds)) / LIK.interactions_red(kin) ;
        FPRbyKinaseSVD(kin,2:end)      = cumsum(LIK.indivnonhits_red(kin,SVDinds)) / LIK.nonhits_red(kin);

        YbyKinase(kin,2:end) = cumsum(LIK.indivhits_red(kin,NNinds)) - cumsum(LIK.indivnonhits_red(kin,NNinds));
        YbyKinaseSVD(kin,2:end) = cumsum(LIK.indivhits_red(kin,SVDinds)) - cumsum(LIK.indivnonhits_red(kin,SVDinds));

        % count the hits that make it into the top X%
        Lthis.hitsintop(kin) = sum(LIK.indivhits_red(kin,NNinds(1:Lthis.topind)));

        % count nonhits that we get to retrieve y% of the true hits
        zz = find(cumsum(LIK.indivhits_red(kin,NNinds))>=Lthis.hitind(kin),1); % index of all ligands that returns the desired number of hits
        Lthis.FPfortophits(kin) = zz - Lthis.hitind(kin);

    end

elseif REDO_PCA

    for kin=1:NLK

        % indices of active ligands sorted by predicted affinities
        % for REDO_PCA map these back into the full set of ligands
        NNinds  = Lthis.TrainLig(LpcaS.predsort(kin,:));
        SVDinds = Lthis.TrainLig(LpcaS.predsortSVD(kin,:));

        TPRbyKinase(kin,2:end)      = cumsum(LIK.indivhits(kin,NNinds)) / Lthis.FilteredInteractions(kin) ;
        FPRbyKinase(kin,2:end)      = cumsum(LIK.indivnonhits(kin,NNinds)) / Lthis.FilteredNonHits(kin);

        TPRbyKinaseSVD(kin,2:end)   = cumsum(LIK.indivhits(kin,SVDinds)) / Lthis.FilteredInteractions(kin) ;
        FPRbyKinaseSVD(kin,2:end)   = cumsum(LIK.indivnonhits(kin,SVDinds)) / Lthis.FilteredNonHits(kin);

        YbyKinase(kin,2:end) = cumsum(LIK.indivhits(kin,NNinds)) - cumsum(LIK.indivnonhits(kin,NNinds));
        YbyKinaseSVD(kin,2:end) = cumsum(LIK.indivhits(kin,SVDinds)) - cumsum(LIK.indivnonhits(kin,SVDinds));

        % Recall: count the hits that make it into the top X%
        Lthis.hitsintop(kin) = sum(LIK.indivhits(kin,NNinds(1:Lthis.topind)));

        % count nonhits that we get to retrieve y% of the true hits
        zz = find(cumsum(LIK.indivhits(kin,NNinds))>=Lthis.hitind(kin),1); % index of all ligands that returns the desired number of hits
        Lthis.FPfortophits(kin) = zz - Lthis.hitind(kin);

        for imethod = 1:4
            % indices of ligands sorted by predicted affinity
            % mapped into the full set of ligands (1:NLK)
            Ginds = Lthis.TrainLig(RSA{imethod}.predsort(kin,:));
            % TPR and FPR
            RSA{imethod}.TPRbyKinase(kin,2:end) = cumsum(LIK.indivhits(kin,Ginds))    / Lthis.FilteredInteractions(kin);
            RSA{imethod}.FPRbyKinase(kin,2:end) = cumsum(LIK.indivnonhits(kin,Ginds)) / Lthis.FilteredNonHits(kin);
            % the ROC curve is TPR (y axis) vs. FPR (x axis)

            % Youden is TP - FP (count, not rate)
            RSA{imethod}.YbyKinase(kin,2:end)   = cumsum(LIK.indivhits(kin,Ginds))  - cumsum(LIK.indivnonhits(kin,Ginds));
            
            % count hits in top X% of predicted affinities
            RSA{imethod}.hitsintop(kin) = sum(LIK.indivhits(kin,Ginds(1:Lthis.topind)));

            % nonhits to recover Y% of hits
            zz = find(cumsum(LIK.indivhits(kin,Ginds))>=Lthis.hitind(kin),1); % # of all ligands to get the hits we want
            RSA{imethod}.FPfortophits(kin) = zz - Lthis.hitind(kin); % how many nonhits among these

            % Area under the ROC curve 
            RSA{imethod}.AUC = 0.5*sum((RSA{imethod}.FPRbyKinase(:,2:end) - RSA{imethod}.FPRbyKinase(:,1:end-1)).* (RSA{imethod}.TPRbyKinase(:,1:end-1)+RSA{imethod}.TPRbyKinase(:,2:end)),2);

        end

    end

    % operations that did not need a loop over kinases
    for imethod=1:4
        RSA{imethod}.RecallTopX = RSA{imethod}.hitsintop ./ Lthis.FilteredInteractions;
        RSA{imethod}.PrecisionTopX = RSA{imethod}.hitsintop / Lthis.topind;% #hits above the cut / total above the cut
        RSA{imethod}.PrecisionTopY = Lthis.hitind ./ (Lthis.hitind + RSA{imethod}.FPfortophits);
    end

else % default case - neither REDO_PCA nor LASSO_MODE

    for kin=1:NLK
        TPRbyKinase(kin,2:end)      = cumsum(LIK.indivhits(kin,Lthis.predsort(kin,:))) / LIK.interactions(kin) ;
        FPRbyKinase(kin,2:end)      = cumsum(LIK.indivnonhits(kin,Lthis.predsort(kin,:))) / LIK.nonhits(kin);

        TPRbyKinaseSVD(kin,2:end)   = cumsum(LIK.indivhits(kin,Lthis.predsortSVD(kin,:))) / LIK.interactions(kin) ;
        FPRbyKinaseSVD(kin,2:end)   = cumsum(LIK.indivnonhits(kin,Lthis.predsortSVD(kin,:))) / LIK.nonhits(kin);

        YbyKinase(kin,2:end) = cumsum(LIK.indivhits(kin,Lthis.predsort(kin,:))) - cumsum(LIK.indivnonhits(kin,Lthis.predsort(kin,:)));
        YbyKinaseSVD(kin,2:end) = cumsum(LIK.indivhits(kin,Lthis.predsortSVD(kin,:))) - cumsum(LIK.indivnonhits(kin,Lthis.predsortSVD(kin,:)));

        % Recall: count the hits that make it into the top X% 
        Lthis.hitsintop(kin) = sum(LIK.indivhits(kin,Lthis.predsort(kin,1:Lthis.topind)));

        % count nonhits that we get to retrieve y% of the true hits
        zz = find(cumsum(LIK.indivhits(kin,:))>=Lthis.hitind(kin),1); % index of all ligands that returns the desired number of hits
        Lthis.FPfortophits(kin) = zz - Lthis.hitind(kin);

    end

end

% compute Area Under the Curve by the trapezoid rule
LIK.AUC_NN = 0.5*sum((FPRbyKinase(:,2:end) - FPRbyKinase(:,1:end-1)).* (TPRbyKinase(:,1:end-1)+TPRbyKinase(:,2:end)),2);
LIK.AUC_SVD = 0.5*sum((FPRbyKinaseSVD(:,2:end) - FPRbyKinaseSVD(:,1:end-1)).* (TPRbyKinaseSVD(:,1:end-1)+TPRbyKinaseSVD(:,2:end)),2);



%% identify individual optimal cutoffs and max Y using YbyKinase etc.

[LIKindiv.MaxY, MyMaxInd ] = max(YbyKinase,[],2);% max value and its index

if LASSO_MODE
    LIKindiv.BestCut = diag(RedLpca.sort_checkNN(:,MyMaxInd));% corresponding cutoff
elseif REDO_PCA
    LIKindiv.BestCut = diag(LpcaS.sort_checkNN(:,MyMaxInd));% corresponding cutoff
else
    LIKindiv.BestCut = diag(Lthis.sort_checkNN(:,MyMaxInd));% corresponding cutoff
end



LIKindiv.TP = diag(TPRbyKinase(:,MyMaxInd)).* LIK.interactions;% true positive count
LIKindiv.FP = diag(FPRbyKinase(:,MyMaxInd)).* LIK.nonhits;% false positive count

LIKindiv.indivTP = (Lthis.checkNN > repmat(LIKindiv.BestCut,1,NLC)).*LIK.indivhits;
LIKindiv.indivFP = (Lthis.checkNN > repmat(LIKindiv.BestCut,1,NLC)).*LIK.indivnonhits;

LIKindiv.TPtoo = sum(LIKindiv.indivTP,2);
LIKindiv.FPtoo = sum(LIKindiv.indivFP,2);

LIKindiv.TPN = LIKindiv.TP ./ LIK.interactions;% true positive rate
LIKindiv.FPN = LIKindiv.FP ./ LIK.interactions;
% alternative way to get TP, FP


% note -- can't use the second output of sort() as a *matrix* of indices
%TPRtoo = [zeros(NLK,1),cumsum(LIK.indivhits(Lpca.predsort),2)./repmat(LIK.interactions,1,NLC)];


%% now do this for all interactions of all kinases

if ~LASSO_MODE & ~REDO_PCA
    LIKall.preds = reshape(Lthis.checkNN,1,[]); % NN predicted values
    LIKall.preds_SVD = reshape(Lthis.checkSVD,1,[]); % SVD predicted values
    LIKall.hits = reshape(LIK.indivhits,1,[]); % hit  flag 1 if hit, 0 if not
    LIKall.nonhits = reshape(LIK.indivnonhits,1,[]); % nonhit flag, same idea
elseif LASSO_MODE
    LIKall.preds = reshape(RedLpca.checkNN,1,[]); % predicted values
    LIKall.preds_SVD = reshape(RedLpca.check,1,[]); % predicted values
    LIKall.hits = reshape(LIK.indivhits_red,1,[]); % hit  flag 1 if hit, 0 if not
    LIKall.nonhits = reshape(LIK.indivnonhits_red,1,[]); % nonhit flag, same idea
elseif REDO_PCA
    % Lthis.checkNN is full size but the inactive ligand columns are
    % identically zero
    LIKall.preds = reshape(Lthis.checkNN(:,Lthis.ActiveLigands),1,[]); % NN predicted values
    LIKall.preds_SVD = reshape(Lthis.checkSVD(:,Lthis.ActiveLigands),1,[]); % SVD predicted values
    LIKall.hits = reshape(LIK.indivhits(:,Lthis.ActiveLigands),1,[]); % hit  flag 1 if hit, 0 if not
    LIKall.nonhits = reshape(LIK.indivnonhits(:,Lthis.ActiveLigands),1,[]); % nonhit flag, same idea
end

LIKall.hitcount = sum(LIKall.hits);
LIKall.nonhitcount = sum(LIKall.nonhits);

% do this separately for the training group
% and separately for the rest (val + test group)
if ~LASSO_MODE & ~REDO_PCA
    LIKtrain.hits = reshape(Linds{1,1},1,[]); % hit  flag 1 if hit, 0 if not
    LIKtrain.nonhits = reshape(Linds{1,2},1,[]); % nonhit flag, same idea
    LIKrest.hits = reshape(Linds{2,1},1,[]); % hit  flag 1 if hit, 0 if not
    LIKrest.nonhits = reshape(Linds{2,2},1,[]); % nonhit flag, same idea
elseif LASSO_MODE
    LIKtrain.hits = reshape(Linds{1,1}(:,SLG),1,[]); % hit  flag 1 if hit, 0 if not
    LIKtrain.nonhits = reshape(Linds{1,2}(:,SLG),1,[]); % nonhit flag, same idea
    LIKrest.hits = reshape(Linds{2,1}(:,SLG),1,[]); % hit  flag 1 if hit, 0 if not
    LIKrest.nonhits = reshape(Linds{2,2}(:,SLG),1,[]); % nonhit flag, same idea
elseif REDO_PCA
    LIKtrain.hits = reshape(Linds{1,1}(:,Lthis.ActiveLigands),1,[]); % hit  flag 1 if hit, 0 if not
    LIKtrain.nonhits = reshape(Linds{1,2}(:,Lthis.ActiveLigands),1,[]); % nonhit flag, same idea
    LIKrest.hits = reshape(Linds{2,1}(:,Lthis.ActiveLigands),1,[]); % hit  flag 1 if hit, 0 if not
    LIKrest.nonhits = reshape(Linds{2,2}(:,Lthis.ActiveLigands),1,[]); % nonhit flag, same idea
end
LIKtrain.hitcount = sum(LIKtrain.hits);
LIKtrain.nonhitcount = sum(LIKtrain.nonhits);
LIKrest.hitcount = sum(LIKrest.hits);
LIKrest.nonhitcount = sum(LIKrest.nonhits);

% sort by predicted values
[LIKall.sortvals, LIKall.sortpos] = sort(LIKall.preds,'descend');
[LIKall.sortvals_SVD, LIKall.sortpos_SVD] = sort(LIKall.preds_SVD,'descend');
% 
% if ~LASSO_MODE
%     CumTPR_ini = zeros(1,NLK*NLC+1);
% else
%     CumTPR_ini = zeros(1,NLK*SNLC+1);
% end

for imethod=1:4
    for iset=1:3
        % TPR /FPR for the ROC curves
        % initialize
        RSA{imethod}.TPRbySet{iset} = zeros(1,NLK*Lthis.ActiveLigandCount+1);
        RSA{imethod}.FPRbySet{iset} = zeros(1,NLK*Lthis.ActiveLigandCount+1);
        % cumulative sums of hits and nonhits
        RSA{imethod}.TPRbySet{iset}(2:end) = cumsum(RG.HitsbySet{iset}(RSA{imethod}.allsort_pos)) / RG.HitCountbySet(iset);
        RSA{imethod}.FPRbySet{iset}(2:end) = cumsum(RG.NonHitsbySet{iset}(RSA{imethod}.allsort_pos)) / RG.NonHitCountbySet(iset);
        % area under the ROC curve computed by trapezoid rule:
        % sum( (x_(k+1) - x_k ) * (y_(k+1) _ y(k) ) / 2 )
        RSA{imethod}.ROCbySet(iset) = 0.5 * sum(...
            (RSA{imethod}.FPRbySet{iset}(2:end) - RSA{imethod}.FPRbySet{iset}(1:end-1)) .* ...
            (RSA{imethod}.TPRbySet{iset}(2:end) + RSA{imethod}.TPRbySet{iset}(1:end-1)) );
        %ROC_ALL = 0.5*sum((CumFPR_all(2:end) - CumFPR_all(1:end-1)).* (CumTPR_all(1:end-1)+CumTPR_all(2:end)),2);
    end
end

CumTPR_ini = zeros(1,NLK*Lthis.ActiveLigandCount+1);

CumTPR_all = CumTPR_ini;
CumFPR_all = CumTPR_ini;
CumTPR_all_SVD = CumTPR_ini;
CumFPR_all_SVD = CumTPR_ini;

CumTPR_all(2:end) = cumsum(LIKall.hits(LIKall.sortpos)) / LIKall.hitcount;
CumFPR_all(2:end) = cumsum(LIKall.nonhits(LIKall.sortpos)) / LIKall.nonhitcount;
CumTPR_all_SVD(2:end) = cumsum(LIKall.hits(LIKall.sortpos_SVD)) / LIKall.hitcount;
CumFPR_all_SVD(2:end) = cumsum(LIKall.nonhits(LIKall.sortpos_SVD)) / LIKall.nonhitcount;

% training kinase set only ..
CumTPR_train = CumTPR_ini;
CumFPR_train = CumTPR_ini;
CumTPR_train_SVD = CumTPR_ini;
CumFPR_train_SVD = CumTPR_ini;

CumTPR_train(2:end) = cumsum(LIKtrain.hits(LIKall.sortpos)) / LIKtrain.hitcount;
CumFPR_train(2:end) = cumsum(LIKtrain.nonhits(LIKall.sortpos)) / LIKtrain.nonhitcount;
CumTPR_train_SVD(2:end) = cumsum(LIKtrain.hits(LIKall.sortpos_SVD)) / LIKtrain.hitcount;
CumFPR_train_SVD(2:end) = cumsum(LIKtrain.nonhits(LIKall.sortpos_SVD)) / LIKtrain.nonhitcount;

% not training kinase set only ..
CumTPR_rest = CumTPR_ini;
CumFPR_rest = CumTPR_ini;
CumTPR_rest_SVD = CumTPR_ini;
CumFPR_rest_SVD = CumTPR_ini;

CumTPR_rest(2:end) = cumsum(LIKrest.hits(LIKall.sortpos)) / LIKrest.hitcount;
CumFPR_rest(2:end) = cumsum(LIKrest.nonhits(LIKall.sortpos)) / LIKrest.nonhitcount;
CumTPR_rest_SVD(2:end) = cumsum(LIKrest.hits(LIKall.sortpos_SVD)) / LIKrest.hitcount;
CumFPR_rest_SVD(2:end) = cumsum(LIKrest.nonhits(LIKall.sortpos_SVD)) / LIKrest.nonhitcount;



ROC_ALL = 0.5*sum((CumFPR_all(2:end) - CumFPR_all(1:end-1)).* (CumTPR_all(1:end-1)+CumTPR_all(2:end)),2);

ROC_train = 0.5*sum((CumFPR_train(2:end) - CumFPR_train(1:end-1)).* (CumTPR_train(1:end-1)+CumTPR_train(2:end)),2);
ROC_rest  = 0.5*sum((CumFPR_rest(2:end)  - CumFPR_rest(1:end-1)).* (CumTPR_rest(1:end-1)+CumTPR_rest(2:end)),2);

% ... and for SVD
ROC_ALL_SVD = 0.5*sum((CumFPR_all_SVD(2:end) - CumFPR_all_SVD(1:end-1)).* (CumTPR_all_SVD(1:end-1)+CumTPR_all_SVD(2:end)),2);

ROC_train_SVD = 0.5*sum((CumFPR_train_SVD(2:end) - CumFPR_train_SVD(1:end-1)).* (CumTPR_train_SVD(1:end-1)+CumTPR_train_SVD(2:end)),2);
ROC_rest_SVD  = 0.5*sum((CumFPR_rest_SVD(2:end)  - CumFPR_rest_SVD(1:end-1)).* (CumTPR_rest_SVD(1:end-1)+CumTPR_rest_SVD(2:end)),2);

if ~REDO_PCA
    fprintf('ROC: all %.4f|%.4f train %.4f|%.4f rest %.4f|%.4f ; ',...
        ROC_ALL,ROC_ALL_SVD,ROC_train,ROC_train_SVD,ROC_rest,ROC_rest_SVD);
end

if LASSO_MODE
    LassoNote = sprintf(' -- LASSO:  %d/%d ligands [%.2f%%] | %d/%d interactions [%.2f%%] ',...
        SNLC,NLC,SNLC/NLC*100,...
        sum(LIK.interactions_red),sum(LIK.interactions),sum(LIK.interactions_red)/sum(LIK.interactions)*100);
end


%% Youden etc. - look at TP - FP counts / rates and choose optimal cutoff[s]


Youden = CumTPR_all(2:end)*LIKall.hitcount - CumFPR_all(2:end)*LIKall.nonhitcount;
YoudenN = CumTPR_all(2:end)- CumFPR_all(2:end);
Youden_SVD = CumTPR_all_SVD(2:end)*LIKall.hitcount - CumFPR_all_SVD(2:end)*LIKall.nonhitcount;
YoudenN_SVD = CumTPR_all_SVD(2:end)- CumFPR_all_SVD(2:end);

% best cutoff for the ensemble based on Youden (TPC - FPC)
MaxY = max(Youden);
MaxY_SVD = max(Youden_SVD);

MaxYInd = find(Youden==MaxY,1);
BestCut = LIKall.sortvals(MaxYInd);
MaxYInd_SVD = find(Youden_SVD==MaxY_SVD,1);
BestCut_SVD = LIKall.sortvals(MaxYInd_SVD);

MaxYN = max(YoudenN);
MaxYN_SVD = max(YoudenN_SVD);
MaxYIndN = find(YoudenN==MaxYN,1);
BestCutN = LIKall.sortvals(MaxYIndN);

%% identify and count TP FP etc. with a specific cutoff value

% looking for entries in Lthis.checkNN that are above a chosen cutoff, and are
% hits / nonhits

% TP / FP flags by kinase x ligand group
LCut.indivTP = (Lthis.checkNN > BestCut).* LIK.indivhits;
LCut.indivFP = (Lthis.checkNN > BestCut).* LIK.indivnonhits;
LCut.indivTP_SVD = (Lthis.checkSVD > BestCut_SVD).* LIK.indivhits;
LCut.indivFP_SVD = (Lthis.checkSVD > BestCut_SVD).* LIK.indivnonhits;

% TP / FP counts by kinase
LCut.TP = sum(LCut.indivTP,2);
LCut.FP = sum(LCut.indivFP,2);
LCut.TP_SVD = sum(LCut.indivTP_SVD,2);
LCut.FP_SVD = sum(LCut.indivFP_SVD,2);

% TP / FP normalized to the interaction count
LCut.TPN = LCut.TP ./ LIK.interactions;
LCut.FPN = LCut.FP ./ LIK.interactions;
LCut.TPN_SVD = LCut.TP_SVD ./ LIK.interactions;
LCut.FPN_SVD = LCut.FP_SVD ./ LIK.interactions;
%% choose kinases of interest for the various plots
% choose kinases to analyze
%interesting = [321,396,113,328,395];%147, 175, 116, 40, 184];%330, 135, 242,252, 356];%, 68 ,138 188,139];%, 212];%221, 370, 374];% 434,350,202];%85%396, 42
%interesting = [142, 200,276,214,396]; % for  alt_single_200_1000_500_200_40_65740.mat
%interesting = [10,31,73,115,144]; % test kinases with a strange AUC
%interesting = [350,202,423,394,239]; % test kinases with a strange AUC
%interesting = [349,44,196,165,431]; % test kinases with a strange AUC
%interesting = [350,202,434,281,208,29];%212,202,85,155]; % kinase selection for run 96396
%InterestingKinases = [423,191,76,200,350,202];%,434];%,281,208,29];%212,202,85,155]; % kinase selection for run 96396
%InterestingKinases = [336,396,29,221,413];% used for run 75004  -set2
InterestingKinases = [221,396,29,88,353];%,163,413];%[332,366,32,162,423,251];% used for run 74839 - set2
ione = InterestingKinases(1);
itwo = InterestingKinases(2);

ttsone = sprintf('Kinase no. %d [%d] %s',LigandInfoKinases(ione),ione,KGOT.Name{LigandInfoKinases(ione)});
ttstwo = sprintf('Kinase no. %d [%d] %s',LigandInfoKinases(itwo),ione,KGOT.Name{LigandInfoKinases(itwo)});

testlabel = {'train','val','test'};

%% Look at ranking how many: {TP in the top x% | FP above y% of hits}

% use: Lthis.sort_checkP -- each kinases pred values sorted descending
%   or LPredHits.all -- array with only the hit values at the resp. pos.

% old code - commented out, to be removed
% how many hits are in the top x% [...] of all predictions
% Lthis.topfraction = 0.05; % fraction x% we want
% Lthis.topind = ceil(Lthis.topfraction*Lthis.ActiveLigandCount); % index of that fraction
% 
% Lthis.hitsintop = zeros(NLK,1); % will hold the answer (# hits in top x%)
% 
% for ikin = 1:NLK
%     if LASSO_MODE      
%         Lthis.hitsintop(ikin) = sum(LIK.indivhits_red(ikin,RedLpca.predsort(ikin,1:Lthis.topind)));
%     elseif REDO_PCA
%         NNinds  = Lthis.TrainLig(LpcaS.predsort(kin,:)); % sorted active ligands mapped to the full set
%         Lthis.hitsintop(ikin) = sum(LIK.indivhits(ikin,NNinds(1:Lthis.topind)));
%     else
%         Lthis.hitsintop(ikin) = sum(LIK.indivhits(ikin,Lthis.predsort(ikin,1:Lthis.topind)));
%     end
% end
%
% how many false positives for y% of the true ones
% %Lthis.hitfraction = 0.5;% fraction y% we want
% %Lthis.hitind = ceil(Lthis.hitfraction * LIK.interactions); % corresp index (among hits, sorted descending)
% Lthis.hitind = ceil(Lthis.hitfraction * Lthis.FilteredInteractions); % corresp index (among hits, sorted descending)
% Lthis.FPfortophits = zeros(NLK,1); % will hold the answers (# FP above y% of hits)
% for ikin=1:NLK
%     if Lthis.hitind(ikin) > 0
%     % find the [lowest] value of the top y% of hits
%     zz = find(Lthis.PredHits.all(ikin,:),Lthis.hitind(ikin)); 
%     zval = Lthis.PredHits.all(ikin,zz(end));
%     % count the nonhits above this value
%     Lthis.FPfortophits(ikin) = length(find(Lthis.PredNonHits.all(ikin,:)>=zval));
%     else
%         % fprintf('No recorded interactions for %s \n',kinlabel(ikin) );
%     end
% end

% brief outputs
% fraction of hits that are in the top 5% [or 1%] of predicted affinities (means)
LK = {find(LIK.TVC==1),find(LIK.TVC>2)};% 1=train 2=test ONLY

if ~REDO_PCA
    fprintf('ROC (mean by kin/by hit) tr %.4f|%.4f test %.4f|%.4f ;',...
        mean(LIK.AUC_NN(LK{1}),"omitnan"), sum(LIK.AUC_NN(LK{1}).*Lthis.FilteredInteractions(LK{1}),'omitnan') / sum(Lthis.FilteredInteractions(LK{1}),'omitnan'),...
        mean(LIK.AUC_NN(LK{2}),"omitnan"), sum(LIK.AUC_NN(LK{2}).*Lthis.FilteredInteractions(LK{2}),'omitnan') / sum(Lthis.FilteredInteractions(LK{2}),'omitnan'));
end

% first value is averaged by kinase -- number of interactions in top /
% total interactions for that kinase 
% note that some kinases have no "filtered" interactions i.e. when their
% ligands are not active (not hit by the train kinase set)
RecallTopX = Lthis.hitsintop ./ Lthis.FilteredInteractions; % hits in top x% / all hits [TP + FN]
PrecisionTopX = Lthis.hitsintop / Lthis.topind; % hits in top x% / total top x% [TP + FP]
% second value is by interaction count ie all interactions in the top
% fraction (for all kinases) / the corresponding total
if ~REDO_PCA
    fprintf('[top %d%% of total] Recall(mean/global): train %.4f|%.4f test %.4f|%.4f ',...
        100*Lthis.topfraction,...
        mean(RecallTopX(LK{1}),"omitnan"),...
        sum(Lthis.hitsintop(LK{1}))/sum(Lthis.FilteredInteractions(LK{1})),...
        mean(RecallTopX(LK{2}),"omitnan"),...
        sum(Lthis.hitsintop(LK{2}))/sum(Lthis.FilteredInteractions(LK{2})));

    fprintf('Precision(mean/median): train %.4f|%.4f test %.4f|%.4f ;',...
        mean(PrecisionTopX(LK{1}),"omitnan"),median(PrecisionTopX(LK{1}),"omitnan"),...
        mean(PrecisionTopX(LK{2}),"omitnan"),median(PrecisionTopX(LK{2}),"omitnan"));
end
z2 = Lthis.hitind ./ (Lthis.hitind + Lthis.FPfortophits); % TP / (TP + FP) [aka Precision] to capture y% of hits for each kinase
% not clear what this does - trying to get a global equivalent to the FPR
% [z2sort, z2sortind] = sort(z2,'descend'); % sort the kinases by descending precision
% z2ind_train = find(cumsum(Lthis.FilteredInteractions(z2sortind(LK{1})))>Lthis.hitfraction*sum(Lthis.FilteredInteractions(LK{1})),1); % index for 50% of all hits
% z2ind_rest =  find(cumsum(Lthis.FilteredInteractions(z2sortind(LK{2})))>Lthis.hitfraction*sum(Lthis.FilteredInteractions(LK{2})),1); % index for 50% of all hits
% z2med_train = z2sort(z2ind_train); % corresponding FP rate
% z2med_rest = z2sort(z2ind_rest);
% 
% fprintf('[top %d%% of hits] Prec (TP/(TP+FP)): train %.4f|%.4f rest %.4f|%.4f ;', Lthis.hitfraction*100,...
%     median(z2(LK{1}),'omitnan'), z2med_train,...
%     median(z2(LK{2}),'omitnan'), z2med_rest);

if ~REDO_PCA
    fprintf('[top %d%% of hits] Precision (mean/median): train %.4f|%.4f test %.4f|%.4f ;', Lthis.hitfraction*100,...
        mean(z2(LK{1}),'omitnan'), median(z2(LK{1}),'omitnan'),...
        mean(z2(LK{2}),'omitnan'), median(z2(LK{2}),'omitnan'));

    fprintf(' net:%s |%s\n',netstring,RunNoteString); 
end

%% line output [to be] used in ResultSummaries
if REDO_PCA
    % output using the RSA and RG constructs
    for im=1:4
        fprintf('%12s ',RSA{im}.label);
        fprintf('ROC (gl) trn %.4f val %.4f tst %.4f ', RSA{im}.ROCbySet(1:3));
        if exist('ff','var')
            fprintf(ff,'%12s ',RSA{im}.label);
            fprintf(ff, ' %.4f %.4f %.4f ', RSA{im}.ROCbySet(1:3));
        end
        % ROC by kinase
        fprintf('(mn) by kin|hit ');
        setlabel ={'trn','val','tst'};
        for is=1:3
            inds = find(LIK.TVC==is);
            inds = find(LIK.TVC==is & Lthis.FilteredInteractions > 0);
            fprintf('%s %.4f|%.4f ',...
                setlabel{is},...
                mean(RSA{im}.AUC(inds),'omitnan'),...
                sum(RSA{im}.AUC(inds).*Lthis.FilteredInteractions(inds),'omitnan')/sum(Lthis.FilteredInteractions(inds)));
            if exist('ff','var')
                fprintf(ff, '%.4f %.4f ',...
                    mean(RSA{im}.AUC(inds),'omitnan'),...
                    sum(RSA{im}.AUC(inds).*Lthis.FilteredInteractions(inds),'omitnan')/sum(Lthis.FilteredInteractions(inds)));
            end
        end
        % Precision / Recall
        fprintf('[top %d%% of all] Rec(mn/gl):',100*Lthis.topfraction);
        for is=1:3
            %inds = find(LIK.TVC==is);
            inds = find(LIK.TVC==is & Lthis.FilteredInteractions > 0);
            fprintf('%s %.4f|%.4f ',...
                setlabel{is},...
                mean(RSA{im}.RecallTopX(inds),'omitnan'),...
                sum(RSA{im}.hitsintop(inds),'omitnan')/sum(Lthis.FilteredInteractions(inds)));
            if exist('ff','var')
                fprintf(ff, '%.4f %.4f ',...
                    mean(RSA{im}.RecallTopX(inds),'omitnan'),...
                    sum(RSA{im}.hitsintop(inds),'omitnan')/sum(Lthis.FilteredInteractions(inds)));
            end
        end        
        fprintf(' Prc(mn/mdn):');
        for is=1:3
            inds = find(LIK.TVC==is);
            inds = find(LIK.TVC==is & Lthis.FilteredInteractions > 0);
            fprintf('%s %.4f|%.4f ',...
                setlabel{is},...
                mean(RSA{im}.PrecisionTopX(inds),'omitnan'),...
                median(RSA{im}.PrecisionTopX(inds),'omitnan'));
            if exist('ff','var')
                fprintf(ff, '%.4f %.4f ',...
                    mean(RSA{im}.PrecisionTopX(inds),'omitnan'),...
                    median(RSA{im}.PrecisionTopX(inds),'omitnan'));
            end
        end
        fprintf('[top %d%% of hits]',100*Lthis.hitfraction)
        fprintf(' Prc(mn/mdn):');
        for is=1:3
            inds = find(LIK.TVC==is);
            inds = find(LIK.TVC==is & Lthis.FilteredInteractions > 0);
            fprintf('%s %.4f|%.4f ',...
                setlabel{is},...
                mean(RSA{im}.PrecisionTopY(inds),'omitnan'),...
                median(RSA{im}.PrecisionTopY(inds),'omitnan'));
            if exist('ff','var')
                fprintf(ff, '%.4f %.4f ',...
                    mean(RSA{im}.PrecisionTopY(inds),'omitnan'),...
                median(RSA{im}.PrecisionTopY(inds),'omitnan'));
            end
        end
        fprintf('Net: %s File: %s',netstring,mfn);
        fprintf('\n');
        if exist('ff','var')
            fprintf(ff, '%s %s \n',netstring,mfn);
        end
    end
end
%%

if NO_FIGURES, return; end

% *** PLOTS ***

%% Items for all plots

% Colors for train, val, test groups
AllColor = [0.3010 0.7450 0.9330]; % light blue
TrainColor = [0 0.4470 0.7410];% medium / dark blue
ValColor = [0.9290 0.6940 0.1250];% medium orange
TestColor = [0.8500 0.3250 0.0980];% medium red
TVColor = 0.5*(ValColor + TestColor);

TVCColors = [TrainColor; ValColor; TestColor];

% Utility function
kinlabel = @(i) sprintf('Kin.%d [%d] %s %s',...
    LigandInfoKinases(i),...
    i,...
    KGOT.Name{LigandInfoKinases(i)},...
    testlabel{LIK.TVC(i)});

% Labels for kinase sets and prediction methods
ShortMethodLabel = {'NN','SVD','Lin','Combo'};

%% figure 100: TPR , FPR , TPR - FPR (Youden) vs. cutoff
PLOT_THIS =  true;%false;%

if PLOT_THIS

    FigNo=100;
    figure(FigNo)
    clf
    % true & false positives as a function of the cutoff value
    % for a single kinase

    if ~LASSO_MODE & ~REDO_PCA
        plot(Lthis.sort_checkNN(ione,:),TPRbyKinase(ione,2:end)*LIK.interactions(ione)-FPRbyKinase(ione,2:end)*LIK.nonhits(ione),'LineWidth',2); hold on;
        plot(Lthis.sort_checkNN(ione,:),TPRbyKinaseSVD(ione,2:end)*LIK.interactions(ione)-FPRbyKinaseSVD(ione,2:end)*LIK.nonhits(ione),'LineWidth',2); hold on;
        plot(Lthis.sort_checkNN(ione,:),TPRbyKinase(ione,2:end)*LIK.interactions(ione),'LineWidth',2); hold on;
        plot(Lthis.sort_checkNN(ione,:),FPRbyKinase(ione,2:end)*LIK.nonhits(ione),'LineWidth',2); hold on;
        plot([0,1],LIK.interactions(ione)*[1,1],'k--');
        legend('Youden [NN+SVD]',' -- -- [SVD only]','true positive rate (TPR)','false positive rate (FPR)','number of hits')
    elseif LASSO_MODE
        plot(RedLpca.sort_checkNN(ione,:),TPRbyKinase(ione,2:end)); hold on;
        plot(RedLpca.sort_checkNN(ione,:),FPRbyKinase(ione,2:end)); hold on;
    elseif REDO_PCA
        Fig100_m1 = 3; % main set
        Fig100_m2 = 2; % comparison set
        
        plot(RSA{Fig100_m1}.sort_check(ione,:),RSA{Fig100_m1}.TPRbyKinase(ione,2:end)*Lthis.FilteredInteractions(ione)-RSA{Fig100_m1}.FPRbyKinase(ione,2:end)*Lthis.FilteredNonHits(ione),'LineWidth',2); hold on;
        plot(RSA{Fig100_m2}.sort_check(ione,:),RSA{Fig100_m2}.TPRbyKinase(ione,2:end)*Lthis.FilteredInteractions(ione)-RSA{Fig100_m2}.FPRbyKinase(ione,2:end)*Lthis.FilteredNonHits(ione),'LineWidth',2); hold on;
        plot(RSA{Fig100_m1}.sort_check(ione,:),RSA{Fig100_m1}.TPRbyKinase(ione,2:end)*Lthis.FilteredInteractions(ione),'LineWidth',2); hold on;
        plot(RSA{Fig100_m1}.sort_check(ione,:),RSA{Fig100_m1}.FPRbyKinase(ione,2:end)*Lthis.FilteredNonHits(ione),'LineWidth',2); hold on;
        plot([0,1],Lthis.FilteredInteractions(ione)*[1,1],'k--');
        legend( ...
            ['Youden (',RSA{Fig100_m1}.label,')'],...
            [' -- -- ',RSA{Fig100_m2}.label],...
            'true positive rate (TPR)',...
            'false positive rate (FPR)',...
            sprintf('number of hits (of %d active ligands)',Lthis.ActiveLigandCount));
    end

    xlabel 'Cutoff'
    ylabel 'Count'
    title 'TP, FP and Youden statistic for one kinase'
    subtitle(kinlabel(ione),'interpreter','none')

    if REDO_PCA
        xlim([min(RSA{Fig100_m1}.check(ione,LIK.indivhits(ione,Lthis.ActiveLigands)>0)), max(RSA{1}.check(ione,LIK.indivhits(ione,Lthis.ActiveLigands)>0))]+[-0.05, 0.05]);
        ylim([-0.05, 1.15]*Lthis.FilteredInteractions(ione))
    else
        xlim([min(Lthis.checkNN(ione,LIK.indivhits(ione,:)>0)), max(Lthis.checkNN(ione,LIK.indivhits(ione,:)>0))]+[-0.05, 0.05]);
        ylim([-0.05, 1.15]*LIK.interactions(ione))
    end

    FormatThisFigure
    if SAVEFIGURES
        FigureName = [SaveDir, 'Fig100_YoudenOneKinase_', netfileID];
        SaveThisFigure5x5;
    end

    clear Fig100_m1 Fig100_m2 FigNo
end



%% figure 101: TPR vs FPR (ROC curve[s] by kinase)
PLOT_THIS =  true;%false;%
if PLOT_THIS

    FigNo = 101;
    Fig101_m = 1;
    figure(FigNo)
    clf

    USE_LOGSCALE = false;

    ShortMethodLabel = {'NN','SVD','Lin','Combo'};

    KinSuffix = '';

    % ROC curve: normalized hit and vs. counts
    % up to given predicted affinity values, sorted by affinity value
    legstr={};%build the legend as we add lines
    for j=1:length(InterestingKinases)%5%100
        myK = InterestingKinases(j);
        if REDO_PCA
            plot(RSA{Fig101_m}.FPRbyKinase(myK,:),RSA{Fig101_m}.TPRbyKinase(myK,:),'-','LineWidth',2);hold on;
        else
            plot(FPRbyKinase(myK,:),TPRbyKinase(myK,:),'-','LineWidth',2);hold on;
        end

        legstr = [legstr,sprintf('%s %.3f',kinlabel(myK) ,RSA{Fig101_m}.AUC(myK))];

        KinSuffix = [KinSuffix,sprintf('_%d',myK)];
    end
    xlabel 'FP Rate'
    ylabel 'TP Rate'
    title 'ROC Curve'
    %subtitle(['selected kinases ',RSA{Fig101_m}.label]);
    subtitle(['selected kinases ',ShortMethodLabel{Fig101_m}],'interpreter','none');
    if USE_LOGSCALE
        set(gca,'xscale','log');
        set(gca,'yscale','log');
    else
        xlim([-0.1 1.1]);
        ylim([-0.1 1.1]);
    end
    legend(legstr,'interpreter','none');
    legend('location','southeast')
    FormatThisFigure
    if SAVEFIGURES
        FigureName = [SaveDir, 'Fig101_ROCbyKinase_',ShortMethodLabel{Fig101_m}, netfileID, KinSuffix];
        FigW=6.5;
        FigH=6.5;
        SaveThisFigure;
    end

    clear FigNo Fig101_m
end
%% figure 102: histograms of area under ROC curve for individual kinases
PLOT_THIS = true;
if PLOT_THIS
    figure(102); clf;
    BAR_PLOT = true;%false;%
    if REDO_PCA
        Fig102_m=3;
        This = RSA{Fig102_m}.AUC;
    else
        This = LIK.AUC_NN ;
    end


    gridd = 0:0.05:1;

    NormType = 'percent';%'count';%

    if BAR_PLOT

        SHOW3 = true;%false;% % show train, val, test separately, or just train & (val+test)

        htrain = histcounts(This(LIK.TVC==1),gridd,'normalization',NormType);
        hval = histcounts(This(LIK.TVC==2),gridd,'normalization',NormType);
        htest = histcounts(This(LIK.TVC==3),gridd,'normalization',NormType);

        hvt = histcounts(This(LIK.TVC>1),gridd,'normalization',NormType);% val+test
        htot = histcounts(This(LIK.TVC>0),gridd,'normalization',NormType);% train+val+test

        gridc = 0.5*(gridd(1:end-1)+gridd(2:end));%centered grid points to match the histcount function

        if SHOW3

            % show the total as stairs
            h1=histogram('binedges',gridd,'bincounts',htot,'displaystyle','stairs','linewidth',2,'EdgeColor','k');hold on;


            b2=bar(gridc,[htest;hval;htrain]);hold on;

            b2(1).FaceColor=TestColor;%[0.9290 0.6940 0.1250];
            b2(2).FaceColor=ValColor;%[0.8500 0.3250 0.0980];
            b2(3).FaceColor=TrainColor;%[0.3010 0.7450 0.9330];

            legend('all','test','validation','train')
        else

            % show the total as  bars
            h1 = histogram('binedges',gridd,'bincounts',htot,...
                'displaystyle','bar',...
                'linewidth',1,...
                'EdgeColor','k',...
                'FaceColor',AllColor);hold on;

            b2=bar(gridc,[hvt;htrain]);hold on;
            b2(1).FaceColor=TVColor;%[0.9290 0.6940 0.1250];
            b2(2).FaceColor=TrainColor;%[0.8500 0.3250 0.0980];

            legend('all','test + validate','train')
        end
    else

        DisplayStyle = 'bar';
        htot = histogram(This,gridd,'normalization',NormType,'displaystyle',DisplayStyle); hold on;
        hvt = histogram(This(LIK.TVC>1),gridd,'normalization',NormType,'displaystyle','stairs','LineWidth',2); hold on;
        htest = histogram(This(LIK.TVC==3),gridd,'normalization',NormType,'displaystyle','stairs','LineWidth',2); hold on;
        legend('all','val+test','test')
    end



    legend('location','northwest')
    xlabel 'Area under the ROC curve'
    ylabel(['Kinases (',NormType,')']);

    if ~REDO_PCA
        title('AUC Distribution by Kinase')
    else
        title(['AUC Distribution by Kinase (',RSA{Fig102_m}.label,')'] )
        ylim([0 60]);
    end
    subtitle(['net: ',netstring],'interpreter','none');

    FormatThisFigure
    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig102_AUC_histo', netfileID ];
        SaveThisFigure5x5

    end

    clear Fig102_m
end

%% figure 103: scatter plot of AUC and number of hits
PLOT_THIS = true;
USE_SVD = false;%true;%
FigNo = 2103;
if PLOT_THIS
    figure(FigNo);clf;

    This = LIK.AUC_NN;
    if USE_SVD
        This = LIK.AUC_SVD;
    end

    if REDO_PCA
        Fig103_m = 4;
        
        if USE_SVD
            Fig103_m=2;
        end
        This = RSA{Fig103_m}.AUC;
    end

    for itvc=1:3 % loop over train, val, test
        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        %s=scatter(LIK.interactions(LIK.TVC==itvc),This(LIK.TVC==itvc),...
        s=scatter(Lthis.FilteredInteractions(LIK.TVC==itvc),This(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',TVCColors(itvc,:),...
            'MarkerFaceAlpha',0.5);

        hold on
        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = row;
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';
        s.DataTipTemplate.DataTipRows(1).Label = 'Hits';
        s.DataTipTemplate.DataTipRows(2).Label = 'AUC';
        s.DataTipTemplate.DataTipRows(2).Format = '%.4f';

    end

    box on
    set(gca,'xscale','log')
    ylabel('AUC (area under ROC curve)');
    xlabel('Number of recorded interactions');
    legend('train','validation','test')
    tstring = 'ROC area vs. Interaction Count, by kinase ';
    if REDO_PCA
        tstring = {tstring, RSA{Fig103_m}.label};
    elseif USE_SVD 
        tstring = [tstring,' - SVD prediction']; 
    end
    title(tstring)
    subtitle(['net: ',netstring],'interpreter','none');
    ylim([0.4 1]);

    legend('location','southeast');

    FormatThisFigure;
    legend('boxon');

    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig103_AUCvsIntCount_fixedscale_', netfileID ];
        SaveThisFigure5x5
    end
end
%% figure 105: TPR vs FPR (global ROC curve)

PLOT_THIS = true;
if PLOT_THIS
    figure(105)
    clf
    Fig105_m = 1;

    USE_LOGSCALE = false;

    % ROC curve: normalized hit and vs. counts
    % up to given predicted affinity values, sorted by affinity value
    legstr = cell(1,3);

    if REDO_PCA

        setnames = {'train','val','test'};
        for iset=1:3
            plot(RSA{Fig105_m}.FPRbySet{iset}, RSA{Fig105_m}.TPRbySet{iset},'-','LineWidth',2,'Color',TVCColors(iset,:)); hold on;
            legstr{iset} = sprintf('%6s - AUC:%.3f',setnames{iset},RSA{Fig105_m}.ROCbySet(iset));
        end

        % OPTIONAL: add dashed lines for another method (SVD or LinReg)
        % legend('noupdate')
        % for iset=1:3
        %     plot(RSA{3}.FPRbySet{iset}, RSA{3}.TPRbySet{iset},'--','Color',TVCColors(iset,:)); hold on;
        %     legstr{iset} = sprintf('%6s - AUC:%.3f',setnames{iset},RSA{1}.ROCbySet(iset));
        % end

    else

        plot(CumFPR_all,CumTPR_all,'-','LineWidth',2,'Color',AllColor); hold on;
        legstr{1} = sprintf('all   - AUC:%.3f',ROC_ALL);

        plot(CumFPR_train,CumTPR_train,'-','LineWidth',2,'Color',TrainColor); hold on;
        legstr{2} = sprintf('train - AUC:%.3f',ROC_train);

        plot(CumFPR_rest,CumTPR_rest,'-','LineWidth',2,'Color',TVColor); hold on;
        legstr{3} = sprintf('val + test - AUC:%.3f',ROC_rest);

    end




    xlabel 'FP Rate'
    ylabel 'TP Rate'
    if USE_LOGSCALE
        set(gca,'xscale','log');
        set(gca,'yscale','log');
    else
        xlim([-0.1 1.1]);
        ylim([-0.1 1.1]);
    end
    title('ROC Curve')

    if REDO_PCA
        %title(['ROC Curve ',RSA{Fig105_m}.label]);
        title(['ROC Curve ',ShortMethodLabel{Fig105_m}]);
    else
        title('ROC Curve');
    end

    %subtitle(['net: ',netstring],'interpreter','none');
    legend(legstr,'interpreter','none');
    legend('location','southeast')
    FormatThisFigure

    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig105_ROC_global', netfileID ];
        FigW=6.5;
        FigH=6.5;
        SaveThisFigure

    end

end

%% figure 205: TPR vs FPR (ROC curve, SVD only)
PLOT_THIS = false;
if PLOT_THIS
    figure(205)
    clf

    USE_LOGSCALE = false;

    % ROC curve: normalized hit and vs. counts
    % up to given predicted affinity values, sorted by affinity value
    legstr = cell(1,3);

    plot(CumFPR_all_SVD,CumTPR_all_SVD,'-','LineWidth',2,'Color',AllColor); hold on;
    legstr{1} = sprintf('all   - AUC:%.2f',ROC_ALL_SVD);

    plot(CumFPR_train_SVD,CumTPR_train_SVD,'-','LineWidth',2,'Color',TrainColor); hold on;
    legstr{2} = sprintf('train - AUC:%.2f',ROC_train_SVD);

    plot(CumFPR_rest_SVD,CumTPR_rest_SVD,'-','LineWidth',2,'Color',TVColor); hold on;
    legstr{3} = sprintf('val + test - AUC:%.2f',ROC_rest_SVD);



    xlabel 'FP Rate (SVD only)'
    ylabel 'TP Rate (SVD only)'
    if USE_LOGSCALE
        set(gca,'xscale','log');
        set(gca,'yscale','log');
    else
        xlim([-0.1 1.1]);
        ylim([-0.1 1.1]);
    end
    title('ROC Curve (SVD)')
    subtitle(['net: ',netstring],'interpreter','none');
    legend(legstr,'interpreter','none');
    legend('location','southeast')
    FormatThisFigure

    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig205_ROC_global', netfileID ];
        SaveThisFigure

    end

end
%% figure 106: Youden statistic, global (true positive - false positive rate)
PLOT_THIS = true;%false;%
NORMALIZE_THIS = false;%true;

FIG106_LOGSCALE = false;%true;%

FIG106_SHOW_SVD = false;%true;

if PLOT_THIS

    figure(106);
    clf
    if NORMALIZE_THIS
        plot(LIKall.sortvals, max(YoudenN,0),'LineWidth',2); hold on;
        if FIG106_SHOW_SVD, plot(LIKall.sortvals, max(YoudenN_SVD,0),'LineWidth',2);  end
        plot(LIKall.sortvals,CumTPR_all_SVD(2:end),'LineWidth',2);
        plot(LIKall.sortvals,CumFPR_all_SVD(2:end),'LineWidth',2);
        ylabel 'Youden (TPR - FPR)'
        subtitle '(normalized)'
        legend('Youden (NN+SVD)','Youden (SVD only)','TP','FP')
    else
        plot(LIKall.sortvals, max(Youden,0),'LineWidth',2); hold on;
        if FIG106_SHOW_SVD,  plot(LIKall.sortvals, max(Youden_SVD,0),'LineWidth',2); end
        plot(LIKall.sortvals,CumTPR_all_SVD(2:end)*LIKall.hitcount,'LineWidth',2);
        plot(LIKall.sortvals,CumFPR_all_SVD(2:end)*LIKall.nonhitcount,'LineWidth',2);
        plot([0,1],[1,1]*LIKall.hitcount,'k--')
        ylabel 'True Pos - False Pos (COUNTS)'
        %subtitle 'counts'

        if FIG106_SHOW_SVD
            legend('Youden (NN+SVD)','Youden (SVD only)','TP','FP','total hits')
        else
            legend('Youden (NN+SVD)','TP','FP','total hits')
        end
        if FIG106_LOGSCALE
            set(gca,'yscale','log');
        else
            % set the scale to show just the max hits
            ylim([-0.05, 1.15]*LIKall.hitcount)
        end

        xlim([-0.2 1.2])
    end
    xlabel 'Cutoff'
    title 'True Positives - False Positives'

    FormatThisFigure;

    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig106_YoudenGlobal_', netfileID ];
        FigW=6.5;
        FigH=6.5;
        SaveThisFigure

    end

end
%ylim([0 LIKall.hitcount]);
%% figure 123: scatter plot of AUC and number of hits, NN vs SVD
PLOT_THIS = true;
FigNo = 3123;
Method1 = 2;
Method2 = 4;
if PLOT_THIS
    figure(FigNo);clf;

    for itvc=1:3

        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        % s=scatter(...
        %     LIK.AUC_SVD(LIK.TVC==itvc),...
        %     LIK.AUC_NN(LIK.TVC==itvc),...
        %     LIK.interactions(LIK.TVC==itvc),...
        %     'filled',...
        %     'MarkerEdgeColor','k', ...
        %     'MarkerEdgeAlpha',0.5,...
        %     'MarkerFaceAlpha',0.5);

        s=scatter(...
            RSA{Method1}.AUC(LIK.TVC==itvc),...
            RSA{Method2}.AUC(LIK.TVC==itvc),...
            LIK.interactions(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);


        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(2).Label = 'AUC (NN )';
        s.DataTipTemplate.DataTipRows(1).Label = 'AUC (SVD)';

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        s.DataTipTemplate.DataTipRows(1).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(2).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on
    %set(gca,'xscale','log')
    ylabel([RSA{Method2}.label,' predicted']);
    xlabel([RSA{Method1}.label,' predicted']);

    legend('train','validation','test')
    legend('autoupdate','off')




    sll = plot([0,1],[0,1],':');
    %xlim([min(LIK.AUC_SVD(LIK.TVC==itvc))-0.05, 1.05]);
    %ylim([min(LIK.AUC_NN(LIK.TVC==itvc))-0.05, 1.05]);
    xlim([min(RSA{Method1}.AUC(LIK.TVC==itvc))-0.05, 1.05]);
    ylim([min(RSA{Method2}.AUC(LIK.TVC==itvc))-0.05, 1.05]);
    xlim([-0.05 1.05]);
    ylim([-0.05 1.05]);
    %axis equal
    tstring = 'AUC (area under ROC curve) - NN vs. SVD';

    title(tstring)
    subtitle(['net: ',netstring],'interpreter','none');

    FormatThisFigure;
    legend('location','northwest')
    if SAVEFIGURES

        FigureName=[SaveDir, 'AUC_NNvsSVD_', netfileID ];
        SaveThisFigure
    end
end

%% figure 125: TP vs FP 
%  scatter plot of true and false positives by kinase
%  using a specific cutoff

PLOT_THIS = false;%true;%
FigNo = 125;
USE_LOGSCALE = false;%true;
if PLOT_THIS
    figure(FigNo);clf;

    ThisCut = BestCut; % cutoff value to be used
    for itvc=1:3

        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        s=scatter(...
            LCut.FPN(LIK.TVC==itvc),...
            LCut.TPN(LIK.TVC==itvc),...
            LIK.interactions(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);

        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(1).Label = 'FP';
        s.DataTipTemplate.DataTipRows(1).Value = LCut.FP(LIK.TVC==itvc);
        s.DataTipTemplate.DataTipRows(2).Label = 'TP';
        s.DataTipTemplate.DataTipRows(2).Value = LCut.TP(LIK.TVC==itvc);

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        s.DataTipTemplate.DataTipRows(1).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(2).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on

    xlabel('False Positives / All hits');
    ylabel('True Positives / All hits (=TPR)');
    legend('train','validation','test')
    legend('autoupdate','off')

    if ~USE_LOGSCALE

        sll = plot([0,1],[0,1],':');
        xlim([-0.05, 0.05]+[min(LCut.FPN), max(LCut.FPN)]);
        ylim([-0.05, 0.05]+[min(LCut.TPN), max(LCut.TPN)]);

    else
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([1.0e-3 3])
        ylim([1.0e-3 3])
        plot([1.0e-3 3],[1.0e-3 3],':');
    end
    %axis equal
    tstring = 'TP and FP as fraction of all hits';
    title(tstring)
    ststring = sprintf('Global cutoff: %.3f',BestCut);
    subtitle(ststring,'interpreter','none');

    FormatThisFigure;
    legend('location','northeast')
    if SAVEFIGURES

        FigureName=[SaveDir, 'TPFPCut_', netfileID ];
        SaveThisFigure
    end
end


%% figure 225: TP vs FP , SVD only
%  scatter plot of true and false positives by kinase
%  using a specific cutoff

PLOT_THIS = false;%true;%
FigNo = 225;
USE_LOGSCALE = false;%true;
if PLOT_THIS
    figure(FigNo);clf;

    for itvc=1:3

        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        s=scatter(...
            LCut.FPN_SVD(LIK.TVC==itvc),...
            LCut.TPN_SVD(LIK.TVC==itvc),...
            LIK.interactions(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);

        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(1).Label = 'FP';
        s.DataTipTemplate.DataTipRows(1).Value = LCut.FP_SVD(LIK.TVC==itvc);
        s.DataTipTemplate.DataTipRows(2).Label = 'TP';
        s.DataTipTemplate.DataTipRows(2).Value = LCut.TP_SVD(LIK.TVC==itvc);

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        s.DataTipTemplate.DataTipRows(1).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(2).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on

    xlabel('False Positives / All hits');
    ylabel('True Positives / All hits (=TPR)');
    legend('train','validation','test')
    legend('autoupdate','off')


    if ~USE_LOGSCALE

        sll = plot([0,1],[0,1],':');
        xlim([-0.05, 0.05]+[min(LCut.FPN), max(LCut.FPN)]);
        ylim([-0.05, 0.05]+[min(LCut.TPN), max(LCut.TPN)]);

    else
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([1.0e-3 3])
        ylim([1.0e-3 3])
        plot([1.0e-3 3],[1.0e-3 3],':');
    end
    %axis equal
    tstring = 'TP and FP as fraction of all hits using SVD only';
    title(tstring)
    ststring = sprintf('Global cutoff: %.3f',BestCut_SVD);
    subtitle(ststring,'interpreter','none');

    FormatThisFigure;
    legend('location','northeast')
    if SAVEFIGURES

        FigureName=[SaveDir, 'TPFPCut_SVD', netfileID ];
        SaveThisFigure
    end
end

%% figure 310: (also 1310) counts of hits in the top x% vs. FP in the top y% of hits

PLOT_THIS = true;%false;%

FIG310_WEIGHT_BY_KINASE = true;%false;%
FigNo = 310 + 1000*(1-FIG310_WEIGHT_BY_KINASE);

SubtitleString = {'Weighted by Hits','Weighted by Kinase'};

if PLOT_THIS
    figure(FigNo);clf;
    USE_LOGSCALE = true;%false;%

    Lstring1 = sprintf('Recall (%% of all hits) in top %d%% predicted affinities',100*Lthis.topfraction);
    Lstring2 = sprintf('Precision (TP / (TP+FP) to capture %d%% of hits',100*Lthis.hitfraction);

    % big scatter plot
    subplot(3,3,[1,2,4,5]);

    for itvc=1:3 % loop over groups

        % NOTE: x axis is TP / (TP+FN) aka RECALL
        %       y axis is TP / (TP+FP) aka PRECISION 

        % omit kinases with no hits (the train / validation set might 
        % have these because the respective ligands do not interact 
        % with the train kinases)
        LKinaseIndices = find(LIK.TVC==itvc & Lthis.FilteredInteractions>0);
        %LKinaseIndices = find(LIK.TVC==itvc);

        LKI{itvc}=LKinaseIndices;
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        % %old (ligands not filtered )
        % s=scatter(...
        %     100*Lthis.hitsintop(LKinaseIndices)./LIK.interactions(LKinaseIndices),... % recall
        %     Lthis.hitind(LKinaseIndices) ./ (Lthis.hitind(LKinaseIndices) + Lthis.FPfortophits(LKinaseIndices)),...
        %     LIK.interactions(LKinaseIndices),...
        %     'filled',...
        %     'MarkerEdgeColor','k', ...
        %     'MarkerEdgeAlpha',0.5,...
        %     'MarkerFaceAlpha',0.5);

        % updated - use the computed recall / precision vectors
        s=scatter(...
            100*RSA{1}.RecallTopX(LKinaseIndices),... % recall
            RSA{1}.PrecisionTopY(LKinaseIndices),...
            Lthis.FilteredInteractions(LKinaseIndices),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);

        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(1).Label = sprintf('Hit count in top %d%% of preds:',100*Lthis.topfraction);
        s.DataTipTemplate.DataTipRows(1).Value = Lthis.hitsintop(LKinaseIndices);
        s.DataTipTemplate.DataTipRows(2).Label = sprintf('FP count above %d%% of hits:',100*Lthis.hitfraction);
        s.DataTipTemplate.DataTipRows(2).Value = Lthis.FPfortophits(LKinaseIndices);

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        %s.DataTipTemplate.DataTipRows(1).Format = '%.1f';
        %s.DataTipTemplate.DataTipRows(2).Format = '%.1f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on

    xlabel(Lstring1);
    ylabel(Lstring2);
    legend('train','validation','test')
    legend('autoupdate','off')
    legend('boxon')


    if USE_LOGSCALE
        set(gca,'yscale','log')
    end

    tstring = 'Hit Identification Efficiency';
    title(tstring)
    subtitle '';
    %subtitle(SubtitleString{1+FIG310_WEIGHT_BY_KINASE});


    FormatThisFigure;
    legend('location','southeast')
    legend boxon;

    % (1/2) small histograms
    subplot(3,3,7:8)
    
    
    if FIG310_WEIGHT_BY_KINASE
        histstring = '% of kinases';
    else
        histstring = '% of known interactions ("hits")';
    end


    %yyaxis right
    rrh = 0:10:100;
    rrhc = 0.5*(rrh(1:end-1)+rrh(2:end));
    
    % data to plot (1) -- % of hits that are in the top X% of predicted
    % affinities for each kinase (X% of predicted affinities for that
    % specfic kinase)
    ddx = 100*Lthis.hitsintop./LIK.interactions;
    ddx1 = ddx(LKI{1});
    ddx2 = ddx(union(LKI{2},LKI{3}));

    % update to use computed vectors and omit moot kinases (filtered
    % interactions)
    LKItotal = find(Lthis.FilteredInteractions);
    ddx = 100*RSA{1}.RecallTopX(LKItotal);
    ddx1 = 100*RSA{1}.RecallTopX(LKI{1});
    ddx2 = 100*RSA{1}.RecallTopX(union(LKI{2},LKI{3}));

    if FIG310_WEIGHT_BY_KINASE
        % histogram counts for all, train , test+val
        % hcx0 = histcounts(ddx,rrh);
        % hcx1 = histcounts(ddx(LKI{1}),rrh);
        % hcx2 = histcounts(ddx(union(LKI{2},LKI{3})),rrh);
        % updated
        hcx0 = histcounts(ddx,rrh);
        hcx1 = histcounts(ddx1,rrh);
        hcx2 = histcounts(ddx2,rrh);
    else
        % weighted histograms (each kinase x its known hit count)
        % hcx0 = whistcounts(ddx,LIK.interactions,rrh);
        % hcx1 = whistcounts(ddx(LKI{1}),LIK.interactions(LKI{1}),rrh);
        % hcx2 = whistcounts(ddx(union(LKI{2},LKI{3})),LIK.interactions(union(LKI{2},LKI{3})),rrh);
        % updated
        hcx0 = whistcounts(ddx,Lthis.FilteredInteractions(LKItotal),rrh);
        hcx1 = whistcounts(ddx1,Lthis.FilteredInteractions(LKI{1}),rrh);
        hcx2 = whistcounts(ddx2,Lthis.FilteredInteractions(union(LKI{2},LKI{3})),rrh);
    end
    % percentages if needed
    hcx1p = 100*hcx1 / sum(hcx1);
    hcx2p = 100*hcx2 / sum(hcx2);
    

    % plot the total as stairs (non-shaded)
    histogram('BinEdges',rrh,'BinCounts',hcx0,...
        'normalization','percentage',...
        'displaystyle','stairs',...
        'linewidth',2,...
        'edgecolor','k');
    hold on;

    

    %yyaxis left
    % plot the two subgroups as bars, also normalized to percents
    bp=bar(rrhc,[hcx1p;hcx2p]);
    bp(1).FaceColor=TrainColor;%[0.9290 0.6940 0.1250];
    bp(2).FaceColor=TVColor;%[0.8500 0.3250 0.0980];
    
    legend('all','train','val+test')
    
    legend('location','northwest')
    xlabel(Lstring1);
    ylabel(histstring);

    FormatThisFigure

    % two small histograms: 2/2
    subplot(3,3,[3,6])
    %FIG310b_WEIGHT_BY_KINASE = true;%false;%

    %yyaxis right
    rry = [kron([1e-4,1e-3,1e-2,0.1],[1,2,5]),1];% logbins
    %rry = [kron([1e-4,1e-3,1e-2,0.1],[1,10^(1/3),10^(2/3)]),1];% "equal"
    rryc = sqrt(rry(1:end-1).*rry(2:end)); % log bins, geometric mean

    % old
    ddy = Lthis.hitind ./ (Lthis.hitind + Lthis.FPfortophits);
    ddy1 = ddy(LKI{1});
    ddy2 = ddy(union(LKI{2},LKI{3}));

    %updated to use computed vectors and account for omitted ligands
    ddy = RSA{1}.PrecisionTopY(LKItotal);
    ddy1 = RSA{1}.PrecisionTopY(LKI{1});
    ddy2 = RSA{1}.PrecisionTopY(union(LKI{2},LKI{3}));

    if FIG310_WEIGHT_BY_KINASE
        hcy0 = histcounts(ddy,rry);
        hcy1 = histcounts(ddy1,rry);
        hcy2 = histcounts(ddy2,rry);
    else
        % old
        % hcy0 = whistcounts(ddy,LIK.interactions,rry);
        % hcy1 = whistcounts(ddy1,LIK.interactions(LKI{1}),rry);
        % hcy2 = whistcounts(ddy2,LIK.interactions(union(LKI{2},LKI{3})),rry);

        % updated
        hcy0 = whistcounts(ddy,Lthis.FilteredInteractions(LKItotal),rry);
        hcy1 = whistcounts(ddy1,Lthis.FilteredInteractions(LKI{1}),rry);
        hcy2 = whistcounts(ddy2,Lthis.FilteredInteractions(union(LKI{2},LKI{3})),rry);
    end

    hcy1p = 100*hcy1 / sum(hcy1);
    hcy2p = 100*hcy2 / sum(hcy2);

    %histogram(100*Lthis.hitsintop./LIK.interactions,rrh);
    %histogram(log10(ddy), ...  %'BinMethod','auto', ...
    histogram(...
        'BinCounts',hcy0, ...  %'BinMethod','auto', ...
        'Binedges',log10(rry),...
        'normalization','percentage',...
        'displaystyle','stairs',...
        'linewidth',2,...
        'edgecolor','k',...
        'orientation','horizontal');
    hold on;
    xlabel(histstring);

    %yyaxis left
    bp=barh(log10(rryc),[hcy2p;hcy1p]);
    bp(2).FaceColor=TrainColor;%[0.9290 0.6940 0.1250];
    bp(1).FaceColor=TVColor;%[0.8500 0.3250 0.0980];
    legend('all','val+test','train')

    legend('location','southeast')
    %set(gca,'yscale','log');
    ylabel(Lstring2);

    ylim(log10([rry(1) rry(end)]))
    %yticks(log10(rry));
    %yt = cell(length(rry)); for iy=1:length(rry), yt{iy} = sprintf('%.2e',rry(iy)); end
    yticks(-4:0);
    yt = cell(5); for iy=1:5, yt{iy} = sprintf('10^{%d}',iy-5); end
    %yminorticks(rry);
    yticklabels(yt)
    %ylabel 'Kinase count'

    FormatThisFigure

    if SAVEFIGURES

        FigureName=[SaveDir, sprintf('Fig%d_PrecisionRecall_Combo_', FigNo), netfileID ];
        SaveThisFigure

    end

end
%% figure 325: TP vs FP by kinase
%  scatter plot of true and false positives by kinase
%  using by-kinase cutoffs

PLOT_THIS = false;%true;%
FigNo = 325;
USE_LOGSCALE = false;%true;
if PLOT_THIS
    figure(FigNo);clf;

    for itvc=1:3

        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        s=scatter(...
            LIKindiv.FPN(LIK.TVC==itvc),...
            LIKindiv.TPN(LIK.TVC==itvc),...
            LIK.interactions(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);

        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(1).Label = 'FP';
        s.DataTipTemplate.DataTipRows(1).Value = LIKindiv.FP(LIK.TVC==itvc);
        s.DataTipTemplate.DataTipRows(2).Label = 'TP';
        s.DataTipTemplate.DataTipRows(2).Value = LIKindiv.TP(LIK.TVC==itvc);

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        s.DataTipTemplate.DataTipRows(1).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(2).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on

    xlabel('False Positives / All hits');
    ylabel('True Positives / All hits (=TPR)');
    legend('train','validation','test')
    legend('autoupdate','off')
    legend('boxon')



    if ~USE_LOGSCALE

        plot([0,1],[0,1],':');
        %xlim([-0.05, 0.05]+[min(LCut.FPN), max(LCut.FPN)]);
        %ylim([-0.05, 0.05]+[min(LCut.TPN), max(LCut.TPN)]);

    else
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([1.0e-3 3])
        ylim([1.0e-3 3])
        plot([1.0e-3 3],[1.0e-3 3],':');
    end
    %axis equal
    tstring = 'TP and FP as fraction of all hits';
    title(tstring)
    ststring = sprintf('by-kinase cutoffs');
    subtitle(ststring,'interpreter','none');

    FormatThisFigure;
    legend('location','southeast')
    legend boxon;
    if SAVEFIGURES

        FigureName=[SaveDir, 'TPFPCutbyKinase', netfileID ];
        SaveThisFigure
    end
end
%% figure 444 scatter plot of by kinase cutoffs vs best TP etc
PLOT_THIS = false;%true;%
FigNo = 444;
USE_LOGSCALE = false;%true;
if PLOT_THIS
    figure(FigNo);clf;

    for itvc=1:3

        LKinaseIndices = find(LIK.TVC==itvc);
        TrueKinaseIndices = LigandInfoKinases(LKinaseIndices);

        s=scatter(...
            LIKindiv.BestCut(LIK.TVC==itvc),...
            LIKindiv.TP(LIK.TVC==itvc)./LIK.interactions(LIK.TVC==itvc),...
            LIK.interactions(LIK.TVC==itvc),...
            'filled',...
            'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.5,...
            'MarkerFaceAlpha',0.5);

        hold on

        % data tips
        s.DataTipTemplate.DataTipRows(1).Label = 'Cut';
        s.DataTipTemplate.DataTipRows(1).Value = LIKindiv.BestCut(LIK.TVC==itvc);
        s.DataTipTemplate.DataTipRows(2).Label = 'TP';
        s.DataTipTemplate.DataTipRows(2).Value = LIKindiv.TP(LIK.TVC==itvc);

        row = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        row2 = dataTipTextRow('Ind',LKinaseIndices);
        s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('KinaseName',KGOT.Name(TrueKinaseIndices));
        s.DataTipTemplate.DataTipRows(end+1) = row2;
        s.DataTipTemplate.Interpreter = 'none';

        s.DataTipTemplate.DataTipRows(1).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(2).Format = '%.3f';
        s.DataTipTemplate.DataTipRows(3).Label = 'Hits';

    end

    box on

    xlabel('Cut threshold (by kinase)');
    ylabel('True Positives / All hits (=TPR)');
    legend('train','validation','test')
    legend('autoupdate','off')


    if ~USE_LOGSCALE

        sll = plot([0,1],[0,1],':');
        xlim([-0.05, 0.05]+[min(LCut.FPN), max(LCut.FPN)]);
        ylim([-0.05, 0.05]+[min(LCut.TPN), max(LCut.TPN)]);

    else
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([1.0e-3 3])
        ylim([1.0e-3 3])
        plot([1.0e-3 3],[1.0e-3 3],':');
    end
    %axis equal
    tstring = 'By-kinase cutoff values and TPR';
    title(tstring)


    FormatThisFigure;
    legend('location','northeast')
    if SAVEFIGURES

        FigureName=[SaveDir, 'TPFPCutbyKinase', netfileID ];
        SaveThisFigure
    end
end

%% figure 1346 -- mean and spread of predicted values by training status
PLOT_THIS = true;%false;%
if PLOT_THIS
    figure(1346);
    clf



    % kinase indices for each kinase in three passes (train, validate. test)
    TrueKinaseIndices = {tr.trainInd, tr.valInd, tr.testInd};
    % x indices for each kinase in three passes (train, validate. test)
    xinds = {LIK.unsort(tr.trainInd), LIK.unsort(tr.valInd), LIK.unsort(tr.testInd)};
    % color strings
    cstr = {'b','r','m'};
    cstr2 = {'c','g','g'};
    % marker strings
    mstr = {'o','d','s'};

    % formatting tweaks
    FS = [10,8,6,6]; % fonts used in FormatThisFigure
    MyLW = 0.6;
    SizeFac = 0.2;%0.4;
    offset = 0.05; % slight horizontal offset for the non-hit markers

    npl=5;
    titlestring = ['NN Predicted Affinities ' netstring];



    MyPredHits = Lthis.PredHits;
    MyPredNonHits = Lthis.PredNonHits;
    MyInteractions = Lthis.FilteredInteractions;

    for ipl=1:5 % loop over groups of ~100 kinases

        subplot(npl,1,ipl)

        % plot the hits with error bars and size indications
        for igr=1:3
            s=errorbar(xinds{igr},MyPredHits.mean(TrueKinaseIndices{igr}),MyPredHits.low(TrueKinaseIndices{igr}),MyPredHits.high(TrueKinaseIndices{igr}),[cstr{igr} '.'],'LineWidth',MyLW); hold on;
        end

        % the legend freezes after this
        if ipl==5
            legend('train','validation','test');
            %legend('location','southoutside');
            legend('orientation','horizontal');
            legend('autoupdate','off');
        end

        for igr=1:3
            s=scatter(xinds{igr},MyPredHits.mean(TrueKinaseIndices{igr}),SizeFac*(0.01+MyInteractions(TrueKinaseIndices{igr})),cstr{igr},mstr{igr},'filled'); hold on;
            row = dataTipTextRow('Ind',TrueKinaseIndices{igr});
            s.DataTipTemplate.DataTipRows(end+1) = row;
        end


        % non-hits
        for igr=1:3
            errorbar(xinds{igr}+offset,MyPredNonHits.mean(TrueKinaseIndices{igr}),MyPredNonHits.low(TrueKinaseIndices{igr}),MyPredNonHits.high(TrueKinaseIndices{igr}),[cstr2{igr} '.'],'LineWidth',MyLW,'Marker','none'); hold on;
            s=scatter(xinds{igr}+offset,MyPredNonHits.mean(TrueKinaseIndices{igr}),SizeFac*(0.01+MyInteractions(TrueKinaseIndices{igr})),cstr2{igr},mstr{igr}); hold on;
            row = dataTipTextRow('Ind',TrueKinaseIndices{igr});
            s.DataTipTemplate.DataTipRows(end+1) = row;
        end


        % lines for guidance
        plot([1,NLK],[0,0],'k:');
        plot([1,NLK],[1,1],'k:');
        ylim([-0.2 1.2])
        xlim([100*(ipl-1) 100*ipl])
        box on;

        FS = [10,8,8,8]; % fonts used in FormatThisFigure
        FormatThisFigure

        if(ipl==1)
            title(titlestring,'interpreter','none')
        end


        FormatThisFigure;

                if ipl<5
            legend('off');
        end
    end % loop over subplots (100 kinases per plot)



    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig1346_', netfileID ];

        if RAW_SVD
            FigureName = ['rawSVD_' FigureName];
        end

        SaveThisFigureLarge

    end
end
%% figure 1345 -- boxplot of hits and nonhits (alt to 1346)
PLOT_THIS = false;%true;%
if PLOT_THIS
    figure(1345)
    clf;
    zz1 = Lthis.checkNN;
    zz1(Linds{1,1}==0)=nan;
    zz2 = Lthis.checkNN;
    zz2(Linds{2,1}==0)=nan;
    subplot(2,1,1)
    boxplot(zz1');
    xlim([0 100])
    title 'Training set'
    subplot(2,1,2)
    boxplot(zz2');
    xlim([0 100])
    title 'Test/Validation set'
end

%% figure 3 - Paper Fig R2 - [predicted] affinities vs index 

PLOT_THIS = true;%false;%
if PLOT_THIS

    FigNo = 3;
    figure(FigNo); clf;
    USELOGSCALE_FIG3 = false;%true;%
    FIG3_SHOW_PCA = false;%
    FIG3_DETAIL = true;%true;%false;%

    rr = [500 1000];


    trainstatus={'train','validation','test'};

    for ifr=1:length(InterestingKinases) % loop over selected kinases

        % features of this kinase
        myK = InterestingKinases(ifr); % kinase index in the set of 455 ie.LigandInfoKinases
        myH = find(LCAMTred(myK,:)>=4); % hits for this kinase
        myNH = find(LCAMTred(myK,:)<4); % nonhits for this kinase

        myTVC = LIK.TVC(myK); % train status -- 1,2,3 / train,val,test

        % myVec
        clear 'myVec' 'normfacc' 'zz'

        myVec.ref = LCAMTs(myK,:)';
        myVec.SVD = Lthis.checkSVD(myK,:)';
        myVec.NN = Lthis.checkNN(myK,:)';

        % rescale so that the higher values are close to 1
        % choose the scale factor to set the mean of the top 1% of predicted afinities to 1 
        zz=sort(Lthis.checkNN(myK,:),'descend');
        normfacc = mean(zz(1:Lthis.hitind(myK))); % so the mean of the top fraction of hits is at 1.0
        myVec.NNrescale = Lthis.checkNN(myK,:)' ;%/ normfacc;

        clear 'zz' 

        %WW = {myVec.ref,myVec.SVD,myVec.NN};
        myVec = struct2table(myVec);

        % strings for labelling
        sKinaseNo = sprintf('Kin.%d [%d]', LigandInfoKinases(myK),myK);%
        sKinaseName = KGOT.Name{LigandInfoKinases(myK)};

        sTrainStatus = trainstatus{myTVC};
        sSetLabel = {'actual','SVD','NN+SVD'};

        
        subplot(length(InterestingKinases),1,ifr)

        yyaxis left
        %plot(myVec.ref,'o-','Linewidth',1,'markersize',10); hold on;
        stem(myVec.ref,'o-','Linewidth',1,'markersize',10); hold on;
        if FIG3_SHOW_PCA
            plot(myVec.SVD,'+-','Linewidth',1,'markersize',8);
        end
        %plot(myVec.NN,'x-','Linewidth',0.5,'markersize',8)
        %plot(myVec.NNrescale,'x-','Linewidth',0.5,'markersize',8)
        yyaxis right
        stem(myVec.NNrescale,'x-','Linewidth',0.5,'markersize',8)


        if FIG3_SHOW_PCA
            pcalabel = sprintf('PCA (%d vectors)',nmaxx2);
            legend('full',pcalabel,'predicted')
        else
            legend('true','predicted')
        end
        titlestring = [sKinaseNo ' ' sKinaseName ' ' sTrainStatus];
        title(titlestring,'interpreter','none');

        if FIG3_DETAIL
            xlim(rr);%
        else
            xlim([0 5300]);
        end

        if USELOGSCALE_FIG3
            set(gca,'yscale','log');
            ylim([1e-3 10]);
        else
            yyaxis left
            ylim([-0.5 1.5]);
            yyaxis right
            ylim([-0.5 1.5]*normfacc);

        end

    end

    if SAVEFIGURES
        FigureName = sprintf('%sFig%d_AffinitySpectra_IndivKinases_%s',...
            SaveDir,FigNo,netfileID);

        if FIG3_DETAIL
            FigureName = [FigureName, '_detail'];
        end
        %FigureName=[SaveDir, 'TPFPCutbyKinase', netfilename ];
        FigW = 11;
        FigH = 11;
        SaveThisFigureLarge     
    end
end
%% figure 900 - scatter plots of actual / SVD / NN predictions
% plots to browse
PLOT_THIS = true;%false;%
if PLOT_THIS

    FigNo = 900;
    maxkin_900 = 3;

    figure(FigNo); clf;


    trainstatus={'train','validation','test'};

    nkin_900 = min(length(InterestingKinases),maxkin_900);

    for ifr=1:nkin_900 % loop over selected kinases

        
        
        % features of this kinase
        myK = InterestingKinases(ifr); % kinase index in the set of 455 ie.LigandInfoKinases
        %myH = find(LCAMTred(myK,:)>=4); % hits for this kinase
        %myNH = find(LCAMTred(myK,:)<4); % nonhits for this kinase
        

        myTVC = LIK.TVC(myK); % train status -- 1,2,3 / train,val,test

        if(exist('myVec','var')), clear 'myVec'; end

        if REDO_PCA
            %Method1 = 4;
            %Method2 = 3;
            Fig900_m = [2,1,3,4];
            %myVec.ref = Lthis.check(myK,Lthis.ActiveLigands);
            myVec.ref = LCAMTs(myK,Lthis.ActiveLigands);
            %myVec.SVD = RSA{Method1}.check(myK,:);
            %myVec.NN = RSA{Method2}.check(myK,:);
            %myVec.Lin = RSA{3}.check(myK,:);
            %myVec.Combo = RSA{4}.check(myK,:);
            myHs = find(LIK.indivhits(myK,Lthis.ActiveLigands));
            myNHs = find(LIK.indivnonhits(myK,Lthis.ActiveLigands));
        else
            myH = find(LIK.indivhits(myK,:)); % hits for this kinase
            myNH = find(LIK.indivnonhits(myK,:)); % nonhits for this kinase

            myHs = intersect(myH,Lthis.ActiveLigands);
            myNHs = intersect(myNH,Lthis.ActiveLigands);
            myVec.ref = LCAMTs(myK,:)';
            myVec.SVD = Lthis.checkSVD(myK,:)';
            myVec.NN = Lthis.checkNN(myK,:)';
        end


        %WW = {myVec.ref,myVec.SVD,myVec.NN};
        myVec = struct2table(myVec);

        % strings for labelling
        titlestring = kinlabel(myK);
        %sKinaseNo = sprintf('Kin.%d [%d]', LigandInfoKinases(myK),myK);%
        %sKinaseName = KGOT.Name{LigandInfoKinases(myK)};
        %sTrainStatus = trainstatus{myTVC};

        %
        if REDO_PCA
            %sSetLabel = {'actual',RSA{Fig900_m(1)}.label,RSA{Fig900_m(2)}.label};
            ShortMethodLabel = {'NN','SVD','Lin','Combo'};
        else
            sSetLabel = {'actual','SVD','NN+SVD'};
        end
        
        Nsp = 2;
        for sp=1:Nsp % subplot index (two rows)

            subplot(Nsp,nkin_900,ifr + (sp-1)*nkin_900)

            if REDO_PCA
                WW1 = RSA{Fig900_m(1)}.check(myK,:); % common x axis
                if sp==1 % first row from the top
                    WW2 = myVec.ref;
                    ystring = 'Actual';
                else
                    WW2 = RSA{Fig900_m(sp)}.check(myK,:);
                    %ystring = RSA{Fig900_m(sp)}.label;
                    ystring = ShortMethodLabel{Fig900_m(sp)};
                end
                %xstring = RSA{Fig900_m(1)}.label;
                xstring = ShortMethodLabel{Fig900_m(1)};
            else
                if sp==1 % sets are {ref,SVD,NN}
                    sset = [2,3];
                elseif sp==2
                    sset = [2,1];
                end

                WW1 = table2array(myVec(:,sset(1)));
                WW2 = table2array(myVec(:,sset(2)));
                xstring = sSetLabel{sset(1)};
                ystring = sSetLabel{sset(2)};
            end

            scatter(WW1(myNHs),WW2(myNHs),5,'k',"filled"); hold on
            scatter(WW1(myHs),WW2(myHs),5,'r',"filled"); hold on


            plot([0,1],[0,1],'k:');

            subtitle(titlestring,'interpreter','none');
            %title(sKinaseNo,'interpreter','none')

            xlabel(xstring);
            ylabel(ystring);


            xlim([-0.2 1.2]);
            ylim([-0.2 1.2]);
            legend('NH','hits')
            legend('location','southeast')
            box on

            FS=[10, 10, 10, 10];
            FormatThisFigureNew;

        end

    end % loop over "interesting" kinases

    if SAVEFIGURES
        FigureName = sprintf('%sFig%d_AffinityScatter_IndivKinases_%s',...
            SaveDir,FigNo,netfileID);
        %FigureName=[SaveDir, 'TPFPCutbyKinase', netfilename ];
        SaveThisFigureLarge       
    end
    
end
%% fig 8804 -- histograms of predicted activities by component
PLOT_THIS = true;%false;%
if PLOT_THIS
    figure(8804)
    clf;

    AffinityLabel = 'Scaled Affinity';
    
    ScaleExtension = 0.49;%0.1;%
    AFF_SCALE_8804 = [0, 1] + ScaleExtension*[-1,1];
    LOG_SCALE_8804 = true;false;%

    NormType = 'count';
    if ~LOG_SCALE_8804, NormType = 'probability'; end
    gridd = AFF_SCALE_8804(1):0.02:AFF_SCALE_8804(2);
    tts1 = 'Training';
    tts2 = 'Other';

    FS = [10,8,8,8]; % fonts used in FormatThisFigure

    if REDO_PCA

        METHOD_8804(1) = 1; % 2=SVD, 1=NN, 3=Lin+, 4=Combo
        METHOD_8804(2) = 3;
        for iset=1:2
            DATA_8804{iset} = zeros(size(Lthis.check));% NN=1
            DATA_8804{iset}(:,Lthis.ActiveLigands) = RSA{METHOD_8804(iset)}.check; % SVD=2
        end
    else
        DATA_8804{1} = Lthis.checkSVD;
        DATA_8804{2} = Lthis.checkNN;
    end

    %SetLabel = {'training','test/validation'};
    SetLabel = {'training','validation','test'};

    for igroup = 1:3
    subplot(3,2,2*igroup-1)

    histogram(DATA_8804{1}(:,:),gridd,'normalization',NormType,'displaystyle','stairs')
    hold on;
    %histogram(Lthis.check(TrainSet,:),gridd,'normalization',NormType)
    histogram(DATA_8804{1}(Linds3{igroup,1}==1),gridd,'normalization',NormType)
    histogram(DATA_8804{1}(Linds3{igroup,2}==1),gridd,'normalization',NormType)
    if REDO_PCA
        title([RSA{METHOD_8804(1)}.label, ' Predicted Group Affinities'])
    else
        title 'SVD Predicted Group Affinities'
    end
    subtitle([netstring ' ' SetLabel{igroup} ' set'],'interpreter','none')
    xlabel(AffinityLabel);
    ylabel(NormType)
    legend('all',[SetLabel{igroup} ' hits'],[SetLabel{igroup} ' non hits'])
    legend('interpreter','none')
    if LOG_SCALE_8804, set(gca,'yscale','log'); end
    FS = [10,8,8,8]; % fonts used in FormatThisFigure
    FormatThisFigure

    subplot(3,2,2*igroup)
    %NormType = 'probability';
    histogram(DATA_8804{2}(:,:),gridd,'normalization',NormType,'displaystyle','stairs')
    hold on;
    %histogram(DATA_8804{2}(TrainSet,:),gridd,'normalization',NormType)
    histogram(DATA_8804{2}(Linds3{igroup,1}==1),gridd,'normalization',NormType)
    histogram(DATA_8804{2}(Linds3{igroup,2}==1),gridd,'normalization',NormType)
    if REDO_PCA
        title([RSA{METHOD_8804(2)}.label, ' Predicted Group Affinities'])
    else
        title 'NN Predicted Group Affinities'
    end
    subtitle([netstring SetLabel{igroup} ' set'],'interpreter','none')
    %subtitle([netstring ' training set'],'interpreter','none')
    xlabel(AffinityLabel);
    ylabel(NormType)
    %legend('all',tts1)
    legend('all',[SetLabel{igroup} ' hits'],[SetLabel{igroup} ' non hits'])
    legend('interpreter','none')
    if LOG_SCALE_8804, set(gca,'yscale','log'); end
    FS = [10,8,8,8]; % fonts used in FormatThisFigure
    FormatThisFigure

   
    end

    % save
    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig8804_AffinityHistograms_', netfileID];
        if RAW_SVD
            FigureName = ['rawSVD_' FigureName];
        end

        FigureName = [FigureName, sprintf('_%d',myK)];

        SaveThisFigureLargePortrait

    end

    clear NormType; % figure specific switch[es]
    clear LOG_SCALE_8804 AFF_SCALE_8804
    %clear DATA_8804 METHOD_8804 

end % figure 8804

%% fig 8830 -- histograms of predicted affinities by kinase
PLOT_THIS = true;%false;%
if PLOT_THIS
    figure(8830)
    clf;

    AffinityLabel = 'Scaled Affinity';
    
    ScaleExtension = 0.49;%0.1;%
    AFF_SCALE_8830 = [0, 1] + ScaleExtension*[-1,1];
    LOG_SCALE_8830 = true;false;%

    MAX_KINASES_8830 = 3;

    NormType = 'count';
    if ~LOG_SCALE_8830, NormType = 'probability'; end
    gridd = AFF_SCALE_8830(1):0.02:AFF_SCALE_8830(2);
    tts1 = 'Training';
    tts2 = 'Other';

    FS = [10,8,8,8]; % fonts used in FormatThisFigure

    if REDO_PCA
        METHOD_8830(1) = 1; % 2=SVD, 1=NN, 3=Lin+, 4=Combo
        METHOD_8830(2) = 3;
        for iset=1:2
            DATA_8830{iset} = zeros(size(Lthis.check));% NN=1
            DATA_8830{iset}(:,Lthis.ActiveLigands) = RSA{METHOD_8830(iset)}.check; % SVD=2
        end
    else
        DATA_8830{1} = Lthis.checkSVD;
        DATA_8830{2} = Lthis.checkNN;
    end

    
    % loop over kinase list
    nkin_8830 = min(length(InterestingKinases),MAX_KINASES_8830);
    for ikin=1:nkin_8830
        for icol=1:2

            subplot(nkin_8830,2,2*ikin-2+icol)
            % histogram for a selected kinase

            myK = InterestingKinases(ikin); % kinase index in the set of 455 ie.LigandInfoKinases
            % filtered hits / nonhits of this kinase, absolute ligand index
            myH =  Lthis.ActiveLigands(find(LIK.indivhits(myK,Lthis.ActiveLigands))); % hits for this kinase, filtered
            myNH = Lthis.ActiveLigands(find(LIK.indivnonhits(myK,Lthis.ActiveLigands))); % nonhits for this kinase, filtered




            H = cell(1,3);
            H{1}=histogram(DATA_8830{icol}(:,:),gridd,'normalization',NormType,'displaystyle','stairs'); hold on;
            H{2}=histogram(DATA_8830{icol}(myK,myH),gridd,'normalization',NormType);
            H{3}=histogram(DATA_8830{icol}(myK,myNH),gridd,'normalization',NormType);
            H{2}.FaceColor = 'r';
            H{2}.FaceAlpha = 0.5;
            H{3}.FaceAlpha = 0.5;
            if REDO_PCA
                title([RSA{METHOD_8830(icol)}.label, ' Predicted Group Affinities'])
            else
                title 'SVD Predicted Group Affinities'
            end
            %subtitle(ttsone,'interpreter','none')
            subtitle(kinlabel(myK),'interpreter','none');
            xlabel(AffinityLabel);
            ylabel(NormType)

            % nstring = kinlabel(myK);
            % xf = 0.35; yf = 0.95;
            legend('all kinases','hits','nonhits');
            legend('interpreter','none')
            if LOG_SCALE_8830, set(gca,'yscale','log'); end
            FS = [10,8,8,8]; % fonts used in FormatThisFigure
            FormatThisFigure



        end
    end

    if SAVEFIGURES

        FigureName=[SaveDir, 'Fig8830_AffinityHistograms_', netfileID];
        if RAW_SVD
            FigureName = ['rawSVD_' FigureName];
        end

        %FigureName = [FigureName, sprintf('_%d',myK)];

        SaveThisFigureLargePortrait

    end

    clear NormType; % figure specific switch[es]
    clear LOG_SCALE_8830 AFF_SCALE_8830
    clear DATA_8830 METHOD_8830

end

%% stop here

return

%% sanity checks ...
figure(12233);
clf
% look at LCAMTred size NT x NLC
plot(1:NLC,min(LCAMTred,[],1));
hold on;
plot(1:NLC,max(LCAMTred,[],1));

%% more old plots
figure(12234);
clf
% look at LCAMTred size NT x NLC
plot(1:NLC,min(LCAMTs,[],1));
hold on;
plot(1:NLC,max(LCAMTs,[],1));
plot(1:NLC,mean(LCAMTs,1));
ylim([-0.5 1.5]);

figure(12235);
clf
subplot(4,1,1);
pcolor(Linds{1,1});
shading flat
colormap jet
title 'Train x Hits'

subplot(4,1,2);
pcolor(Linds{1,2});
shading flat
title 'Train x Non-hits'

subplot(4,1,3);
pcolor(Linds{2,1});
shading flat
title 'Other x Hits'

subplot(4,1,4);
pcolor(Linds{2,2});
shading flat
title 'Other x Non-hits'

%% utility functions (to be used only locally)

function whc = whistcounts(values,weights,binedges)
% works like the regular histcount(values,binedges) 
% but each entry in vlaues is counted with the corresponding weight

% caveats - assume all sizes match, i.e.:
% values(N,1) ; weights(N,1) ; binedges(N+1,1)

% binedges assumed strictly increasing
% weights should be positive

N = length(values); % of items in the sample
M = length(binedges) - 1; % of bins

if ...
        min(size(values)) >1 || min(size(weights)) >1 || min(size(binedges)) >1 || ...
        length(weights) ~= N || ... % weights should match 
        M < 1 || ... % at least one bin
        min(binedges(2:end) - binedges(1:end-1)) <= 0 % bin edges strictly increasing
    fprintf('error whistcounts\n')
    return
end

whc = zeros(1,M);

% set up a bin label vector
binvec = zeros(1,N);

% loop over the bins and label each entry that is above the lower edge
for ibin = 1:M+1
    binvec(values>=binedges(ibin)) = ibin;
end
% the top bin also includes values matching the upper edge
binvec(values==binedges(M+1)) = M;

% bin=0 or bin = M+1 means outside the range


% convert to counts
for ibin = 1:M
    %whc(ibin) = length(find(bin==ibin));% plain hist counts
    whc(ibin) = sum(weights(binvec==ibin));% weighted
end



return

end