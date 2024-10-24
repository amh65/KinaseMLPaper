% creates the LIK struct / table to identify hits, nonhits, etc.
LIK_VERBOSE = false;%true;%

if LIK_VERBOSE, fprintf('**LIKsetup ..'); end

%% list of kinases with some relevant group info
NLK = length(LigandInfoKinases);

LIK.level = zeros(NLK,1); % Tchem or Tclin
LIK.level(strcmp(KGOT.TDL(LigandInfoKinases),'Tchem'))=2;
LIK.level(strcmp(KGOT.TDL(LigandInfoKinases),'Tclin'))=1;

% add train / val / test info to LIK
if ~exist("USE_TRACE_SPLIT","var")
    % this is for backward compatibility, LIKsetup before 7/16/2024 relied
    % on an existing trace to assign the train / val /test sets
    % this is still the default if the switch is not found (and a tr object
    % exists)
    if exist('tr','var')
        USE_TRACE_SPLIT = true;
        fprintf('LIKsetup warning: USE_TRACE_SPLIT not found; assuming TRUE, will use existing tr object\n');
    else
        USE_TRACE_SPLIT=false;
        fprintf('LIKsetup warning: USE_TRACE_SPLIT not found; assuming FALSE b/c no tr was found\n');
    end
end

if USE_TRACE_SPLIT && ~exist('tr','var')
    USE_TRACE_SPLIT=false;
    fprintf('LIKsetup warning: USE_TRACE_SPLIT was TRUE but no trace found; changing to FALSE\n');
end

if USE_TRACE_SPLIT 
    LIK.TVC = zeros(NLK,1); % training (1), validation (2) or test (3)
    LIK.TVC(tr.trainInd)=1;
    LIK.TVC(tr.valInd)=2;
    LIK.TVC(tr.testInd)=3;
else 
    LIK.TVC = CurrentTVC;
end
%LIK.lind = (1:NLK)'; % index of original kinase (in KGOT)
%LIK.trueind=LigandInfoKinases'; never used this 


% part 1 ends ; part 2 follows
%% training and other vector indices
TrainSet = find(LIK.TVC==1);%find(LIK.TVC<2); TODO -- make sure this is consistent with previous runs 08-03-2024
OtherSet = find(LIK.TVC>=2);
ThreeSets = cell(3,1);% replaces TrainSet and OtherSet
for iset=1:3
    ThreeSets{iset} = find(LIK.TVC==iset); 
end

%% identify indices by training / other set and hit / non-hit
% binary flags on NLK x NLG arrays indicating hits and kinases sets 
% 1 if both criteria are true ,zero otherwise
% useful in building count-based stats

Linds = cell(2,2);
%Linds higlights four groups of indices:
%  Linds{:,1} are hits , {:,2} are non-hits
%  Linds{1,:} are training, {2,:} are "other" (eg validation and test)
for iset=1:2
    for ihit=1:2
        myset = zeros(size(LCAMTred));
        if iset==1 % training set
            myset(TrainSet,:)=1;
        elseif iset==2 % other set
            myset(OtherSet,:)=1;
        end

        % only the relevant entries remain 1
        if ihit==1 % hits remain 1
            myset(LCAMTred<4)=0;
        elseif ihit==2 % non-hits remain 1
            myset(LCAMTred>=4)=0;
        end
        Linds{iset,ihit}=myset;
    end
end
clear myset ihit iset

% flag hits and nonhits
LIK.indivhits = LCAMTred>=4;
LIK.indivnonhits = ~LIK.indivhits;

Linds3 = cell(3,2); % replaces Linds which only had two sets
% LindsThree{:,1} show hits , {:,2} show non-hits
% LindsThree{j,:} matches TVC (j=1 train, j=2, val, j=3 test)
blankmap = false(size(LCAMTred));
for iset=1:3
    for ihit=1:2
        Linds3{iset,ihit} = blankmap;
    end
    Linds3{iset,1}(ThreeSets{iset},:) = LIK.indivhits(ThreeSets{iset},:);
    Linds3{iset,2}(ThreeSets{iset},:) = LIK.indivnonhits(ThreeSets{iset},:);
end
%clear myhitmap blankmap
clear blankmap
clear ihit iset


% flag hits and nonhits
%LIK.indivhitstoo = Linds{1,1}+Linds{2,1};
%LIK.indivnonhitstoo = Linds{1,2}+Linds{2,2};

% count hits and non-hits for each kinase
% LIK.interactions = sum(Linds{1,1}+Linds{2,1},2);
% LIK.nonhits = sum(Linds{1,2}+Linds{2,2},2);
LIK.interactions = sum(LIK.indivhits,2);
LIK.nonhits = sum(LIK.indivnonhits,2);

% % sort by number of interactions (descending)
[~, LIK.sort] = sort(LIK.interactions,'descend');% LIK.sort(1) is the index with the highest # of interactions (350 in this case)
[~, LIK.unsort] = sort(LIK.sort);% LIK.unsort is the inverse mapping, gives the position in the sorted order for a given index, eg. LIK.unsort(350)=1

% NOTE: LASSO_MODE might not be present when LIKsetup is called

if exist('LASSO_MODE','var')
    if LASSO_MODE
        LIK.indivhits_red = LIK.indivhits(:,SLG);
        LIK.indivnonhits_red = LIK.indivnonhits(:,SLG);

        LIK.interactions_red = sum(LIK.indivhits_red,2);
        LIK.nonhits_red = sum(LIK.indivnonhits_red,2);

        [~, LIK.sort_red] = sort(LIK.interactions_red,'descend');% LIK.sort(1) is the index with the highest # of interactions (350 in this case)
        [~, LIK.unsort_red] = sort(LIK.sort_red);% LIK.unsort is the inverse mapping, gives the position in the sorted order for a given index, eg. LIK.unsort(350)=1
    end
end

if LIK_VERBOSE, fprintf('. done\n'); end

clear LIK_VERBOSE