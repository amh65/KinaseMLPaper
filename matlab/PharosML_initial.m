% replication code #1 in the sequence
% reads the relevant data files and generates the structures used for
% model training and analysis
% NOTE -- the PCA for ligand affinities is based only on the TRAIN set
%         this is done in a separate script
%         the Lpca structure here is optional 

%% processing options (truly OPTIONAL)

% construct separate activity matrices by activity type
MAKE_AFFINITY_MATRIX_BY_ACTIVITY_TYPE = false;%true;

% perform PCA for the full ligand affinity matrix
% [PCA is typically done using only the train set, this is for comparison]
DO_FULL_LPCA = false;

%% Settings
%  keep the input filenames and processing options here
clear PML; % this will hold the inut file names and folder[s]
PML.DataPath = '../input_files/SourceData/';

% identifying info on kinase targets
PML.KinaseTargetFile = [PML.DataPath,'TCRDv613_KinaseDomainAnnotations.xlsx'];
% NOTE: column 1098 (the last one) in this file is not to be used
%       it contains a string of all binary entries '1000100.. '


% ligand - kinase interactions from a Pharos query
PML.LigandQueryFile = [PML.DataPath,'LigandQuery.xlsx'];

% ligand similarity data from Steve
PML.LigandSMILESFile = [PML.DataPath,'LigandSMILES.csv'];
PML.LigandClusterFile = [PML.DataPath,'LigandClusters.xlsx'];



%% load data on targets / kinases (name, ID, structure)
% kinase: basic and strucutral annotation data

fprintf('Reading Kinase (target) data from: %s  ..', PML.KinaseTargetFile);

KT = readtable(PML.KinaseTargetFile);
% add meaningful names where needed
KT.Properties.VariableNames{1} = 'TCRD_ID';
KT.Properties.VariableNames{2} = 'Name';
KT.Properties.VariableNames{3} = 'Description';
KT.Properties.VariableNames{4} = 'UniProt';
KT.Properties.VariableNames{5} = 'HGNC_Symbol';
KT.Properties.VariableNames{6} = 'TDL';


[TIndList, ia, ~ ]= unique(KT.TCRD_ID); % list of target IDs from the kinase list
% [C ia ic] = unique(A) --> A=C(ic), C=A(ia)

KTsort = KT(ia,:); % sorted by TCRD_ID, consistent with the data vectors TargetVec etc

NT = length(TIndList);



TargetVec = table2array(KTsort(:,7:1097)); % list of binary features (structural keys style)
% NOTE: there is an extra [spurious] column in the excel version of the
% data file: 'TCRDv613_KinaseDomainAnnotations.xlsx'; 
% we do not include it in the binary vector TargetVec 

fprintf('. done\n');

% Kinases: gene ontology (GO) annotations

PML.GOAFile = [PML.DataPath,'TCRDv613_GOAnnotations.xlsx'];

KGOT = readtable(PML.GOAFile);
% add meaningful names where needed
KGOT.Properties.VariableNames{1} = 'TCRD_ID';
KGOT.Properties.VariableNames{2} = 'Name';
KGOT.Properties.VariableNames{3} = 'Description';
KGOT.Properties.VariableNames{4} = 'UniProt';
KGOT.Properties.VariableNames{5} = 'HGNC_Symbol';
KGOT.Properties.VariableNames{6} = 'TDL';


[GOTIndList, ia, ic ]= unique(KGOT.TCRD_ID); % list of target IDs from the kinase list
% [C ia ic] = unique(A) --> A=C(ic), C=A(ia)

KGOTsort = KGOT(ia,:);

NGT = length(GOTIndList);


GOTargetVec = table2array(KGOTsort(:,7:end)); % list of binary features (structural keys style)

% Joint kinase feature vectors

JVec = [TargetVec, GOTargetVec];

% kinase identifying info, should be the same for KTsort and KGOTsort
KTBasicData = KTsort(:,1:6);

% important data
% PML -- input files used
% KTBasicData -- kinase identifying data,
%                matching KTsort and KGOT / KGOTsort
%                sorted by the TCRD_ID of kinases,
% JVec -- joint kinase feature vectors, from [TargetVec, GOTargetVec]
% NT -- number of kinases (635)

%return
%% read data on kinase target x ligand interactions
% TODO: clean ligand + properties list
%       ligand activitiy matrices (IC50 and the others)
%       target-target distance using common ligands

fprintf('Reading Ligand x Kinase interaction data from: %s  ..', PML.LigandQueryFile);

TKL = readtable(PML.LigandQueryFile); % this is a list of (target, ligand) pairs

% kinase-target list from the ligand query file
[IDList,ia,ic] = unique(TKL.id); % list of target IDs from the ligand list
% IDList has unique IDs;
% ia are indices in the original column / array, ie. IDList = TKL.id(ia)
% ic are indices in the unique list ie.  TKL.id = IDList(ic) )

if max(abs(IDList-TIndList)) > 0
    fprintf('**Fatal Error**\n\nTarget lists differ between target file and ligand file!\n');
    return
else
    % if IDList is consistent with the original kinase list, use the
    % sorting to label rows in the interaction file
    TKL.idindex = ic;
end

clear IDList ia ic

fprintf('. done\n');

%return

%% ligand lists

fprintf('Creating ligand list  ..');

% clean the entries of TKL, remove those with
% badly identified ligands and missing activity
TKLClean = TKL;

% optional: require PubChemID and / or ChEMBLID
%TKLClean = TKLClean(~isnan(TKLClean.LigandPubChemID),:);
%TKLClean = TKLClean(~strcmp(TKLClean.LigandChEMBLID,""),:);

% require SMILES (necessary for the ligand-ligand similarity analysis)
TKLClean = TKLClean(~strcmp(TKLClean.LigandSMILES,""),:);

% skip entries that don't list an TKL.activity value and type
TKLClean = TKLClean(~isnan(TKLClean.LigandActivity),:);
TKLClean = TKLClean(~strcmp(TKLClean.LigandActivityType,""),:);


% slow but safe -- make a list of SMILES and then find all corresponding
% identifying info: PubChemID, ChemBLID, name

[C,ials,icls] = unique(TKLClean.LigandSMILES,'first');
% C has the unique SMILES strings;
% ia are indices in TL; i.e. C = TL.LigandSMILES(ia)
% ic indices in C that point back to TL i.e. TL.LigandSMILES = C(ic)

TKLClean.LigandIndex = icls;

Ligand = TKLClean(ials,[4,5,6,7]);
% so for any [allowed] index ii ..
% TLC.LigandSMILES(ia(ii)) = Ligand.LigandSMILES(ii)
% LigandSMILES(ic(ii)) = TLC.LigandSMILES(ii)
NL = size(Ligand,1);

fprintf('. done\n');

%return
%% reference SMILES list
fprintf('Reading Ligand SMILES list from: %s  ..', PML.LigandSMILESFile);

LSMILES = readtable(PML.LigandSMILESFile);
% UniqueSMILES = unique(LSMILES);

% LSMILES is slightly larger than Ligand.LigandSMILES because there are a
% couple of ligands that are listed in the ligand x kinase interaction file
% but the corresponding activity is blank

% label the entries of Ligand with their index in LSMILES
for ilig = 1:NL
    is = find(strcmp(LSMILES.SMILES,Ligand.LigandSMILES{ilig}),1);
    Ligand.SMILESind(ilig) = is;
end

fprintf('. done\n');

%% ready made ligand clusters from a file

fprintf('Reading Ligand Clusters (based on SMILES) from: %s  ..', PML.LigandClusterFile);

LC = readtable(PML.LigandClusterFile);

NLC = length(unique(LC.Cluster));


% add a cluster index to the Ligand struct
% this maps each ligand to the ligand cluster corresponding to its SMILES
% index (i.e. the ligand group corresponding to the SMILES).
for ilig = 1:NL
    %is = find(strcmp(LSMILES.SMILES,Ligand.LigandSMILES{ilig}),1);
    Ligand.Cluster(ilig) = LC.Cluster(Ligand.SMILESind(ilig));
end

% also add a ligand cluster index to the interaction table
TKLClean.LCIndex = Ligand.Cluster(TKLClean.LigandIndex);

fprintf('. done\n');

%% Ligand x Kinase activity matrix

if MAKE_AFFINITY_MATRIX_BY_ACTIVITY_TYPE

    fprintf('Building Ligand x Kinase activity matrices ..')

    % group by activity type
    activities = unique(TKLClean.LigandActivityType);
    NA = size(activities,1);
    AM = cell(NA,1);

    % Ligand x Kinase interaction matrices
    LigandInteractionCount = zeros(NL,NA);


    for iact=1:NA % loop over activity types

        % filter for the current activity type
        v = find(strcmp(TKLClean.LigandActivityType,activities{iact}));

        % collect target x ligand activity data

        % AM will hold the recorded activitiy values of this type (ia)
        AM{iact} = zeros(NT,NL);

        % will count the distinct targets this ligand iteracts with
        LigandTargetCount = zeros(NL,1);

        for i=1:NL

            % for each ligand, identify target interactions
            u = find(TKLClean.LigandIndex==i);
            % keep only activity values for this activity type
            u = intersect(u,v);


            if ~isempty(u) % if we found any activity entries

                % look up the index (position) of the corresponding targets
                w = TKLClean.idindex(u);

                % filter out the unique targets for the given ligand
                % note - the interactions are sorted by target ID
                wu = unique(w);

                % this is the number of unique interactions for each ligand
                LigandTargetCount(i)=length(wu);

                % corresponding max activities
                actmax = zeros(length(wu),1);
                for j=1:LigandTargetCount(i)
                    q = find(w==wu(j)); % subset of indices in w that have the same w value
                    actmax(j) = max(TKLClean.LigandActivity(u(q)));
                end


                AM{iact}(wu,i)=actmax;%TL.LigandActivity(u);

            end
        end


        LigandInteractionCount(:,iact) = LigandTargetCount;

    end

    fprintf('. done\n');

end

%% ligand x kinase matrix with no distinction by activity type

% this will hold the highest observed activity value of **any** type
% the metrics are actually comparable and are defined the same way
AMT = zeros(NT,NL);

% will count the distinct targets this ligand iteracts with
LigandTargetCount = zeros(NL,1);

for i = 1:NL
    u = find(TKLClean.LigandIndex==i);
    if ~isempty(u)
        % look up the index (position) of the corresponding targets
        w = TKLClean.idindex(u);

        % filter out the unique targets for the given ligand
        wu = unique(w);

        % this is the number of unique interactions for each ligand
        LigandTargetCount(i)=length(wu);

        % corresponding max activities
        actmax = zeros(length(wu),1);
        for j=1:LigandTargetCount(i)
            q = find(w==wu(j)); % subset of indices in w that have the same w value
            actmax(j) = max(TKLClean.LigandActivity(u(q)));
        end


        AMT(wu,i)=actmax;%TL.LigandActivity(u);

    end


end

%% Ligand Group (cluster) x Kinase activity matrix BY ACTIVIY TYPE

if MAKE_AFFINITY_MATRIX_BY_ACTIVITY_TYPE

    fprintf('Building Ligand Cluster x Kinase activity matrices ..')

    % Ligand Group x Kinase interaction matrices
    LCInteractionCount = zeros(NLC,NA);

    for iact=1:NA % loop over activity types

        % filter for the current activity type
        v = find(strcmp(TKLClean.LigandActivityType,activities{iact}));

        % collect target x ligand activity data

        % LCAM will hold the recorded activitiy values of this type (ia)
        LCAM{iact} = zeros(NT,NLC);

        % will count the distinct targets this ligand iteracts with
        LCTargetCount = zeros(NLC,1);

        for i=1:NLC

            % for each ligand CLUSTER, identify target interactions
            u = find(TKLClean.LCIndex==i);
            % keep only activity values for this activity type
            u = intersect(u,v);


            if ~isempty(u) % if we found any activity entries

                % look up the index (position) of the corresponding targets
                w = TKLClean.idindex(u);

                % filter out the unique targets for the given ligand
                % note - the interactions are sorted by target ID
                wu = unique(w);

                % this is the number of unique interactions for each ligand
                % cluster
                LCTargetCount(i)=length(wu);

                % corresponding max activities
                actmax = zeros(length(wu),1);
                for j=1:LCTargetCount(i)
                    q = find(w==wu(j)); % subset of indices in w that have the same w value
                    actmax(j) = max(TKLClean.LigandActivity(u(q)));
                end


                LCAM{iact}(wu,i)=actmax;%TL.LigandActivity(u);

            end
        end


        LCInteractionCount(:,iact) = LCTargetCount;

    end

    fprintf('.done\n');

    % joint kinase x ligand cluster
    LCM = [LCAM{4}, LCAM{5}, LCAM{6}];

end

%% redo the by ligand cluster activity matrix with no distinction for activity type

fprintf('Building Joint (all activity type) Ligand Cluster x Kinase activity matrix ..')

% LCAMT will hold max recorded kinase x ligand group activity values 
% irrespective of activity type

LCAMT = zeros(NT,NLC);

% will count the distinct targets this ligand iteracts with
LCTargetCount = zeros(NLC,1);

for i=1:NLC

    % for each ligand cluster, identify target interactions
    u = find(TKLClean.LCIndex==i);

    if ~isempty(u) % if we found any activity entries

        % look up the index (position) of the corresponding targets
        w = TKLClean.idindex(u);

        % filter out the unique targets for the given ligand cluster
        wu = unique(w);

        % this is the number of unique interactions for each ligand
        % cluster
        LCTargetCount(i)=length(wu);

        % corresponding max activities
        actmax = zeros(length(wu),1);
        for j=1:LCTargetCount(i)
            q = find(w==wu(j)); % subset of indices in w that have the same w value
            actmax(j) = max(TKLClean.LigandActivity(u(q)));
        end


        LCAMT(wu,i)=actmax;%TL.LigandActivity(u);

    end
end


%LCInteractionCount(:,iact) = LCTargetCount;

fprintf('.done\n');

%% PCA for Kinase vectors

% silence warning about linearly independent components in PCA
% PCA will return a complete basis anyway (dimension of input - 1)
% Warning: 'Columns of X are linearly dependent to within machine precision.
%    Using only the first 354 components to compute TSQUARED.'
warning('off','stats:pca:ColRankDefX');

% built in pca
[Tpca.coeff,Tpca.score,Tpca.latent,Tpca.tsquared,Tpca.explained,Tpca.mu] = pca(JVec);
% see explanation in the Matlab help
% [coeff,score,latent,tsquared,explained,mu] = pca(X,Name,Value)

%% pre-PCA for Ligand affinities

% Kinases [identified by their index] that have ligand binding info
LigandInfoKinases = find(sum(LCAMT,2));

% LCAMT is the ligand cluster based equivalent of AMT
%   each row represents a kinase
%   columns are ligand clusters

% reduced affinity matrix - only rows with nonzero entries
LCAMTred = LCAMT(LigandInfoKinases,:);

% SCALED affinity matrix LCAMTs
% scale each component [column] to [0,1]
MMin = repmat(min(LCAMTred,[],1),size(LCAMTred,1),1);
MMax = repmat(max(LCAMTred,[],1),size(LCAMTred,1),1);
LCAMTs = (LCAMTred - MMin) ./ (MMax - MMin);
LCAMTs(MMax == MMin)=0;

if DO_FULL_LPCA

    % built-in PCA for thh ligand affinity matrix
    [Lpca.coeff,Lpca.score,Lpca.latent,Lpca.tsquared,Lpca.explained,Lpca.mu] = pca(LCAMTs);

    % check .. this is how the original vectors are obtained from the PCA
    Lpca.mean = repmat(mean(LCAMTs,1),size(LCAMTs,1),1);
    Lpca.check = Lpca.score * Lpca.coeff' + Lpca.mean;
    max(max(LCAMTs - Lpca.check))

    % scale AGAIN to [0,1] -- on the PCA side
    Lpca.minscore = repmat(min(Lpca.score,[],1),size(Lpca.score,1),1);
    Lpca.maxscore = repmat(max(Lpca.score,[],1),size(Lpca.score,1),1);
    Lpca.sscore = (Lpca.score - Lpca.minscore) ./ (Lpca.maxscore - Lpca.minscore);

    %NOTE: the NN is typically trained to predict Lpca.sscore = net(Tpca.score)

end

%% important variables (ligand sector)
% TKL -- ligand query file input (raw)
% TKLClean -- same as TKL, only entries that have nonempty
%             Ligand{ SMILES, Activity, ActivityType }
%             has .LigandIndex and .idindex [for kinases]
%             .idindex matches the sort order of all 635 kinases
%             .LigandIndex refers to ligands sorted by SMILES
% Ligand -- table with unique ligands in TKLClean,
%           identified and sorted by  SMILES
%           .Cluster (LC / LG) membership is appended by matching SMILES
%           the ligand SMILES list comes from one file and
%           the cluster assignment from another;
%           however, the SMILES list is really superfluous
%           .SMILESind is the position in LSMILES / LC