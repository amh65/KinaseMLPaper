% constructs SVD vectors etc. using only the training (+validation) sets

RedoPCA_CheckPrintouts = false;

%% built-in PCA using only the TRAIN kinases ("short basis")
%  alternatively we used both T+V kinases -- "long basis"
TrainSet = find(LIK.TVC==1);
OtherSet = find(LIK.TVC>=2);
NTK = length(TrainSet);
NOK = length(OtherSet);

% record these as fields in the LpcaS object
LpcaS.BasisSet = TrainSet; 
LpcaS.NonBasisSet = OtherSet;


% silence warning about linearly independent components in PCA
% PCA will return a complete basis anyway (dimension of input - 1)
% Warning: 'Columns of X are linearly dependent to within machine precision.
%    Using only the first 354 components to compute TSQUARED.'
warning('off','stats:pca:ColRankDefX');

%NOTE: LCAMTs contains only the [455] vectors for documented kinases (=with
%      some known ligand interaction  )
%      unlike LCAMT which has the entire [635] set and blank rows for the
%      non-documented kinases

LCAMTss = LCAMTs(TrainSet,:); % only the TRAIN kinases
[LpcaS.coeff,LpcaS.trainscore,LpcaS.latent,LpcaS.tsquared,LpcaS.explained,LpcaS.mu] = pca(LCAMTss);

NSV = NTK - 1; % number of singular vectors (# of vectors in the train set -1)

% reconstruction check .. this is how the original vectors are obtained from the PCA
LpcaS.trainmean = repmat(mean(LCAMTss,1),size(LCAMTss,1),1);
LpcaS.traincheck = LpcaS.trainscore * LpcaS.coeff' + LpcaS.trainmean;
if RedoPCA_CheckPrintouts
    fprintf('LpcaS reconstruction check (train only): %.4e\n',...
        max(max(abs(LCAMTss - LpcaS.traincheck))));
end

% also, Lpca.coeff(NLC,NTK) are orthonormal vectors...
% Z = LpcaS.coeff' * LpcaS.coeff;
% max(max(Z - eye(size(Z,1))))

% so the scores should be reobtained by scalar product
ScoresToo = (LCAMTss - LpcaS.trainmean) * LpcaS.coeff;
if RedoPCA_CheckPrintouts
    fprintf('LpcaS score check: %.4e\n',...
        max(max(abs(ScoresToo - LpcaS.trainscore))));
end


% get scores for the "other" kinases 
LpcaS.othermean = repmat(mean(LCAMTs(TrainSet,:),1),length(OtherSet),1);% uses the mean of the **train** vectors
LpcaS.otherscore = (LCAMTs(OtherSet,:) - LpcaS.othermean) * LpcaS.coeff;

% assemble the score vectors for both sets
LpcaS.score = zeros(NLK,NSV);
LpcaS.score(TrainSet,:) = LpcaS.trainscore;
LpcaS.score(OtherSet,:) = LpcaS.otherscore;
% also the two sets have different means...
LpcaS.mean = zeros(NLK,NLC);
LpcaS.mean(TrainSet,:) = LpcaS.trainmean;
LpcaS.mean(OtherSet,:) = LpcaS.othermean;

% scale AGAIN to [0,1] -- on the PCA side
LpcaS.minscore = repmat(min(LpcaS.score,[],1),size(LpcaS.score,1),1);
LpcaS.maxscore = repmat(max(LpcaS.score,[],1),size(LpcaS.score,1),1);
LpcaS.sscore = (LpcaS.score - LpcaS.minscore) ./ (LpcaS.maxscore - LpcaS.minscore);
% LpcaS.sscore = (LpcaS.score - min(LpcaS.score,[],1)) ./ ...
%     (max(LpcaS.score,[],1) - min(LpcaS.score,[],1));

% final reconstruction check 
LpcaS.check = LpcaS.score * LpcaS.coeff' + LpcaS.mean;
if RedoPCA_CheckPrintouts
    fprintf('LpcaS overall reconstruction check: %.4e\n',...
    max(max(abs(LCAMTs - LpcaS.check))));
end

% which ligands are [not] hit by the train set ?
LCAMTb = LCAMTs; % binary LCAMT, only 1's and 0's
LCAMTb(LCAMTs>0) = 1;

% identify ligands that interact with the train set
LpcaS.TrainLig = find(sum(LCAMTb(TrainSet,:),1)>0);
NLCS = length(LpcaS.TrainLig);
