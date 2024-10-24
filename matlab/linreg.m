% linear regression between the pca
% predict Lthis.sscore from Xdata

%trlist = find(LIK.TVC==1); % list of training kinases - use TrainSet
%instead

Lthis.sscoreLin = zeros(size(Lthis.sscore)); % will hold the linear predictions

%NTdim = rank(Tpca.score);
NTdim = rank(Tpca.score(LigandInfoKinases,:));
% NTdim = rank(Tpca.score(LigandInfoKinases(TrainSet),:));
nmax1_lin = NTdim; % max possible
%nmax1_lin = nmax1;% dimension of input set from outside (used in the NN)
nmax2_lin = Lthis.NSV;%max possible?
%nmax2_lin = nmaxx2;

%fprintf('linreg: %d %d \n',nmax1_lin,nmax2_lin);

CClin = zeros(nmax1_lin,nmax2_lin); % Xdata is size (NLK,nmax1) but we predict NSV dimensional coordinate vectors
BClin = zeros(NLK,nmax2_lin); % offset

if nmax1_lin ~= nmax1
    XdataLin = Tpca.score(LigandInfoKinases,1:nmax1_lin);
else
    XdataLin = Xdata;
end

% linear regression to predict coordinate vectors
NormBeta = zeros(nmax2_lin,1);
for myind=1:nmax2_lin% needs to be done one dimension at a time
    % lin regression to predict the current component (myind)
    % from ALL components of the input vector
    % using all vectors pairs available in the training set
    [Mdl,FitInfo] = fitrlinear(XdataLin(TrainSet,1:nmax1_lin),Lthis.sscore(TrainSet,myind),...
        'learner','leastsquares',...
        'regularization','lasso');%'svm');%
    BClin(:,myind) = Mdl.Bias;
    CClin(:,myind) = Mdl.Beta;
    % checks
    NormBeta(myind) = norm(Mdl.Beta);
    if(norm(Mdl.Beta)>0.05)
         zz = Mdl.Bias + XdataLin(TrainSet,:)*Mdl.Beta;
         ZisInd = myind;
         ZisMdl = Mdl;
    end
end

%Lthis.sscoreLin = zeros(NLK,Lthis.NSV); % keep it full size but only fill the number of dimensions used
Lthis.sscoreLin(:,1:nmax2_lin) = BClin + XdataLin*CClin;
return
%%
figure(1000);
clf
MyComp = 3;
plot(Lthis.sscore(:,MyComp),Lthis.sscoreLin(:,MyComp),'b.');
hold on
plot(Lthis.sscore(TrainSet,MyComp),Lthis.sscoreLin(TrainSet,MyComp),'ro');
plot([0,1],[0,1],'k--')
%% compute R^2 for all components
Rsq = zeros(NLK,1);
for MyComp=1:nmax2_lin
    MeanScores = mean(Lthis.sscore(:,MyComp));
    VarScores = Lthis.sscore(:,MyComp) - repmat(MeanScores,NLK,1);
    RezScores = Lthis.sscore(:,MyComp) - Lthis.sscoreLin(:,MyComp);
    zzrez = sum(RezScores.^2,1);
    zzvar = sum(VarScores.^2,1);
    Rsq(MyComp) = 1 - zzrez / zzvar;
end