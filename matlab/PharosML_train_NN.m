% replication code #3 in the sequence

% NN setup


if ~exist('ThisNetConfig','var') % ThisNetConfig used to set the sizes from outside
    ThisNetConfig.net = [1000,500,200];%[2000,1000,500];%[150,150];%[
    ThisNetConfig.n1 = 50;%100;%200;%
    ThisNetConfig.n2 = 30;
end

% initialize the net object
net = feedforwardnet(ThisNetConfig.net,'trainscg'); %
% number of principal vectors used for input / output

nmax1 = ThisNetConfig.n1;%
nmax2 = ThisNetConfig.n2;%

% refine training parameter settings 

% stopping criteria
net.trainParam.epochs = 10000; % let it run longer (default is 1000 epochs)
net.trainParam.time = 3600; % time limit ("wall clock time")
net.trainParam.goal = 1e-5;
net.trainParam.max_fail=50;%100;%30;%5;%

% use the same training / val / test set each time
net.divideFcn = 'divideind';
if USE_TRACE_SPLIT % if we use the same training / val / test set for multiple runs
    net.divideParam.trainInd = tr.trainInd;
    net.divideParam.valInd = tr.valInd;
    net.divideParam.testInd= tr.testInd;

    %extra safety -- make sure LIK.TVC matches the trace split if we
    %choose to keep the train/validation/test assignments from the
    %previous train run
    LIK.TVC(tr.trainInd) = 1;
    LIK.TVC(tr.valInd  ) = 2;
    LIK.TVC(tr.testInd ) = 3;

else % preferred mode - use the TVT split in the LIK object
    net.divideParam.trainInd = find(LIK.TVC==1);
    net.divideParam.valInd = find(LIK.TVC==2);
    net.divideParam.testInd= find(LIK.TVC==3);
end


%% choose the input and output data

% Lthis is the generic object used with NN training; 
% for replication purposes Lthis is the same as LpcaS
% this is / was used for alternative PCA / projection schemes
Lthis = LpcaS; 

% for replication purposes .csscore is the same as .sscore 
Lthis.csscore = Lthis.sscore;

% input data
Xdata = Tpca.score(LigandInfoKinases,1:nmax1);

%% train the NN

Ydata = Lthis.csscore(:,1:nmax2); % Y is NLK x p

% TEST set the *train* set Y to 0 to make sure it is not used in the NN
% training
% Ydata(net.divideParam.testInd,:)=0;
% TEST

[net, tr] = train(net,Xdata',Ydata');

%% save into a file with a descriptive name

if SAVE_THIS_RUN

    mfn = sprintf('%d',nmax1);
    for j=1:length(net.layers)-1
        mf2 = sprintf('_%d',net.layers{j}.dimensions);
        mfn = [mfn,mf2];
    end
    mfn2 = sprintf('_%d',nmax2);
    mfn =[mfn mfn2];

    mfn = sprintf('single_%s_%d_%.0fs_%depochs_perf%.3e.mat',mfn,randi(100000),...
        tr.time(end),tr.epoch(end),tr.best_perf);

    mfn = ['redoPCA_',mfn];

    CleanSave;

end
%end
%% get the predicted Y


YY = net(Xdata');
YY = YY';



%% map back to the original variables

Lthis.sscoreNN = YY;

nmaxx2 = nmax2;

% undo the [0,1] scaling on the PCA side
Lthis.scoreNN = Lthis.minscore(:,1:nmaxx2) + ...
    Lthis.sscoreNN .* (Lthis.maxscore(:,1:nmaxx2) -...
    Lthis.minscore(:,1:nmaxx2));

% map back to the original representation
Lthis.checkNN = Lthis.mean + Lthis.scoreNN * Lthis.coeff(:,1:nmaxx2)';
% NOTE: this is comparable to the SCALED LG affinity vectors

%% plots for comparison


%% reconstruct the pre-PCA vectors using only nmax2 principals
% Lthis.check = Lthis.score(:,1:nmax2) * Lthis.coeff(:,1:nmax2)' + Lthis.mean;
Lthis.check = Lthis.score(:,1:nmaxx2) * Lthis.coeff(:,1:nmaxx2)' + Lthis.mean;

%% aside: count the number of ligand cluster interactions for each target
if ~exist('LCAMTb','var')
    LCAMTb = LCAMTs;
    LCAMTb(LCAMTs>0)=1;
end

% lrr is the set of ligands to be shown

%if REDO_PCA
lrr = LpcaS.TrainLig;
%else
%    lrr = 1:NLC;
%end

% hit counts by ligand
TLC = sum(LCAMTb(:,lrr),2);

return;
if ~MAKE_NNTRAIN_PLOTS
    return
end
%% only plots from here on

% Utility function
kinlabel = @(i) sprintf('Kin.%d [%d] %s %s',...
    LigandInfoKinases(i),...
    i,...
    KGOT.Name{LigandInfoKinases(i)},...
    testlabel{LIK.TVC(i)});

figure(3)
clf;

if exist('interesting','var')
    ione = interesting(1);
    itwo = interesting(2);
else
    ione = 446;%155;%160;%155;%85;%108;% 446;% 349;
    itwo = 349;% 350;
end
plot(LCAMTs(ione,lrr),'o-','Linewidth',1,'markersize',10);
hold on;
plot(Lthis.check(ione,lrr),'+-','Linewidth',1,'markersize',8);
plot(Lthis.checkNN(ione,lrr),'x-','Linewidth',0.5,'markersize',8)
pcalabel = sprintf('PCA (%d vectors)',nmax2);
legend('full',pcalabel,'predicted')
tstring = kinlabel(ione);
title(tstring,'interpreter','none');
xlim([0 length(lrr)]);%xlim([500 1000]);%

%%
figure(803)
clf

subplot(2,2,2)
gridd=0:0.01:1;
NormType = 'probability';
histogram(Lthis.checkNN(:,lrr),gridd,'normalization',NormType)
hold on;
histogram(Lthis.checkNN(ione,lrr),gridd,'normalization',NormType)
title 'NN Predicted Group Affinities'
xlabel 'Affinity Measure (log_{10})'
ylabel(NormType)
legend('all',tstring)
legend('interpreter','none')
set(gca,'yscale','log')
subplot(2,2,1)
histogram(Lthis.check(:,lrr),gridd,'normalization',NormType)
hold on;
histogram(Lthis.check(ione,lrr),gridd,'normalization',NormType)
set(gca,'yscale','log')
legend('all',tstring)
legend('interpreter','none')
title 'SVD Predicted Group Affinities'
xlabel 'Affinity Measure (log_{10})'
ylabel(NormType)

subplot(2,2,3)
histogram(TLC,0:2:300)
xlabel 'Number of high affinities'
title 'Kinases and number of LG hits'
ylabel 'Kinase count'
set(gca,'yscale','log')
hold on
histogram(TLC(strcmp(KGOT.TDL(LigandInfoKinases),'Tchem')),0:2:300);
%histogram(TLC(strcmp(KGOT.TDL(LigandInfoKinases),'Tclin')),0:2:300);
legend('both','Tchem')

subplot(2,2,4)
histogram(TLC,0:2:300)
xlabel 'Number of high affinities'
title 'Kinases and number of LG hits'
ylabel 'Kinase count'
set(gca,'yscale','log')
hold on
%histogram(TLC(strcmp(KGOT.TDL(LigandInfoKinases),'Tchem')),0:2:300);
histogram(TLC(strcmp(KGOT.TDL(LigandInfoKinases),'Tclin')),0:2:300);
legend('both','Tclin')

%% Fig 904: scatter plots of kinase affinities (NN vs PCA) by hit count
figure(904)
clf
comprange = [0,10,20,50,300];

for ip=1:4
    subplot(2,2,ip)
    for ii=1:size(Lthis.check,1)
        if TLC(ii) >= comprange(ip) && TLC(ii) < comprange(ip+1)
            scatter(Lthis.check(ii,lrr),Lthis.checkNN(ii,lrr),10,TLC(ii)*ones(1,length(lrr)),"filled")
            hold on
            plot([0,1],[0,1],'k:');
        end
    end
    tstring = sprintf('Between %d and %d hits',comprange([ip,ip+1]));
    title(tstring)
    xlabel 'PCA predicted'
    ylabel 'NN + PCA predicted'
    title(tstring)
    ylim([-0.2 1.2]);
    colormap jet
    colorbar
    box on
end

%%
figure(906)
clf

maxc = 5;
for ipl=1:maxc

    ione = interesting(ipl);

    subplot(2,maxc,ipl)
    thehits = lrr(LIK.indivhits(ione,lrr)>0);
    thenonhits = lrr(LIK.indivnonhits(ione,lrr)>0);
    therest = setdiff(1:NLC,lrr);

    scatter(Lthis.check(ione,thenonhits),Lthis.checkNN(ione,thenonhits),10,'k',"filled")
    hold on
    scatter(Lthis.check(ione,thehits),Lthis.checkNN(ione,thehits),10,'r',"filled")
    if ~isempty(therest)
        scatter(Lthis.check(ione,therest),Lthis.checkNN(ione,therest),10,'g',"filled")
    end
    plot([0,1],[0,1],'k:');

    tstring = kinlabel(ione);
    title(tstring,'interpreter','none')
    xlabel 'PCA predicted'
    ylabel 'NN + PCA predicted'
    ylim([-0.2 1.2]);
    %legend('one','two')
    box on

    subplot(2,maxc,ipl+maxc)



    scatter(Lthis.check(ione,thenonhits),LCAMTs(ione,thenonhits),10,'k',"filled")
    hold on
    scatter(Lthis.check(ione,thehits),LCAMTs(ione,thehits),10,'r',"filled")
    plot([0,1],[0,1],'k:');


    %tstring = sprintf('Between %d and %d hits',comprange([ip,ip+1]));
    % tstring = sprintf('Kinase no. %d [%d] and %d [%d]',...
    %     LigandInfoKinases(ione),ione,...
    %     LigandInfoKinases(itoo),itoo);
    tstring = sprintf('Kinase no. %d [%d]',...
        LigandInfoKinases(ione),ione);%
    title(tstring)
    xlabel 'PCA predicted'
    ylabel 'actual'
    ylim([-0.2 1.2]);
    %legend('one','two')
    box on

end

return

