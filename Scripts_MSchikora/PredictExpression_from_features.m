%% Description
%This script generates the figures for doing the prediction of native
%expression based on the different features.

%% Load data

clear all

cd('~/Google Drive/CareyLab/ExternalData/');
%WT RNAseq expression. Column 4:
FitFlow_Fast_Slow = dataset('file','./vanDijk15/FitFlow__Fast_Slow.tab','delimiter','\t'); 
FitFlow_Fast_Slow.Properties.VarNames = {'ORF','Gene','Location','WTexpression','x18_2308'};
%3'end contribution to expression:
End3_Contr = dataset('file','./Shalem13/Shalem13.tab','delimiter','\t'); 
End3_Contr.Properties.VarNames = {'ORF','End3_YFP','ClonedIntergenicSequence'};
%5'end contribution to expression (use glucose):
End5_Contr = dataset('file','./Keren13/Keren13_Expression.tab','delimiter','\t'); 
End5_Contr.Properties.VarNames(3) = {'End5_YFP_glucose'};

cd('~/Google Drive/CareyLab/random RNA library/Miki/Data/');
%All integrated gDNA random RNA library datasets, cointains many datasets:
load('Integrated_DataSets_2014_Adds.mat');
load('ORFtoSeq.mat');

cd('~/Google Drive/CareyLab/random RNA library/LBC/2016_September__PromoterInfluencesRNALife__CodonBias_Predicts_Keren_Residuals/')
%DataSet of Keren (5'end contribution) with some adds:
End5_Contr_adds = dataset('file','keren_rnaseq_cai.tab','Delimiter','\t');

cd('~/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Data/');
%UPF1 mutant data --> cointains the information for the 3'utr length:
load('./External Data/S_cer_UPF1_mutant_RNAseq.mat');
UPF1del_Dataset = G;

cd('~/Google Drive/CareyLab/random RNA library/Miki/Data/');
%curated DEE datasets:
load('All_RNA_seq_S_cer_integrated_All_DS.mat')

cd('~/Google Drive/CareyLab/ExternalData/DigitalExpressionExplorer/')
%SRA_S dataset:
SRAs_DS = dataset('file','./metadata.tab','Delimiter','\t');  

cd('~/Google Drive/CareyLab/random RNA library/Miki/Data/');
%DS of Julia:
load('JuliaDS_features.mat')

Julia_DS = DS;
%% Join datasets:

%A dataset that cointains the information for all the native genes that are in each of the 5', 3' and WT genes.
NativeDS = join(End3_Contr,End5_Contr,'Keys','ORF');
NativeDS = join(NativeDS,FitFlow_Fast_Slow,'Keys','ORF');
NativeDS = join(NativeDS,ORFseqDS,'Keys','ORF');

NativeDS = join(NativeDS,End5_Contr_adds,'Keys','ORF');
NativeDS = join(NativeDS,UPF1del_Dataset,'Keys','ORF');
NativeDS = join(NativeDS,DS_cor,'Keys','ORF');


%% Add missing features to the NativeDS:

%Correct Seq:
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq{I}; NativeDS.ORF_seq(I) = ORF;
end
%GC_content and length:
NativeDS.GC_content = 0;
NativeDS.ORF_length = 0;
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq(I); ORF = ORF{1};
    NativeDS.GC_content(I) = (length(regexp(ORF,'C')) + length(regexp(ORF,'G')))./length(ORF);
    NativeDS.ORF_length(I) = length(ORF);
end
%%
%tAI:
NativeDS.tAI = cellfun(@calc_tAI , NativeDS.ORF_seq);
%CAI:
    %[~,CAI]=CodonAdaptationIndex( 'ATG');
    %NativeDS.CAI = arrayfun(@(I)CodonAdaptationIndex( NativeDS.Seq{I}, CAI ) , 1:length(NativeDS)  )' ;
%Sequence feats:
NativeDS.GAAAGA = 0;
NativeDS.ACGTTA = 0;
NativeDS.AG1 = 0;
NativeDS.AC2 = 0;
NativeDS.A2 = 0;
for I = 1:length(NativeDS.TranscriptSequence)  
    %GAAAGA    
    NativeDS.GAAAGA(I) = length(regexp(NativeDS.TranscriptSequence{I},'GAAAGA'));
    %ACGTTA
    NativeDS.ACGTTA(I) = length(regexp(NativeDS.TranscriptSequence{I},'ACGTTA'));
    %position in the transcript:
    AG1 = 0;
    AC2 = 0;
    A2 = 0;
    ORF = NativeDS.ORF_seq{I};
    Cp = 1; %Codon Position
    Cn = 1; %Codon Number
    NumberCodons = length(ORF)/3;
    while Cn <= NumberCodons
        Codon = ORF(Cp:Cp+2);  
        if strcmp(Codon(1),'A') == 1 | strcmp(Codon(1),'G') == 1; AG1 = AG1 + 1; end
        if strcmp(Codon(2),'A') == 1 | strcmp(Codon(2),'C') == 1; AC2 = AC2 + 1; end
        if strcmp(Codon(2),'A') == 1; A2 = A2 + 1; end
        Cp = Cp + 3;
        Cn = Cn + 1;
    end
    NativeDS.AG1(I) = AG1;
    NativeDS.AC2(I) = AC2;
    NativeDS.A2(I) = A2;
end
%Convert To Log The Expression:
NativeDS.log2expression = log2(NativeDS.WTexpression);
NativeDS.log2end3_YFP = log2(NativeDS.End3_YFP);
NativeDS.log2end5_YFP = log2(NativeDS.End5_YFP_glucose);

%% Add Nag and Yas Expression

idxNag = find(strcmp(SRAs_DS.SRP_accession,'SRP000227'));
ERRsNAG = SRAs_DS.SRR_accession(idxNag)';
idxNAG_DS = [];
for ERR = ERRsNAG
    idxNAG_DS = [idxNAG_DS ; find(strcmp(NativeDS.Properties.VarNames,ERR))];
end
NativeDS.NAG = log2(mean(double(NativeDS(:,idxNAG_DS))')'.*1e9);
idxYas = find(strcmp(SRAs_DS.SRP_accession,'SRP002790'));
ERRsYAS = SRAs_DS.SRR_accession(idxYas)';
idxYAS_DS = [];
for ERR = ERRsYAS
    idxYAS_DS = [idxYAS_DS ; find(strcmp(NativeDS.Properties.VarNames,ERR))];
end
NativeDS.YAS = log2(mean(double(NativeDS(:,idxYAS_DS))')'.*1e9);

%% Add all ORF features (>4000) to NativeDS:

%create a dataset with the features that Julia Calculated
%FDS = gDNA_BaseCombinationCount(NativeDS, 'ORF_seq', {@basecount @dimercount @codoncount @aacount});
%FDS = SeqCodonPositionNucFreq(NativeDS.ORF_seq, FDS);
%FDS = SeqTrimerFreq(NativeDS.ORF_seq, FDS);
%FDS = SeqHexamerFreq(NativeDS.ORF_seq, FDS); %this takes the longest time
% join with Native DataSet:
%NativeDS = horzcat(NativeDS,FDS);
%% Add missing features to the gDNA library Dataset (WT or WT_2014):

%Sequence feats:
WT_RPL4.GAAAGA = 0;
WT_RPL4.ACGTTA = 0;
WT_RPL4.AG1 = 0;
WT_RPL4.AC2 = 0;
WT_RPL4.A2 = 0;
for I = 1:length(WT_RPL4.trans_seq)  
    %GAAAGA    
    WT_RPL4.GAAAGA(I) = length(regexp(WT_RPL4.trans_seq{I},'GAAAGA'));
    %ACGTTA
    WT_RPL4.ACGTTA(I) = length(regexp(WT_RPL4.trans_seq{I},'ACGTTA'));
    %position in the transcript:
    AG1 = 0;
    AC2 = 0;
    A2 = 0;
    ORF = WT_RPL4.orf_seq{I};
    Cp = 1; %Codon Position
    Cn = 1; %Codon Number
    NumberCodons = length(ORF)/3;
    while Cn <= NumberCodons
        Codon = ORF(Cp:Cp+2);  
        if strcmp(Codon(1),'A') == 1 | strcmp(Codon(1),'G') == 1; AG1 = AG1 + 1; end
        if strcmp(Codon(2),'A') == 1 | strcmp(Codon(2),'C') == 1; AC2 = AC2 + 1; end
        if strcmp(Codon(2),'A') == 1; A2 = A2 + 1; end
        Cp = Cp + 3;
        Cn = Cn + 1;
    end
    WT_RPL4.AG1(I) = AG1;
    WT_RPL4.AC2(I) = AC2;
    WT_RPL4.A2(I) = A2;
end
%%
%Sequence feats:
WT.GAAAGA = 0;
WT.ACGTTA = 0;
WT.AG1 = 0;
WT.AC2 = 0;
WT.A2 = 0;
for I = 1:length(WT.trans_seq)  
    %GAAAGA    
    WT.GAAAGA(I) = length(regexp(WT.trans_seq{I},'GAAAGA'));
    %ACGTTA
    WT.ACGTTA(I) = length(regexp(WT.trans_seq{I},'ACGTTA'));
    %position in the transcript:
    AG1 = 0;
    AC2 = 0;
    A2 = 0;
    ORF = WT.orf_seq{I};
    Cp = 1; %Codon Position
    Cn = 1; %Codon Number
    NumberCodons = length(ORF)/3;
    while Cn <= NumberCodons
        Codon = ORF(Cp:Cp+2);  
        if strcmp(Codon(1),'A') == 1 | strcmp(Codon(1),'G') == 1; AG1 = AG1 + 1; end
        if strcmp(Codon(2),'A') == 1 | strcmp(Codon(2),'C') == 1; AC2 = AC2 + 1; end
        if strcmp(Codon(2),'A') == 1; A2 = A2 + 1; end
        Cp = Cp + 3;
        Cn = Cn + 1;
    end
    WT.AG1(I) = AG1;
    WT.AC2(I) = AC2;
    WT.A2(I) = A2;
end

%Sequence feats:
WT_2014.GAAAGA = 0;
WT_2014.ACGTTA = 0;
WT_2014.AG1 = 0;
WT_2014.AC2 = 0;
WT_2014.A2 = 0;
for I = 1:length(WT_2014.trans_seq)  
    %GAAAGA    
    WT_2014.GAAAGA(I) = length(regexp(WT_2014.trans_seq{I},'GAAAGA'))./(WT_2014.trans_seq{I}./6);
    %ACGTTA
    WT_2014.ACGTTA(I) = length(regexp(WT_2014.trans_seq{I},'ACGTTA'));
    %position in the transcript:
    AG1 = 0;
    AC2 = 0;
    A2 = 0;
    ORF = WT_2014.orf_seq{I};
    Cp = 1; %Codon Position
    Cn = 1; %Codon Number
    NumberCodons = length(ORF)/3;
    while Cn <= NumberCodons
        Codon = ORF(Cp:Cp+2);  
        if strcmp(Codon(1),'A') == 1 | strcmp(Codon(1),'G') == 1; AG1 = AG1 + 1; end
        if strcmp(Codon(2),'A') == 1 | strcmp(Codon(2),'C') == 1; AC2 = AC2 + 1; end
        if strcmp(Codon(2),'A') == 1; A2 = A2 + 1; end
        Cp = Cp + 3;
        Cn = Cn + 1;
    end
    WT_2014.AG1(I) = AG1;
    WT_2014.AC2(I) = AC2;
    WT_2014.A2(I) = A2;
end
%%
% Add nTE:
WT_RPL4.nte = cellfun(@(X)calc_tAI(X,'nTE') , WT_RPL4.orf_seq);




%% Contribution of each of the regions to NATIVE EXPRESSION expression (current Fig1A and 1D):

idx = find(NativeDS.NAG>-Inf);
All_DS = NativeDS(idx,:);

%AllORFFeats = All_DS.Properties.VarNames(69:end);
predictors = {'log2end5_YFP'}; Ylabel = 'Predicted (5'')' ;
%predictors = {'log2end3_YFP'}; Ylabel = 'Predicted (3'')' ;
%predictors = {'log2end3_YFP','log2end5_YFP'}; Ylabel = 'Predicted (5'' + 3'')' ;
%predictors = {'GC_content','ORF_length','tAI','CAI_left','GAAAGA','ACGTTA','AG1','AC2','A2'}; Ylabel = 'Predicted (ORF)';

%predictors = {'GC_content','ORF_length','tAI','CAI','GAAAGA','ACGTTA','AG1','AC2','A2','UTR3length'}; Ylabel = 'Predicted from ORF';
%predictors = {'GC_content','ORF_length','tAI','CAI_left','GAAAGA','ACGTTA','AG1','AC2','A2','log2end5_YFP'}; Ylabel = 'Predicted (5'' + ORF '')';
%predictors = {'GC_content','ORF_length','tAI','CAI','orf_GAAAGA','orf_ACGTTA','AG1','AC2','A2','UTR3length','log2end3_YFP','log2end5_YFP'}; Ylabel = 'Predicted (5'' + ORF + 3'')';
%predictors = horzcat({'GC_content','ORF_length','tAI','CAI','UTR3length','log2end3_YFP','log2end5_YFP'},AllORFFeats); Ylabel = 'Predicted (5'' + ORF + 3''), Overfit';

All_DS.Predicted_5prime = 0;
%All_DS.Predicted_orf = 0;

%Run the prediction 50 times and take the best one
r2 = 0;
for Pred = 1
    Pred
    %Run 10 fold crossvalidation for running the prediction, split the sequences: by >0.75 and <0.75 CAI
    MeasuredMat = [];
    PredictedMat = [];
    for CAIsplit = 1:2
        %Splitting:
        %if CAIsplit == 1; idx = find(All_DS.CAI < 0.75); end
        %if CAIsplit == 2; idx = find(All_DS.CAI >= 0.75); end
        %Without Splitting:
        if CAIsplit == 1; idx = find(All_DS.CAI_left < 10); end
        if CAIsplit == 2; idx = find(All_DS.CAI_left > 10); end
        if isempty(idx) == 0
            DS = All_DS(idx,:);
            K = 10;
            idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
            for I = 1:K
                I
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas
                
                Predicted = predict(mdl,DS(idx_test,:));
                %Measured = DS.log2expression(idx_test);
                Measured = DS.NAG(idx_test);                
                PredictedMat = [PredictedMat ; Predicted];
                MeasuredMat = [MeasuredMat ; Measured]; 
                All_DS.Predicted_5prime(idx_test) = Predicted;
                %All_DS.Predicted_orf(idx_test) = Predicted;
            end
        end
    end
    r = corr(MeasuredMat,PredictedMat,'type','Spearman');
    if r^2 > r2
        X = MeasuredMat; 
        Y = PredictedMat; 
        r2 = r^2;
                
    end
end

Gray = [.7 .7 .7];
fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on;
plot(linspace(min([X;Y]),max([X;Y]),3),linspace(min([X;Y]),max([X;Y]),3),'-','color',[.7 .7 .7])

%plot(X,Y,'o','markersize',4,'markerfacecolor',[0,0.4,0.8],'color','k')
%plot(X,Y,'o','markersize',9,'markerfacecolor','k','color','k')
%plot(X,Y,'o','markersize',6,'color','k')
plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','none','linewidth',0.7)

%plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)
xlabel('Measured expression');
ylabel(Ylabel);
xlim([min([X;Y]),max([X;Y])])
ylim([min([X;Y]),max([X;Y])])
text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r^2)),'fontsize',10)
set(gca,'fontsize',11)

%% Add predictors to dataset:



%% model with interaction

predictors = {'Predicted_orf','Predicted_5prime'}; Ylabel = 'Predicted (5'' + ORF '')';

%Run the prediction 50 times and take the best one
r2 = 0;
for Pred = 1
    Pred
    %Run 10 fold crossvalidation for running the prediction, split the sequences: by >0.75 and <0.75 CAI
    MeasuredMat = [];
    PredictedMat = [];
    for CAIsplit = 1:2
        %Splitting:
        %if CAIsplit == 1; idx = find(All_DS.CAI < 0.75); end
        %if CAIsplit == 2; idx = find(All_DS.CAI >= 0.75); end
        %Without Splitting:
        if CAIsplit == 1; idx = find(All_DS.CAI_left < 10); end
        if CAIsplit == 2; idx = find(All_DS.CAI_left > 10); end
        if isempty(idx) == 0
            DS = All_DS(idx,:);
            K = 2;
            idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
            for I = 1:K
                I
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                mdl = LinearModel.fit(DS(idx_train,:),'interactions', 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas
                
                Predicted = predict(mdl,DS(idx_test,:));
                %Measured = DS.log2expression(idx_test);
                Measured = DS.NAG(idx_test);                
                PredictedMat = [PredictedMat ; Predicted];
                MeasuredMat = [MeasuredMat ; Measured];
            end
        end
    end
    r = corr(MeasuredMat,PredictedMat,'type','Spearman');
    if r^2 > r2
        X = MeasuredMat; 
        Y = PredictedMat; 
        r2 = r^2;
    end
end

Gray = [.7 .7 .7];
fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on;
plot(linspace(min([X;Y]),max([X;Y]),3),linspace(min([X;Y]),max([X;Y]),3),'-','color',[.7 .7 .7])

%plot(X,Y,'o','markersize',4,'markerfacecolor',[0,0.4,0.8],'color','k')
%plot(X,Y,'o','markersize',9,'markerfacecolor','k','color','k')
%plot(X,Y,'o','markersize',6,'color','k')
plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','none','linewidth',0.7)

%plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)
xlabel('Measured expression');
ylabel(Ylabel);
xlim([min([X;Y]),max([X;Y])])
ylim([min([X;Y]),max([X;Y])])
text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r^2)),'fontsize',10)
set(gca,'fontsize',11)


%% Correlation between each prediction:

idx = find(NativeDS.NAG>-Inf);
All_DS = NativeDS(idx,:);

predictors1 = {'GC_content','ORF_length','tAI','CAI_left','GAAAGA','ACGTTA','AG1','AC2','A2'}; Ylabel = 'Predicted (ORF)';
predictors2 = {'log2end3_YFP','log2end5_YFP'};

for TypeP = 1:2
    if TypeP ==1
        predictors = predictors1;
    end
    if TypeP ==2
        predictors = predictors2;
    end
      
%Run the prediction 50 times and take the best one
r2 = 0;
for Pred = 1:3
    Pred
    %Run 10 fold crossvalidation for running the prediction, split the sequences: by >0.75 and <0.75 CAI
    MeasuredMat = [];
    PredictedMat = [];
    for CAIsplit = 1:2
        %Splitting:
        %if CAIsplit == 1; idx = find(All_DS.CAI < 0.75); end
        %if CAIsplit == 2; idx = find(All_DS.CAI >= 0.75); end
        %Without Splitting:
        if CAIsplit == 1; idx = find(All_DS.CAI_left < 10); end
        if CAIsplit == 2; idx = find(All_DS.CAI_left > 10); end
        if isempty(idx) == 0
            DS = All_DS(idx,:);
            K = 10;
            idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
            for I = 1:K
                I
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas
                
                Predicted = predict(mdl,DS(idx_test,:));
                %Measured = DS.log2expression(idx_test);
                Measured = DS.NAG(idx_test);                
                PredictedMat = [PredictedMat ; Predicted];
                MeasuredMat = [MeasuredMat ; Measured];
            end
        end
    end
    r = corr(MeasuredMat,PredictedMat,'type','Spearman');
    if r^2 > r2
        X = MeasuredMat; 
        Y = PredictedMat; 
        r2 = r^2;
    end
end
if TypeP==1; X = PredictedMat; end
if TypeP==2; Y = PredictedMat; end

end

Gray = [.7 .7 .7]
fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on;
plot(linspace(min([X;Y]),max([X;Y]),3),linspace(min([X;Y]),max([X;Y]),3),'-','color',[.7 .7 .7])

%plot(X,Y,'o','markersize',4,'markerfacecolor',[0,0.4,0.8],'color','k')
%plot(X,Y,'o','markersize',9,'markerfacecolor','k','color','k')
%plot(X,Y,'o','markersize',6,'color','k')
plot(X,Y,'o','markersize',6,'markerfacecolor','k','color','none','linewidth',0.7)

%plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)
xlabel('Predicted (ORF)');
ylabel('Predicted (5''+3'')');
xlim([min([X;Y]),max([X;Y])])
ylim([min([X;Y]),max([X;Y])])
text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r^2)),'fontsize',10)
set(gca,'fontsize',11)

%% Contribution of the ORF to expression (RandomRNA lib) (current 1C)

All_DS = WT_2014;
%All_DS = Julia_DS;

COLORS = [[39,101,186]./255; [0.7 0.7 0.7]];

%predictors = {'orf_GC','insert_l','tAI','orf_stop_codon_loc','utr3_l','premature_stop','GAAAGA','ACGTTA','AG1','AC2','A2'};
predictors = {'orf_GC','insert_l','tAI','utr3_l','AG1','AC2','A2'};

%All_DS = WT_RPL4; 
%predictors = {'insert_GC','insert_l','tAI','CAI','GAAAGA','ACGTTA','AG1','AC2','A2','orf_pct_stop_loc','utr3_l','premature_stop'};
%predictors = { 'utr3_l'};

fh = figure('units','centimeters','position',[5 5 5.5 5.5]); 
hold on;

for D = 1:2
    D
    if D == 1
        idx_len = All_DS.insert_l <= 400 & All_DS.insert_l >= 100 ; 
        colors = flipud(parula(100));
        colors = colors(5:50,:);


    end
    if D == 2
        idx_len = All_DS.insert_l <= 400 & All_DS.insert_l >= 100  & All_DS.premature_stop==0;
        colors = flipud(bone(100));
        colors = colors(5:50,:);
    end
    A = length(find(idx_len))
    DS = All_DS(idx_len,:);

    %Run the prediction 50 times and take the best one
    r2 = 0;
    for Pred = 1
        %Run 10 fold crossvalidation for running the prediction, split the sequences: by >0.75 and <0.75 CAI
        MeasuredMat = [];
        PredictedMat = [];
        for CAIsplit = 1
            %Splitting
            %if CAIsplit == 1; idx = find(All_DS.CAI < 0.75); end
            %if CAIsplit == 2; idx = find(All_DS.CAI >= 0.75); end
            %Without Splitting:
            if CAIsplit == 1; idx = find(DS.CAI < 10); end
            if CAIsplit == 2; idx = find(DS.CAI > 10); end
            if isempty(idx) == 0
                ds = DS(idx,:);
                K = 10;
                idx_cv = crossvalind('Kfold', length(ds.Expression), K);
                for I = 1:K
                    idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                    %mdl = LinearModel.fit(  DS(idx_train,:), 'ResponseVar', 'RPL4', 'PredictorVars', predictors ); TypeData = 'RPL4A_P_R';
                    mdl = LinearModel.fit(ds(idx_train,:), 'ResponseVar', 'Expression', 'PredictorVars', predictors ); TypeData = 'GALL_P_R';
                    Predicted = predict(mdl,ds(idx_test,:));
                    Measured = ds.Expression(idx_test);
                    PredictedMat = [PredictedMat ; Predicted];
                    MeasuredMat = [MeasuredMat ; Measured];
                end
            end
        end
        %r = corr(MeasuredMat,PredictedMat,'type','Spearman');
        r = corr(MeasuredMat,PredictedMat);
        ds.premature_stop
        if r^2 > r2
            X = MeasuredMat; 
            Y = PredictedMat; 
            r2 = r^2;
        end
    end
    plot(linspace(-15,5,3),linspace(-15,5,3),'-','color',[.7 .7 .7])
    %plot(X,Y,'o','markersize',1,'markerfacecolor',COLORS(D,:),'color','none')
    %plot(X,Y,'o','markersize',2,'markerfacecolor',[1,0,0.5],'linewidth',2)
    %dscatter(X,Y)
    %dscatter(X,Y)
    %colormap(colors)
    scatter(X,Y,'o','markerfacecolor',COLORS(D,:),'markeredgecolor',COLORS(D,:),'markerfacealpha',0.2,'markeredgealpha',0.1)
    if D == 2
        %scatter(X,Y,'o','markerfacecolor',COLORS(D,:),'markeredgecolor',COLORS(D,:),'markerfacealpha',1,'markeredgealpha',1,'SizeData',10)
    end

end
xlim([-15,5])
ylim([-15,5])

xlabel('Measured expression');
ylabel('Predicted (ORF)');
%title(TypeData)
%text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r^2)),'fontsize',11)
set(gca,'xticklabel' , { '10^{-3}'  '10^{-2}'  '10^{-1}'  '10^{0}'  '10^{1}'} ,'fontsize',11,'yticklabel' , { '10^{-3}'  '10^{-2}'  '10^{-1}'  '10^{0}'  '10^{1}'} )

%% See if jULIA'S features correlate with mine

predictors = {'orf_GC','insert_l','tAI','CAI','orf_stop_codon_loc','utr3_l','premature_stop','GAAAGA','ACGTTA','AG1','AC2','A2'};
figure; hold on;

X = Julia_DS.A2;
Y = WT_2014.A2;
plot(X,Y,'.')



%% Plot the expression distributions for the GALL library
% mod by LBC June 08

All_DS = WT_2014;
fh = figure('units','centimeters','position',[5 5 8 5 ]);
hold on;
COLORS = [[39,101,186]./255; [0.7 0.7 0.7]];


for type = 1:2
    if type == 2
        idx_len = All_DS.insert_l <= 400 & All_DS.insert_l >= 100 & All_DS.premature_stop == 0;% log10(All_DS.DNA)>0 & All_DS.premature_stop==0; 
    end
    if type == 1
        idx_len = All_DS.insert_l <= 400 & All_DS.insert_l >= 100 ;% log10(All_DS.DNA)>0 & All_DS.premature_stop==0; 
        
    end
    DS = All_DS(idx_len,:);



    X =  DS.Expression  ;
    %  X = All_DS.Expression( ~All_DS.premature_stop );

    bins = linspace( min(X) , max(X) , 50 ) ;

    histogram( X , bins , 'FaceColor' , COLORS(type,:) , 'EdgeAlpha', 0.5 ,'FaceAlpha',1)
    ylabel('Number of inserts');
    set(gca,'fontsize',11)
    xtickl = [0.001 0.01 .1 1 10] ;
    set(gca,'xtick', log2(xtickl) )
    set(gca,'xticklabel' , { '10^{-3}'  '10^{-2}'  '10^{-1}'  '10^{0}'  '10^{1}'} )
    set(gca,'ytick',0:50:1000)
    axis tight ;
    xlim([-11 4.5])
    
end

%% Plot the weights of several predictors in the GALL vs RPL4A libraries:

All_DS = WT_RPL4; 

%zscore predictors:
All_DS.utr3_l_p = zscore(All_DS.utr3_l);
All_DS.CAI_p = zscore(All_DS.CAI);
All_DS.nte_p = zscore(All_DS.nte);
All_DS.tAI_p = zscore(All_DS.tAI);
All_DS.insert_GC_p = zscore(All_DS.insert_GC);
All_DS.insert_l_p = zscore(All_DS.insert_l);

%predictors = { 'utr3_l_p','CAI_p','tAI_p','insert_GC_p','insert_l_p'};
%Labels = {'3''UTR','CAI','tAI','GC','length'};
%Xtick = [1.25,2.25,3.25,4.25,5.25];

predictors = { 'utr3_l_p','tAI_p'};
Labels = {'3''UTR','tAI'};
Xtick = [1.25,2.25,3.25,4.25,5.25];

PTC = [1;0;0;0;0];

idx_len = All_DS.insert_l <= 400 & All_DS.insert_l >= 100 & log10(All_DS.DNA)>0; 
All_DS = All_DS(idx_len,:);

%generate estimates
Ymat = [];
Gmat = [];
P = 1;
for predictor = predictors
    P
    DS = All_DS(All_DS.premature_stop==PTC(P),:);
    %run the model for both libraries:
    for L = [0,0.5]
        for Pred = 1
            MeasuredMat = [];
            PredictedMat = [];
            K = 10;
            idx_cv = crossvalind('Kfold', length(DS.Expression), K);
            for I = 1:K
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                if L == 0
                    mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'Expression', 'PredictorVars', predictor ); TypeData = 'GALL_P_R';
                end
                if L == 0.5
                    mdl = LinearModel.fit(  DS(idx_train,:), 'ResponseVar', 'RPL4', 'PredictorVars', predictor ); TypeData = 'RPL4A_P_R';                    
                end   
                Estimate = mdl.Coefficients.Estimate(2);
                Ymat = [Ymat ; Estimate];
                Gmat = [Gmat ; (P+L)];
            end
        end
    end
    P = P + 1;
end


%plot results
Red = [0.8,0,0]; 
Blue = [0,0,0.8];
colors = [Red;Blue;Red;Blue;Red;Blue;Red;Blue;Red;Blue];
fh = figure('units','centimeters','position',[5 5 8 5]); 

Gvals = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5];
I = 1;
%plot jittered values:
for Gval = Gvals
    hold on;
    G = Gmat(Gmat == Gval);
    Y = Ymat(Gmat == Gval);
    jitterAmount = 0.15;
    jitterValuesX = 2*(rand(size(G))-0.5)*jitterAmount;
    X = G + jitterValuesX;
    plot(X,Y,'o','color',colors(I,:),'markerfacecolor',colors(I,:),'markersize',3)
    I = I + 1;
end

%boxplot:
b = boxplot(Ymat,Gmat,'color','k','positions',Gmat,'width',0.6,'notch','on');
set(b(7,:),'Visible','off')

set(gca,'xtick',Xtick,'xticklabels',Labels,'fontsize',11)
ylabel('Linear model coefficient')
ylim([-1,0.5])