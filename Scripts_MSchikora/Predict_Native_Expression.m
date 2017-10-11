%% Description
%Predict Native expression with and without interactions
%Look for differences in prediction in native genes

%% load data
clear all

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/');
%WT RNAseq expression. Column 4:
FitFlow_Fast_Slow = dataset('file','./vanDijk15/FitFlow__Fast_Slow.tab','delimiter','\t'); 
FitFlow_Fast_Slow.Properties.VarNames = {'ORF','Gene','Location','WTexpression','x18_2308'};
%3'end contribution to expression:
End3_Contr = dataset('file','./Shalem13/Shalem13.tab','delimiter','\t'); 
End3_Contr.Properties.VarNames = {'ORF','End3_YFP','ClonedIntergenicSequence'};
%5'end contribution to expression (use glucose):
End5_Contr = dataset('file','./Keren13/Keren13_Expression.tab','delimiter','\t'); 
End5_Contr.Properties.VarNames(3) = {'End5_YFP_glucose'};

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%All integrated gDNA random RNA library datasets, cointains many datasets:
load('Integrated_DataSets_2014_Adds.mat');
load('ORFtoSeq.mat');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/LBC/2016_September__PromoterInfluencesRNALife__CodonBias_Predicts_Keren_Residuals/')
%DataSet of Keren (5'end contribution) with some adds:
End5_Contr_adds = dataset('file','keren_rnaseq_cai.tab','Delimiter','\t');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Data/');
%UPF1 mutant data --> cointains the information for the 3'utr length:
load('./External Data/S_cer_UPF1_mutant_RNAseq.mat');
UPF1del_Dataset = G;

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%curated DEE datasets:
load('All_RNA_seq_S_cer_integrated_All_DS.mat')

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/DigitalExpressionExplorer/')
%SRA_S dataset:
SRAs_DS = dataset('file','./metadata.tab','Delimiter','\t');  

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%TATA Types:
load('Promoter_Features.mat')
%TFIID/SAGA Types:
load('Yeast_SAGAvsTFIID.mat')


%% merge datasets

%A dataset that cointains the information for all the native genes that are in each of the 5', 3' and WT genes.
NativeDS = join(End5_Contr,FitFlow_Fast_Slow,'Keys','ORF');
NativeDS = join(End5_Contr_adds,NativeDS,'Keys','ORF');
NativeDS = join(NativeDS,DS_cor,'Keys','ORF');
NativeDS = join(NativeDS,ORFseqDS,'Keys','ORF');
% Delete genes that are not in the UPF1del datasets:
BadIDX = []; 
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.ORF(I);
    if isempty(find(strcmp(UPF1del_Dataset.ORF,ORF))) == 1
        BadIDX = [BadIDX ; I];  
    end
end
NativeDS(BadIDX,:) = [];
%Add the UPF1_del dataset:
NativeDS = join(NativeDS,UPF1del_Dataset,'Keys','ORF');


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
%FDS = SeqHexamerFreq(NativeDS.ORF_seq, FDS); %this takes the longest time,
% join with Native DataSet:
%NativeDS = horzcat(NativeDS,FDS);

%% Add 5' and ORF contributions:

idx = find(NativeDS.NAG>-Inf);
All_DS = NativeDS(idx,:);

predictors = {'GC_content','ORF_length','tAI','GAAAGA','ACGTTA','AG1','AC2','A2'}; 

All_DS.Predicted_5prime = 0;
All_DS.Predicted_orf = 0;

for Region = 1:2
    if Region == 1
        predictors = {'log2end5_YFP'}; %5'
    end
    if Region == 2
        predictors = {'GC_content','ORF_length','tAI','CAI_left','GAAAGA','ACGTTA','AG1','AC2','A2'}; %3'
    end   
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
                K = 5;
                idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
                for I = 1:K
                    I
                    idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                    mdl = LinearModel.fit(DS(idx_train,:), 'linear','ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas

                    Predicted = predict(mdl,DS(idx_test,:));
                    %Measured = DS.log2expression(idx_test);
                    Measured = DS.NAG(idx_test);                
                    PredictedMat = [PredictedMat ; Predicted];
                    MeasuredMat = [MeasuredMat ; Measured]; 
                    if Region == 1
                        All_DS.Predicted_5prime(idx_test) = Predicted;
                    end
                    if Region == 2
                        All_DS.Predicted_orf(idx_test) = Predicted;
                    end                    
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
end

%% Avaluate each of the contributions:
predictors = {'Predicted_orf','Predicted_5prime'}; Ylabel = 'Predicted (5'' + ORF '' Directly)';
%predictors = {'GC_content','ORF_length','tAI','CAI_left','GAAAGA','ACGTTA','AG1','AC2','A2','log2end5_YFP'}; Ylabel = 'Predicted (5'' + ORF '' All feats)';

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
            K = 5;
            idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
            for I = 1:K
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                mdl = LinearModel.fit(DS(idx_train,:),'interactions', 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                
                %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas
                
                Predicted = predict(mdl,DS(idx_test,:));
                %Measured = DS.log2expression(idx_test);
                Measured = DS.NAG(idx_test);                
                PredictedMat = [PredictedMat ; Predicted(Predicted > -Inf & Measured > -Inf)];
                MeasuredMat = [MeasuredMat ; Measured(Predicted > -Inf & Measured > -Inf)];
            end
        end
    end
    r = corr(MeasuredMat,PredictedMat,'type','Spearman')
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
r2 = rsquare(X,Y)
text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r2)),'fontsize',10)
set(gca,'fontsize',11)

%% Do the analysis for TATA and TATA-less genes

Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
%predictors = {'GC_content','ORF_length','tAI','GAAAGA','ACGTTA','AG1','AC2','A2'}; Ylabel = 'Predicted (5'' + ORF '' All feats)';
predictors = {'GC_content','ORF_length','tAI','AG1','AC2','A2'}; Ylabel = 'Predicted (ORF)';


for TypeP = 1:2
    
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end

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
            if CAIsplit == 1; idx = find(All_DS.CAI_left < 10 & IDXtata & All_DS.NAG>0); end
            if CAIsplit == 2; idx = find(All_DS.CAI_left > 10 & IDXtata & All_DS.NAG>0); end
            if isempty(idx) == 0
                %Randomized_for_idx = randperm(length(idx),150);
                %idx = idx(Randomized_for_idx);
                DS = All_DS(idx,:);
                K = 5;
                idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
                for I = 1:K
                    idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                    mdl = LinearModel.fit(DS(idx_train,:),'linear', 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy

                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas

                    Predicted = predict(mdl,DS(idx_test,:));
                    %Measured = DS.log2expression(idx_test);
                    Measured = DS.NAG(idx_test);                
                    PredictedMat = [PredictedMat ; Predicted(Predicted > -Inf & Measured > -Inf)];
                    MeasuredMat = [MeasuredMat ; Measured(Predicted > -Inf & Measured > -Inf)];
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
    plot(X,Y,'o','markersize',6,'markerfacecolor',colors(TypeP,:),'color','none','linewidth',0.7)

    %plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)
    xlabel('Measured expression');
    ylabel(Ylabel);
    xlim([min([X;Y]),max([X;Y])])
    ylim([min([X;Y]),max([X;Y])])
    r2 = rsquare(X,Y);
    text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r2)),'fontsize',10)
    set(gca,'fontsize',11,'xtick',[0,5,10])

end

%% Boxplot the residual at different expression levels

Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
%predictors = {'GC_content','ORF_length','tAI','GAAAGA','ACGTTA','AG1','AC2','A2'}; Ylabel = 'Predicted (5'' + ORF '' All feats)';
predictors = {'GC_content','ORF_length','tAI','AG1','AC2','A2'}; Ylabel = 'Predicted (ORF)';

fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on;

for TypeP = 1:2
    
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end

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
            if CAIsplit == 1; idx = find(All_DS.CAI_left < 10 & IDXtata & All_DS.NAG>0); end
            if CAIsplit == 2; idx = find(All_DS.CAI_left > 10 & IDXtata & All_DS.NAG>0); end
            if isempty(idx) == 0
                %subsample
                %Randomized_for_idx = randperm(length(idx),150);
                %idx = idx(Randomized_for_idx);
                DS = All_DS(idx,:);
                K = 2;
                idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
                for I = 1:K
                    idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                    mdl = LinearModel.fit(DS(idx_train,:),'linear', 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy

                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas

                    Predicted = predict(mdl,DS(idx_test,:));
                    %Measured = DS.log2expression(idx_test);
                    Measured = DS.NAG(idx_test);                
                    PredictedMat = [PredictedMat ; Predicted(Predicted > -Inf & Measured > -Inf)];
                    MeasuredMat = [MeasuredMat ; Measured(Predicted > -Inf & Measured > -Inf)];
                    
                    [a,b] = count_unique(round(DS.NAG(idx_train)));
                    %b = boxplot( mdl.Residuals.Raw , round(DS.NAG(idx_train)),'notch','off','color',colors(TypeP,:));
                    %set(b(7,:),'Visible','off');
                    plot( a , b,'.','color',colors(TypeP,:),'markersize',15);
                    
                    
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
    %plot(linspace(min([X;Y]),max([X;Y]),3),linspace(min([X;Y]),max([X;Y]),3),'-','color',[.7 .7 .7])

    %plot(X,Y,'o','markersize',4,'markerfacecolor',[0,0.4,0.8],'color','k')
    %plot(X,Y,'o','markersize',9,'markerfacecolor','k','color','k')
    %plot(X,Y,'o','markersize',6,'color','k')
    %plot(X,Y,'o','markersize',6,'markerfacecolor',colors(TypeP,:),'color','none','linewidth',0.7)

    %plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)

end

    xlabel('Measured expression');
    ylabel('Residual of the fit');
    %xlim([min([X;Y]),max([X;Y])])
    %ylim([-10,5])
    r2 = rsquare(X,Y);
    text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r2)),'fontsize',10)
    %set(gca,'fontsize',11,'xtick',[0,5,10])


%% ecdf

Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on;

for TypeP = 1:2
    
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end
    
    [F,X] = ksdensity(All_DS.YAS(IDXtata));
    plot(X,F,'-','color',colors(TypeP,:),'linewidth',2)

end

xlabel('Expression')


%% Generte 100 fits for confidence intervals for TATA / TATA-less

predictors = {'GC_content','ORF_length','tAI','GAAAGA','ACGTTA','AG1','AC2','A2'}; Ylabel = 'Predicted (5'' + ORF '' All feats)';

PtypeMat = [];
r2Mat = [];
for Run = 1:10
Run
for TypeP = 1:2  
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end

    r2 = 0;
    for Pred = 1
        %Run 10 fold crossvalidation for running the prediction, split the sequences: by >0.75 and <0.75 CAI
        MeasuredMat = [];
        PredictedMat = [];
        for CAIsplit = 1:2
            %Splitting:
            %if CAIsplit == 1; idx = find(All_DS.CAI < 0.75); end
            %if CAIsplit == 2; idx = find(All_DS.CAI >= 0.75); end
            %Without Splitting:
            if CAIsplit == 1; idx = find(All_DS.CAI_left < 10 & IDXtata); end
            if CAIsplit == 2; idx = find(All_DS.CAI_left > 10 & IDXtata); end
            if isempty(idx) == 0
                Randomized_for_idx = randperm(length(idx),225);
                idx = idx(Randomized_for_idx);
                DS = All_DS(idx,:);
                K = 5;
                idx_cv = crossvalind('Kfold', length(DS.WTexpression), K);
                for I = 1:K
                    idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'log2expression', 'PredictorVars', predictors ); %FitFlow
                    mdl = LinearModel.fit(DS(idx_train,:),'linear', 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy
                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'NAG', 'PredictorVars', predictors ); %Nagalaschmy

                    %mdl = LinearModel.fit(DS(idx_train,:), 'ResponseVar', 'YAS', 'PredictorVars', predictors ); %Yas

                    Predicted = predict(mdl,DS(idx_test,:));
                    %Measured = DS.log2expression(idx_test);
                    Measured = DS.NAG(idx_test);                
                    PredictedMat = [PredictedMat ; Predicted(Predicted > -Inf & Measured > -Inf)];
                    MeasuredMat = [MeasuredMat ; Measured(Predicted > -Inf & Measured > -Inf)];
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

    r2 = rsquare(X,Y);
    PtypeMat = [PtypeMat ; TypeP];
    r2Mat = [r2Mat ;  r2];
end
end

% Calculate CI for the differences:

%TATA:
TATAvals = r2Mat(PtypeMat == 1);
[CI_tata] = bootci(1000, @mean, TATAvals)

%TATA-less:
TATAlessvals = r2Mat(PtypeMat == 2);
[CI_tataless] = bootci(1000, @mean, TATAlessvals)