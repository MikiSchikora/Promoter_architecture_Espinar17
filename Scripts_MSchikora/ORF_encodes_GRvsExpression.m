%% Description
%This plots the figure of ORF encoding expression

%% Load data

clear all
B = readtable('~/Google Drive/CareyLab/ExternalData/Brauer08/TableS1.xls');
BASE04 = readtable( '~/Google Drive/CareyLab/ExternalData/Basehoar04/Basehoar04_TATA.csv' ,'FileType','text','Delimiter',';');
BASE04.ORF = BASE04.Gene  ; 
BASE04 = BASE04( : , {'x1_TATA_0_TATA_less' 'ORF'});
PP = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/protein_properties.tab','FileType','text','Delimiter','\t');
%CU = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/s_cerevisiae-codonusage.txt','FileType','text','Delimiter','\t');
T = innerjoin( B , PP ,'Key','ORF');
T = innerjoin( T , BASE04 ,'Key','ORF');
% T cointains all the information
%Load Brauer Raw data

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/Brauer08/');
DS1_Brauer = dataset('file','./DataSet1.tds','delimiter','\t');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/Brauer08/');
DS2_Brauer = dataset('file','./DataSet2.tds','delimiter','\t');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%ORF to Seq:
load('ORFtoSeq.mat');

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

% Load Keren13 Data about each promoter
cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/Keren13/');
KerenS4 = dataset('file','./Keren13_TableS4.txt','delimiter','\t');

cd('~/Google Drive/CareyLab/ExternalData/')
%Brauer 2008:
Brauer08 = dataset('XLSfile','./Brauer08/TableS1.xls');
%Athanasiadou 2017:
Athanasiadou17_Carbon = dataset('XLSfile','./Athanasiadou17/Athanasiadou_et_al_Suppl Table S2 Carbon.xlsx');
Athanasiadou17_Carbon.Properties.VarNames(1) = {'ORF'};
Athanasiadou17_Carbon.Properties.VarNames(8) = {'ExpConst_Carbon'};

Athanasiadou17_Nitrogen = dataset('XLSfile','./Athanasiadou17/Athanasiadou_et_al_Suppl Table S3 Nitrogen.xlsx');
Athanasiadou17_Nitrogen.Properties.VarNames(1) = {'ORF'};
Athanasiadou17_Nitrogen.Properties.VarNames(8) = {'ExpConst_Nitrogen'};

%vanDijk2015
vanDijk15 = dataset('file','./vanDijk15/FitFlowRNAseq.tab','delimiter','\t');
vanDijk15.Properties.VarNames(1) = {'Name'};


%% merge datasets

%A dataset that cointains the informatioN.
NativeDS = DS_cor;
NativeDS = join(NativeDS,ORFseqDS,'Keys','ORF','mergekeys',true);
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
NativeDS = join(NativeDS,UPF1del_Dataset,'Keys','ORF','mergekeys',true);

% Add missing features to the NativeDS:

%Correct Seq:
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq{I}; NativeDS.ORF_seq(I) = ORF;
end

%GC_content and length:
NativeDS.GC_content = 0;
NativeDS.ORF_length = 0;
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq(I); ORF = ORF{1}; ORF = ORF{1};
     (length(regexp(ORF,'C')) + length(regexp(ORF,'G')))./length(ORF);
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
NativeDS.log2expression = log2(NativeDS.UPF1mut_WT_Expression);

% Add Nag and Yas Expression

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

%Add Slope btw Growth and Expression

Tds = table2dataset(T);
NativeDS = join(Tds,NativeDS,'Keys','ORF','mergekeys',true,'type','outer');

% Add Keren 13
NativeDS = join(End5_Contr_adds,NativeDS,'Keys','ORF','type','outer','mergekeys',true);

%Add Athanasiadou17

NativeDS = join(Athanasiadou17_Nitrogen,NativeDS,'Keys','ORF','mergekeys',true,'type','outer');
NativeDS = join(Athanasiadou17_Carbon,NativeDS,'Keys','ORF','mergekeys',true,'type','outer');


%% Correlation between KerenVSgrowthRate and  Brauer SlopeVSGrowthRate

Conditions = KerenS4.Properties.VarNames(2:11)';

GR = log(2)./(double(KerenS4(3,2:11))./60);
MeanPRactivity = log2(double(KerenS4(1,2:11)));

clc
close all;
idx = NativeDS.Glucose > 0;

%idx = NativeDS.Glucose > 0 & NativeDS.BootstrappedP_valueOfSlope < 0.05 & NativeDS.GrowthRateSlope < 0;
All_DS = NativeDS(idx,:);
KerenExpression = NaN(length(All_DS.Glucose),10);
%generate matrix with Keren Expression and genes
I = 1;
for Condition = Conditions'
    if I == 4; Condition = {'Growth39C'}; end
    if I == 6; Condition = {'GluMinusAA'}; end
    if I == 7; Condition = {'GalMinusAA'}; end
            
    idxCond = find(strcmp(All_DS.Properties.VarNames,Condition));
    KerenExpression(:,I) = log2(double(All_DS(:,idxCond)));
    I = I + 1;  
end

%add to All_DS the slope between promoter and growth rate:
All_DS.GrowthRateVSpromoterSlope = 0;
All_DS.GrowthRateVSpromoterSlope_scaled = 0;

for I = 1:length(All_DS.ORF)
    %raw promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:),1);
    All_DS.GrowthRateVSpromoterSlope(I) = Parms(1);
    %scaled promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:)-(MeanPRactivity),1);
    All_DS.GrowthRateVSpromoterSlope_scaled(I) = Parms(1);
end
    fh = figure('units','centimeters','position',[5 5 6.6 6.6]); 
    hold on;

%plot things
Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
for TypeP = 1:2  
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end
    
    %Xval = All_DS.GrowthRateVSpromoterSlope(IDXtata); Xlabel = 'Slope PR-Promoter';
    
    Xval = All_DS.GrowthRateSlope(IDXtata); Xlabel = 'GR-Expresion slope (mRNA)';
    Yval = All_DS.GrowthRateVSpromoterSlope_scaled(IDXtata); Ylabel = 'GR-Expresion slope (Promoter)';
    
    X = Xval(Xval > -Inf & Yval > -Inf & Xval < Inf & Yval < Inf);
    Y = Yval(Xval > -Inf & Yval > -Inf & Xval < Inf & Yval < Inf);
    
    
    %plot(X,Y,'o','color',colors(TypeP,:),'markerfacecolor',colors(TypeP,:))
    s = scatter(X,Y,'o','markerfacecolor',colors(TypeP,:),'markeredgecolor',colors(TypeP,:),'markerfacealpha',0.2,'markeredgealpha',0.1)
    %dscatter(X,Y)
    uistack(s,'bottom')
    [r,p] = corr(X,Y)
    %plot(linspace(min(X),max(X),3),linspace(0,0,3),'-','color',[.7 .7 .7])
    %plot(linspace(0,0,3),linspace(min(Y),max(Y),3),'-','color',[.7 .7 .7])
    plot(linspace(-25,20,3),linspace(-25,20,3),'-','color',[.7 .7 .7])
    
    %title(strcat('r =',{' '},num2str(r,'%.3f'),{' '},'p =',{' '},num2str(p,'%.10f')))
    
    xlim([-22,20])
    ylim([-22,20])
    
    
    xlabel(Xlabel)
    ylabel(Ylabel)
    set(gca,'fontsize',11)

end

%% Histograms of Slope-VS expression

%idx = NativeDS.BootstrappedP_valueOfSlope < 0.05;
%idx = NativeDS.ExpConst_Carbon > 0;
idx = NativeDS.ExpConst_Carbon > -Inf;

All_DS = NativeDS(idx,:);

%All_DS.Slope_GR_mRNA = All_DS.GrowthRateSlope; TypeData = 'Brauer 2008';
%All_DS.Slope_GR_mRNA = All_DS.ExpConst_Carbon; TypeData = 'Athanasiadou 2017 Carbon';
All_DS.Slope_GR_mRNA = All_DS.ExpConst_Nitrogen; TypeData = 'Athanasiadou 2017 Nitrogen';

clc
%plot things
Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
fh = figure('units','centimeters','position',[5 5 6.5 6.5]); 
hold on; grid on;

for TypeP = 1:2  
    if TypeP == 1
        IDXtata = (All_DS.TATA_type==1);
        Label = {'TATA'};
    end
    if TypeP == 2
        IDXtata = (All_DS.TATA_type==0);
        Label = {'TATA-less'};
    end
    
    [F,X] = ksdensity(All_DS.Slope_GR_mRNA(IDXtata)); Xlabel = 'GR-mRNA slope';

    plot(X,F,'-','color',colors(TypeP,:),'linewidth',2)

    xlabel(Xlabel)
    ylabel('fraction of genes')
    set(gca,'fontsize',11)
    
    if TypeP == 1
        TATA_pos = length(All_DS.Slope_GR_mRNA(IDXtata & All_DS.Slope_GR_mRNA>0))
        TATA_neg = length(All_DS.Slope_GR_mRNA(IDXtata & All_DS.Slope_GR_mRNA<0))
        
    end
    if TypeP == 2
        TATAless_pos = length(All_DS.Slope_GR_mRNA(IDXtata & All_DS.Slope_GR_mRNA>0))
        TATAless_neg = length(All_DS.Slope_GR_mRNA(IDXtata & All_DS.Slope_GR_mRNA<0))
        
    end
    
    xlim([-20,20])
   % set(gca,'xtick',[-15,0,15])
    title(TypeData)
end

[h,p] = fishertest([[TATA_pos,TATA_neg];[TATAless_pos,TATAless_neg]])


%% See if ORF features (or tAI) predict Slope better in each promoter cathegory
%Lucas Code

Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];
%predictors = {'GC_content','ORF_length','tAI','AG1','AC2','A2'}; Ylabel = 'Predicted (ORF)';
%predictors = {'GC_content','ORF_length','tAI','AG1','AC2','A2','GrowthRateVSpromoterSlope_scaled'}; Ylabel = 'Predicted (ORF+Promoter)';
predictors = {'GrowthRateVSpromoterSlope_scaled'}; Ylabel = 'Predicted (Promoter)';

%predictors = {'tAI'}; Ylabel = 'Predicted (tAI)';
%predictors = {'UTR3length'}; Ylabel = 'Predicted (UTR3length)';

Conditions = KerenS4.Properties.VarNames(2:11)';

GR = 1./double(KerenS4(3,2:11));
MeanPRactivity = log2(double(KerenS4(1,2:11)));

clc
idx = NativeDS.Glucose > 0 ;
%idx = NativeDS.Glucose > 0 & NativeDS.BootstrappedP_valueOfSlope < 0.05 & NativeDS.GrowthRateSlope < 0;
All_DS = NativeDS(idx,:);
KerenExpression = NaN(length(All_DS.Glucose),10);
%generate matrix with Keren Expression and genes
I = 1;
for Condition = Conditions'
    if I == 4; Condition = {'Growth39C'}; end
    if I == 6; Condition = {'GluMinusAA'}; end
    if I == 7; Condition = {'GalMinusAA'}; end
            
    idxCond = find(strcmp(All_DS.Properties.VarNames,Condition));
    KerenExpression(:,I) = log2(double(All_DS(:,idxCond)));
    I = I + 1;  
end
%add to All_DS the slope between promoter and growth rate:
All_DS.GrowthRateVSpromoterSlope = 0;
All_DS.GrowthRateVSpromoterSlope_scaled = 0;
for I = 1:length(All_DS.ORF_left)
    %raw promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:),1);
    All_DS.GrowthRateVSpromoterSlope(I) = Parms(1);
    %scaled promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:)-(MeanPRactivity),1);
    All_DS.GrowthRateVSpromoterSlope_scaled(I) = Parms(1);
end

fh = figure('units','centimeters','position',[5 5 6.6 6.6]); 
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
    
    DS = All_DS(IDXtata,:);
    
    %define predictors
    GoodIDX = linspace(1,length(DS.ORF_left),length(DS.ORF_left))';
    PredictorsMat = NaN(length(DS.ORF_left),length(predictors));
    I = 1;
    for predictor = predictors
        idx = find(strcmp(DS.Properties.VarNames,predictor));
        PredictorsMat(:,I) = double(DS(:,idx));
        idx = double(DS(:,idx))>-Inf & double(DS(:,idx))<Inf;
        GoodIDX(~idx) = [];
        I = I + 1;
    end
    X = PredictorsMat(GoodIDX,:);
    Y = DS.GrowthRateSlope(GoodIDX,:);
        
    K = 10;    % # of cross-validations
    c = cvpartition( length(X(:,1)) , 'Kfold' , K);
    mdls = cell( K , 1);  % save each CV model
    Ypred = NaN(length(X(:,1)),1);  % cross-validated predicted values
    for I = 1:K
        mdls{I} = fitglm( X(c.training(I),:) ,  Y(c.training(I)) );
        Ypred( c.test(I) ) =  predict( mdls{I} , X(c.test(I),:));
    end
    
    GoodIDX = find(Y > -Inf & Ypred > -Inf);
    Y = Y(GoodIDX);
    Ypred = Ypred(GoodIDX);
    corr( Y , Ypred)
    rsq = rsquare( Y , Ypred) % this for figure
    corr( Y , Ypred)^2

   rsq_cv_mean = mean( cellfun( @(X)X.Rsquared.Ordinary , mdls) )
   rsq_cv_std  = std( cellfun( @(X)X.Rsquared.Ordinary , mdls) )

    X = Y;
    Y = Ypred;
    
    Gray = [.7 .7 .7];
    %plot(linspace(min([X;Y]),max([X;Y]),3),linspace(min([X;Y]),max([X;Y]),3),'-','color',[.7 .7 .7])

    %plot(X,Y,'o','markersize',4,'markerfacecolor',[0,0.4,0.8],'color','k')
    %plot(X,Y,'o','markersize',9,'markerfacecolor','k','color','k')
    %plot(X,Y,'o','markersize',6,'color','k')
    s = scatter(X,Y,'o','markerfacecolor',colors(TypeP,:),'markeredgecolor',colors(TypeP,:),'markerfacealpha',0.3,'markeredgealpha',0.1)
    plot(linspace(-25,20,3),linspace(-25,20,3),'-','color',[.7 .7 .7])
    
    %title(strcat('r =',{' '},num2str(r,'%.3f'),{' '},'p =',{' '},num2str(p,'%.10f')))
    
    xlim([-15,15])
    ylim([-15,15])
    
    uistack(s,'bottom')
    %plot(X,Y,'o','markersize',6,'markerfacecolor',Gray,'color','k','linewidth',0.7)
    xlabel('GR - mRNA slope');
    %ylabel('GR - (ORF + Promoter) slope');
    %ylabel('GR - ORF slope');
    ylabel('GR - Promoter slope');

    r2 = rsquare(X(X>-Inf & Y>-Inf),Y(X>-Inf & Y>-Inf));
    %text(mean(X),mean(Y),strcat('R^2 =',{' '},num2str(r2)),'fontsize',10)
    %set(gca,'fontsize',11,'xtick',[0,5,10],'xtick',[-10,-5],'ytick',[-10,-5])
    set(gca,'fontsize',11,'xtick',[0,5,10],'xtick',[-15,0,15],'ytick',[-15,0,15])

end

%% See if ORF features (or tAI) predict Slope better in each promoter cathegory. Bar Plot
%Lucas Code

Red = [0.8,0,0];
Blue = [39,101,186]./255;
colors = [Red;Blue];

%predictors = {'tAI'}; Ylabel = 'Predicted (tAI)';
%predictors = {'UTR3length'}; Ylabel = 'Predicted (UTR3length)';

Conditions = KerenS4.Properties.VarNames(2:11)';

GR = 1./double(KerenS4(3,2:11));
MeanPRactivity = log2(double(KerenS4(1,2:11)));

clc
idx = NativeDS.Glucose > 0 ;
%idx = NativeDS.Glucose > 0 & NativeDS.BootstrappedP_valueOfSlope < 0.05 & NativeDS.GrowthRateSlope < 0;
All_DS = NativeDS(idx,:);
KerenExpression = NaN(length(All_DS.Glucose),10);
%generate matrix with Keren Expression and genes
I = 1;
for Condition = Conditions'
    if I == 4; Condition = {'Growth39C'}; end
    if I == 6; Condition = {'GluMinusAA'}; end
    if I == 7; Condition = {'GalMinusAA'}; end
            
    idxCond = find(strcmp(All_DS.Properties.VarNames,Condition));
    KerenExpression(:,I) = log2(double(All_DS(:,idxCond)));
    I = I + 1;  
end
%add to All_DS the slope between promoter and growth rate:
All_DS.GrowthRateVSpromoterSlope = 0;
All_DS.GrowthRateVSpromoterSlope_scaled = 0;
for I = 1:length(All_DS.ORF)
    %raw promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:),1);
    All_DS.GrowthRateVSpromoterSlope(I) = Parms(1);
    %scaled promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:)-(MeanPRactivity),1);
    All_DS.GrowthRateVSpromoterSlope_scaled(I) = Parms(1);
end

BarMatMean = NaN(3,2);
STDmat = [];
STDmatMean = [];

%define GR-Expression DataSet
%All_DS.Slope_GR_mRNA = All_DS.GrowthRateSlope; TypeData = 'Brauer 2008';
%All_DS.Slope_GR_mRNA = All_DS.ExpConst_Carbon; TypeData = 'Athanasiadou 2017 Carbon';
%All_DS.Slope_GR_mRNA = All_DS.ExpConst_Nitrogen; TypeData = 'Athanasiadou 2017 Nitrogen';
All_DS.Slope_GR_mRNA = mean([All_DS.ExpConst_Nitrogen,All_DS.ExpConst_Carbon],2); TypeData = 'Athanasiadou 2017 All';


for P = 1:3
    if P == 1
        predictors = {'GrowthRateVSpromoterSlope_scaled'}; 
    end
    if P == 2
        predictors = {'GC_content','ORF_length','tAI'}; 
        %predictors = {'tAI'}; 
        
    end
    if P == 3
        %predictors = {'GC_content','ORF_length','tAI','AG1','AC2','A2','GrowthRateVSpromoterSlope_scaled'}; Ylabel = 'Predicted (ORF+Promoter)';  
        predictors = {'GC_content','ORF_length','tAI','GrowthRateVSpromoterSlope_scaled'}; Ylabel = 'Predicted (ORF+Promoter)';  
    end

    for TypeP = 1:2

        if TypeP == 1
            IDXtata = (All_DS.TATA_type==1);
            Label = {'TATA'};
        end
        if TypeP == 2
            IDXtata = (All_DS.TATA_type==0);
            Label = {'TATA-less'};
        end

        DS = All_DS(IDXtata,:);

        %define predictors
        GoodIDX = linspace(1,length(DS.ORF),length(DS.ORF))';
        PredictorsMat = NaN(length(DS.ORF),length(predictors));
        I = 1;
        for predictor = predictors
            idx = find(strcmp(DS.Properties.VarNames,predictor));
            PredictorsMat(:,I) = double(DS(:,idx));
            idx = double(DS(:,idx))>-Inf & double(DS(:,idx))<Inf;
            GoodIDX(~idx) = [];
            I = I + 1;
        end
        X = PredictorsMat(GoodIDX,:);
        Y = DS.Slope_GR_mRNA(GoodIDX,:);

        K = 10;    % # of cross-validations
        c = cvpartition( length(X(:,1)) , 'Kfold' , K);
        mdls = cell( K , 1);  % save each CV model
        Ypred = NaN(length(X(:,1)),1);  % cross-validated predicted values
        for I = 1:K
            mdls{I} = fitglm( X(c.training(I),:) ,  Y(c.training(I)) );
            Ypred( c.test(I) ) =  predict( mdls{I} , X(c.test(I),:));
        end

        GoodIDX = find(Y > -Inf & Ypred > -Inf);
        Y = Y(GoodIDX);
        Ypred = Ypred(GoodIDX);
        corr( Y , Ypred);
        rsq = rsquare( Y , Ypred); % this for figure
        corr( Y , Ypred)^2;

       rsq_cv_mean = mean( cellfun( @(X)X.Rsquared.Ordinary , mdls) );
       rsq_cv_std  = sem( cellfun( @(X)X.Rsquared.Ordinary , mdls) );
       BarMatMean(P,TypeP) = rsq_cv_mean;
       STDmat = [STDmat;rsq_cv_std];
       STDmatMean = [STDmatMean ; rsq_cv_mean];
       
    end
end

fh = figure('units','centimeters','position',[5 5 6.6 6.6]); 
hold on;

b = bar(BarMatMean);
b(1).FaceColor = Red;
b(2).FaceColor = Blue;
Xvals = [0.85;1.15;1.85;2.15;2.85;3.15];
STDmatMean;
errorbar(Xvals,STDmatMean,STDmat,'linestyle','none','linewidth',2,'color','k')
ylabel('r^2 to predict GR-mRNA')
set(gca,'xtick',[],'fontsize',11)
title(TypeData)

%% Tiny plots of GR vs Promoter

%predictors = {'tAI'}; Ylabel = 'Predicted (tAI)';
%predictors = {'UTR3length'}; Ylabel = 'Predicted (UTR3length)';

Conditions = KerenS4.Properties.VarNames(2:11)';

GR = 1./double(KerenS4(3,2:11));
MeanPRactivity = log2(double(KerenS4(1,2:11)));

clc
idx = NativeDS.Glucose > 0 ;
%idx = NativeDS.Glucose > 0 & NativeDS.BootstrappedP_valueOfSlope < 0.05 & NativeDS.GrowthRateSlope < 0;
All_DS = NativeDS(idx,:);
KerenExpression = NaN(length(All_DS.Glucose),10);
%generate matrix with Keren Expression and genes
I = 1;
for Condition = Conditions'
    if I == 4; Condition = {'Growth39C'}; end
    if I == 6; Condition = {'GluMinusAA'}; end
    if I == 7; Condition = {'GalMinusAA'}; end
            
    idxCond = find(strcmp(All_DS.Properties.VarNames,Condition));
    KerenExpression(:,I) = log2(double(All_DS(:,idxCond)));
    I = I + 1;  
end

%add to All_DS the slope between promoter and growth rate:
All_DS.GrowthRateVSpromoterSlope = 0;
All_DS.GrowthRateVSpromoterSlope_scaled = 0;
for I = 1:length(All_DS.ORF_left)
    %raw promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:),1);
    All_DS.GrowthRateVSpromoterSlope(I) = Parms(1);
    %scaled promoter expression:
    Parms = polyfit(GR,KerenExpression(I,:)-(MeanPRactivity),1);
    All_DS.GrowthRateVSpromoterSlope_scaled(I) = Parms(1);
end


idx = find(All_DS.GrowthRateVSpromoterSlope_scaled < 0);
idx = idx(8);
Promoter = KerenExpression(idx,:)-(MeanPRactivity);
Parms = polyfit(GR,KerenExpression(idx,:)-(MeanPRactivity),1)
Yf = polyval(Parms,GR);

fh = figure('units','centimeters','position',[5 5 3 3]); 
hold on;

plot(GR,Promoter,'.','color',[.7 .7 .7],'markersize',25)
plot(GR,Yf,'-','color',[.7 .7 .7],'linewidth',3)
set(gca,'xtick',[],'ytick',[])


%% Tiny plots of GR vs mRNA


Conditions = {'G0_05','G0_1','G0_15','G0_2','G0_25','G0_3'};
idx = 500;
DS_gr = DS1_Brauer(idx,:);
DS_exp = DS2_Brauer(idx,:);

%define GR and expression
GR = [];
Exp = [];
I = 1;
for Condition = Conditions
    idxGR = find(strcmp(DS_gr.Properties.VarNames,Condition));
    idxExp = find(strcmp(DS_exp.Properties.VarNames,Condition));
    
    GR = [GR ; double(DS_gr(1,idxGR))];
    Exp = [Exp ; double(DS_exp(1,idxExp))];    
    I = I + 1;  
end

Parms = polyfit(GR,Exp,1);
Yf = polyval(Parms,GR);
fh = figure('units','centimeters','position',[5 5 3 3]); 
hold on;

plot(GR,Exp,'.','color','K','markersize',25)
plot(GR,Yf,'-','color','K','linewidth',3)
set(gca,'xtick',[],'ytick',[])

