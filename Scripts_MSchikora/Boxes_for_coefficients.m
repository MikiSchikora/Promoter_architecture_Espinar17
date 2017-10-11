%% Description
%This generates a box with the coefficents for the gDNA library and native
%genes

%It copies a lot from Julia's code

%% load data and define things

cd ~/Google' Drive'/CareyLab/random' RNA library'/Article/

DataDir = '~/Google Drive/CareyLab/random RNA library/Article/Data/';

num_threads = 4;
Blue = [39,101,186]./255;
Red = [0.8,0,0];
ORFs = struct2dataset(fastaread([DataDir 'orf_coding_all.fasta']));
name = {};
for x = 1:length(ORFs)
    a = strsplit(ORFs.Header{x}, ' ');
    name = [name; a{1}];
end
ORFs.Name = name;
load ~/Google' Drive'/CareyLab/random' RNA library'/ExternalData/mRNA_seq_data_sets/mRNA_seq_external_datasets.mat
DS = join(mRNA_seq_external_datasets, ORFs, 'LeftKey', 'ORF_GENE_Name', 'RightKey', 'Name', 'type', 'left','RightVars', 'Sequence');
% % Inner joined DS
% DS = dataset('file', [DataDir 'expr_for_julia_innerjoin.tab'], 'ReadVarNames', true)
Exp_DS = dataset('file', [DataDir 'expr_for_julia_outterjoin.tab'], 'ReadVarNames', true);
%Rs = dataset(NaN(10,1), )
DS = join(DS, Exp_DS(:,{'x1_TATA_0_TATA_less' ,'ORF_GENE_Name'}), 'LeftKey', 'ORF_GENE_Name', 'RightKey', 'ORF_GENE_Name', 'type', 'left', 'RightVars', 'x1_TATA_0_TATA_less');

DS_nat = DS;

%% define DS_lib as the gDNA libraries

cd('~/Google Drive/CareyLab/random RNA library/Miki/Data/');
%All integrated gDNA random RNA library datasets, cointains many datasets:
load('Integrated_DataSets_2014_Adds.mat');
DS_lib = WT_RPL4;
DS_lib.Sequence = DS_lib.orf_seq;


%% Add features 
for D = 1:2
    if D==1; DS = DS_lib; end 
    if D==2; DS = DS_nat; end 
    
    DS.orf_l = cellfun(@length,DS.Sequence);
    DS = DS(~DS.orf_l == 0,:); %Remove seqs without coding sequence or
    DS = DS(~rem(DS.orf_l,3),:); %Remove seqs with sequence length not multiple of 3 (only one gene: Putative aryl-alcohol dehydrogenase AAD6 / YFL056C)
    DS.orf_GC = cellfun(@GCcontent , DS.Sequence);

    % Possible premature stops
    DS.orf_stop_codon_loc = arrayfun(@(x)cell2mat(regexp(x,'*','once')).*3 - 2, nt2aa(DS.Sequence));
    DS.premature_stop = DS.orf_stop_codon_loc./(DS.orf_l-2) < 1;
    DS.orf_pct_stop_loc = ((DS.orf_stop_codon_loc+2)*100)./DS.orf_l;
    DS.nonsense_stops = arrayfun(@cell2mat , regexp(nt2aa(DS.Sequence), '*'), 'UniformOutput', false) ;
    DS.nonsense_stops_number = cell2mat(arrayfun(@(x) size(DS.nonsense_stops{x}, 2)-1, 1:length(DS),'UniformOutput', false)');

    % tAI and CAI
    DS.tAI = cellfun(@calc_tAI , DS.Sequence );

    % base, dimers, codon, aa
    FDS = gDNA_BaseCombinationCount(DS, 'Sequence', {@basecount @dimercount @codoncount @aacount}, 'var_prefix', 'orf_');
    FDS = SeqHexamerFreq(DS.Sequence, FDS, 'hexamers', {'ACGTTA','GAAAGA'}, 'var_prefix', 'orf_');
    FDS = SeqCodonPositionNucFreq(DS.Sequence, FDS, 'var_prefix', 'orf_');


    DS.A2 = FDS.orf_A2;
    DS.AC2 = FDS.orf_AC2;
    DS.AG1 = FDS.orf_AG1;
    DS.GC3 = FDS.orf_CG3;
    
    DS.ACGTTA = FDS.orf_ACGTTA;
    DS.GAAAGA = FDS.orf_GAAAGA;

    if D==1; DS_lib = DS; end 
    if D==2; DS_nat = DS; end 

end

DS_nat.utr3_l = 0;
DS_nat.Expression = DS_nat.Nagalaskhmi;

%% Prepare predictors and response variable

predictors = {'A2','AC2','AG1','tAI','orf_l','orf_gc','utr3_l','premature_stop','GC3'};
ps = 0.5; % pseudocounts


DS = DS_nat;
rnaseq = log10(DS.Expression + ps);
DS.RNAseq = rnaseq;
Library = 'Native Genes';

% idx = DS_lib.insert_l > 100 & DS_lib.insert_l <= 400;
% DS = DS_lib(idx,:);
% rnaseq = DS.RPL4;
% DS.RNAseq = rnaseq;
% Library = 'RPL4A lib';
% 
% idx = DS_lib.insert_l > 100 & DS_lib.insert_l <= 400;
% DS = DS_lib(idx,:);
% rnaseq = DS.Expression;
% DS.RNAseq = rnaseq;
% Library = 'GALL lib';
% 

X_direcion = [DS.A2, DS.AC2, DS.AG1, DS.GC3];
X_tai = [DS.tAI];
X_insert = [DS.orf_l, DS.orf_GC, DS.utr3_l, DS.premature_stop ];
X = zscore(double([X_direcion, X_tai, X_insert]));
Y = rnaseq;
X = X(~isnan(Y),:);
Y = Y(~isnan(Y),:);
N = length(Y);
% CV Model

K = 10;

idx_cv = crossvalind('Kfold', length(X), K);
S_cv = struct('Rsquared',NaN(K,1), 'CoefficientsEstimate', mat2dataset(NaN(K, length(predictors)+1 )), 'Pvalue',  mat2dataset(NaN(K, length(predictors)+1 )));
cv_Ypred = [];
cv_Ytest = [];

for I=1:K
    idx_test = idx_cv == I; idx_train = ~idx_test;
    R_glm = GeneralizedLinearModel.fit(X(idx_train,:), Y(idx_train), 'linear', 'VarNames', [predictors 'Expression']);
    S_cv.CoefficientsEstimate(I,:) = mat2dataset(R_glm.Coefficients.Estimate');
    S_cv.Pvalue(I,:) = mat2dataset(R_glm.Coefficients.pValue');
    Ypred = predict(R_glm, X(idx_test,:));
    r =  nancorr(Ypred, Y(idx_test));
    S_cv.Rsquared(I) = r^2;
    
    cv_Ypred = [cv_Ypred; Ypred];
    cv_Ytest = [cv_Ytest; Y(idx_test)];
    
end

fh = figure('units','centimeters','position',[5 5 6.5 6.5]);
hold on;
[c,p] = corr( cv_Ytest , cv_Ypred ,'rows','complete');
%plot(cv_Ytest, cv_Ypred,'ok','MarkerFaceColor',[1 0 0], 'MarkerEdgeColor', [.2 .2 .2], 'DisplayName', sprintf('Correlation=%0.03f',c) );
dscatter(cv_Ytest, cv_Ypred);
text(min(cv_Ytest)+0.1, max(cv_Ypred)-0.05, ['R^{2} = ' num2str(c^2, '%.3f')],'FontSize', 12);
mean(S_cv.Rsquared)
ylabel('Predicted (ORF model)')
xlabel('mRNA level (log_{10}(FPKM))')
title(Library)
set(gca,'fontsize',11)

%% Generate the boxes

predictors = {'A2','AC2','AG1','tAI','orf_l','orf_gc','utr3_l'};
ps = 0.5; % pseudocounts
clc
fh = figure('units','centimeters','position',[5 5 6.5 6.5]);
hold on;

for D = [1,4,5]
    if D == 1   
        DS = DS_nat;
        rnaseq = log10(DS.Expression + ps);
        DS.RNAseq = rnaseq;
        Library = 'Native Genes';
        color = [.7 .7 .7];
    end
    if D == 2   
        idx = DS_nat.x1_TATA_0_TATA_less == 1;
        DS = DS_nat(idx,:);
        rnaseq = log10(DS.Expression + ps);
        DS.RNAseq = rnaseq;
        Library = 'Native Genes';
        color = 'm';
    end
    if D == 3  
        idx = DS_nat.x1_TATA_0_TATA_less == 0;
        DS = DS_nat(idx,:);
        rnaseq = log10(DS.Expression + ps);
        DS.RNAseq = rnaseq;
        Library = 'Native Genes';
        color = 'c';
    end    
    if D == 4       
        idx = DS_lib.insert_l > 100 & DS_lib.insert_l <= 400;
        DS = DS_lib(idx,:);
        rnaseq = DS.RPL4;
        DS.RNAseq = rnaseq;
        Library = 'RPL4A lib';
        color = Blue;
    end
    if D == 5
        idx = DS_lib.insert_l > 100 & DS_lib.insert_l <= 400;
        DS = DS_lib(idx,:);
        rnaseq = DS.Expression;
        DS.RNAseq = rnaseq;
        Library = 'GALL lib';
        color = Red;
    end
    X_direcion = [DS.A2, DS.AC2, DS.AG1];
    X_tai = [DS.tAI];
    X_insert = [DS.orf_l, DS.orf_GC, DS.utr3_l];
    X = zscore(double([X_direcion, X_tai, X_insert]));
    Y = rnaseq;
    X = X(~isnan(Y),:);
    Y = Y(~isnan(Y),:);
    N = length(Y);
    % CV Model
    K = 10;
    idx_cv = crossvalind('Kfold', length(X), K);
    S_cv = struct('Rsquared',NaN(K,1), 'CoefficientsEstimate', mat2dataset(NaN(K, length(predictors)+1 )), 'Pvalue',  mat2dataset(NaN(K, length(predictors)+1 )));
    cv_Ypred = [];
    cv_Ytest = [];
    for I=1:K
        idx_test = idx_cv == I; idx_train = ~idx_test;
        R_glm = GeneralizedLinearModel.fit(X(idx_train,:), Y(idx_train), 'linear', 'VarNames', [predictors 'Expression']);
        S_cv.CoefficientsEstimate(I,:) = mat2dataset(R_glm.Coefficients.Estimate');
        S_cv.Pvalue(I,:) = mat2dataset(R_glm.Coefficients.pValue');
        Ypred = predict(R_glm, X(idx_test,:));
        r =  nancorr(Ypred, Y(idx_test));
        S_cv.Rsquared(I) = r^2;
        cv_Ypred = [cv_Ypred; Ypred];
        cv_Ytest = [cv_Ytest; Y(idx_test)];
    end
    
    
    %Generate Boxes
    boxplot(1,1)
    Coefficients = mean(double(S_cv.CoefficientsEstimate));
    Coefficients = Coefficients(2:end);
    for I = 1:length(Coefficients)
        Coefficient = Coefficients(I);
        plot(Coefficient,I,'o','color',color,'markerfacecolor',color)
        plot(linspace(-2,2,3),linspace(I,I,3),'-','linestyle',':','color',[.7 .7 .7])
    end
    D
    [c,p] = corr( cv_Ytest , cv_Ypred ,'rows','complete');
    r2 = c^2

    
    
end
predictorsNames = {'A2','AC2','AG1','tAI','length','GC','3''utr l'};
plot(linspace(0,0,7),linspace(0,8,7),'-','color',[.7 .7 .7])

set(gca,'ytick',linspace(1,7,7),'yticklabels',predictorsNames,'xtick',[-0.5,0,0.2],'xticklabels',[-0.5,0,0.2],'fontsize',11)
xlim([-2,2])
ylim([0,8])


