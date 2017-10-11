%% Description
%Here  we investigate how CAI and NMD affect expression in the different libraries

%% Load Data

clear all
cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
load('Integrated_DataSets.mat');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%All integrated gDNA random RNA library datasets, cointains many datasets:
load('Integrated_DataSets_2014_Adds.mat');
load('ORFtoSeq.mat');


%% Prepare data
DSL = All_DS; Year = '2015'; %Integrated sequences
%DSL = WT_2014; Year = '2014';
%% Add features

% Correct Expression by Several Parameters, WT
predictors = {'insert_l','orf_GC'};
mdl = LinearModel.fit(DSL, 'ResponseVar', 'Expression', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.Expr_corrected = Expr_corrected;

% Correct Expression by Several Parameters, UPF1
mdl = LinearModel.fit(DSL, 'ResponseVar', 'UPF1', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.UPF1_corrected = Expr_corrected;

% Correct Expression by Several Parameters, RPL4
mdl = LinearModel.fit(DSL, 'ResponseVar', 'RPL4', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.RPL4_corrected = Expr_corrected;

% Correct Expression by Several Parameters, DHH1
predictors = {'insert_l','orf_GC'};
mdl = LinearModel.fit(DSL, 'ResponseVar', 'DHH1', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.DHH1_corrected = Expr_corrected;

% Correct Expression by Several Parameters, XRN1
predictors = {'insert_l','orf_GC'};
mdl = LinearModel.fit(DSL, 'ResponseVar', 'XRN1', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.XRN1_corrected = Expr_corrected;

% Correct Expression by Several Parameters, WT_2014
predictors = {'insert_l','orf_GC'};
mdl = LinearModel.fit(DSL, 'ResponseVar', 'WT2014', 'PredictorVars', predictors );
Expr_corrected = mdl.Residuals.Raw;
DSL.WT2014_corrected = Expr_corrected;

%NSS addings
DSL = gDNALib_AddSeqFeatures__DownstreamAUGs(DSL);
DSL.NScorLength = DSL.nonsense_stops_number ./ DSL.insert_l;
DSL.nonsense_stops_number_class = DSL.nonsense_stops_number;
more_than_10_idx = find(DSL.nonsense_stops_number > 10);
DSL.nonsense_stops_number_class(more_than_10_idx) = 11;

% Add distances between stops (in terms of number of codons in between)
for I = 1:length(DSL.nonsense_stops)
MatNSpos_Prep = DSL.nonsense_stops(I); %This way it has to fit
MatNSpos = MatNSpos_Prep{1,1};
DistMat = [];
for J = 1:length(MatNSpos) - 1
DistJ =  MatNSpos(J + 1) - MatNSpos(J) - 1;
DistMat = [DistMat , DistJ];
end
DSL.Distance_btw_NSs{I} = DistMat;
DSL.mean_Distance_btw_NSs(I) = mean(DistMat);
DSL.median_Distance_btw_NSs(I) = median(DistMat);
DSL.sum_Distance_btw_NSs(I) = sum(DistMat);
DSL.max_Distance_btw_NSs(I) = max(DistMat);
DSL.min_Distance_btw_NSs(I) = min(DistMat);
DSL.std_Distance_btw_NSs(I) = std(DistMat);
DSL.sem_Distance_btw_NSs(I) = sem(DistMat);
end  

% Types of Stop
for I = 1:length(DSL.orf_seq)
ORF = DSL.orf_seq(I); 
ORF = ORF{1,1};
DSL.StopCodon(I) = cellstr(ORF(length(ORF)-2:length(ORF)));
if strcmp('TGA',DSL.StopCodon(I)) == 1;DSL.StopID(I) = 1;end;
if strcmp('TAA',DSL.StopCodon(I)) == 1;DSL.StopID(I) = 2;end;
if strcmp('TAG',DSL.StopCodon(I)) == 1;DSL.StopID(I) = 3;end;
end

% Add Transcript Sequence from ATG on
for I = 1:length(DSL.trans_seq)
trans_seq = DSL.trans_seq{I};    
DSL.no3utr_seq(I) = cellstr(trans_seq(112:length(trans_seq)));        
end

% Add all numbers and IDs of Stop Codons
for I = 1:length(DSL.no3utr_seq)
Seq = DSL.no3utr_seq{I}; 
NSpos = (DSL.nonsense_stops{I}.*3) - 2;
StopIDmat = [];
for J = NSpos
Stop = cellstr(Seq(J:J+2));
if strcmp('TGA',Stop) == 1;StopID = 1;end;
if strcmp('TAA',Stop) == 1;StopID = 2;end;
if strcmp('TAG',Stop) == 1;StopID = 3;end;
StopIDmat = [StopIDmat , StopID];
end
DSL.StopIDs{I} = StopIDmat;
DSL.NSS_diversity(I) = std(StopIDmat);
DSL.NSS_prevalent(I) = mode(StopIDmat);
end

% Add ORF coverage and NSS position
DSL.pctg_ORF = (DSL.orf_stop_codon_loc ./ (DSL.insert_l + 5)) .* 100; % +5 refers to the initial nucleotides of the ORF
DSL.pctg_ORF(DSL.pctg_ORF > 100) = 100;

%% Make many plots in which we see 3'utr and CAI vs expression, in the different libraries (current 1F,G and 2C,D)

%idx_len = DSL.insert_l <= 400 & DSL.insert_l >= 100;
idx_len = DSL.insert_l <= 10000 & DSL.insert_l >= 100 & DSL.premature_stop == 0; %for the NMD plots

DS = DSL(idx_len,:);
clc
%X = (DS.utr3_l); Xlabel = '3''utr length'; win_size = rem;
%X = DS.CAI; Xlabel = 'Codon Adaptation Index (CAI)';
X = DS.tAI; Xlabel = 'tRNA Adaptation Index (tAI)'; win_size = 380;

Red = [0.8,0,0];
Blue = [39,101,186]./255;
%Ygal = DS.Expression; Ylabel = 'Expression';
%Yrpl = DS.RPL4; 
%Ymut = DS.UPF1;
%Ymut = DS.DHH1;

%Ygal = (DS.Expr_corrected); Ylabel = 'Expression';
%Yrpl = (DS.RPL4_corrected); 
%Ymut = (DS.UPF1_corrected);
%Ymut = (DS.DHH1_corrected);
Ygal = zscore(DS.Expr_corrected); Ylabel = 'Expression';
Yrpl = zscore(DS.RPL4_corrected); 
Ymut = zscore(DS.UPF1_corrected);
%Ymut = zscore(DS.DHH1_corrected);
Yxrn1 = zscore(DS.XRN1_corrected);
%Y2014 = zscore(DS.WT2014_corrected);


fh = figure('units','centimeters','position',[5 5 9 6]); 
hold on;

%MUT:
%[Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Ymut,win_size);
%plot(Xs,smooth(Ys),'-o','color',[.7 .7 .7],'markerfacecolor',[.7 .7 .7],'linewidth',3,'markerfacecolor',[.7 .7 .7],'markersize',7)
%plot(Xs,smooth(Ys),'o','color','k','markersize',9,'linewidth',2)

%GALL
%[Xs,Ys,E,All_X_gal,All_Y_gal] = RNAlib_Smooth_and_plot(X,Ygal,win_size);
[Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Ygal,win_size);

plot(Xs,smooth(Ys),'-','color',Red,'markerfacecolor',Red,'linewidth',2,'markerfacecolor',Red,'markersize',4)
scatter(X,Ygal,'o','markeredgecolor',Red,'markerfacecolor',Red,'markerfacealpha',0.1,'markeredgealpha',0.1,'linewidth',0.1,'SizeData',10)
%s = shadedErrorBar(Xs,smooth(Ys),E,{'color',Red,'linewidth',2});
%s.patch.FaceAlpha = 0.5;

%RPL4:
%[Xs,Ys,E,All_X_rpl,All_Y_rpl] = RNAlib_Smooth_and_plot(X,Yrpl,win_size);
[Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Yrpl,win_size);

plot(Xs,smooth(Ys),'-','color',Blue,'markerfacecolor',Blue,'linewidth',2,'markerfacecolor',Blue,'markersize',4)
scatter(X,Yrpl,'o','markeredgecolor',Blue,'markerfacecolor',Blue,'markerfacealpha',0.1,'markeredgealpha',0.1,'linewidth',0.1,'SizeData',10)
%s = shadedErrorBar(Xs,smooth(Ys),E,{'color',Blue,'linewidth',2});
%s.patch.FaceAlpha = 0.5;
% add asterisks:
% for Xt = unique(All_X_gal)'
%     idx = find(All_X_gal == Xt);
%     [~,p] = ttest2(All_Y_gal(idx),All_Y_rpl(idx));
%     if p <= 0.05
%         %text(Xt,mean([median(All_Y_gal(idx)),median(All_Y_rpl(idx))]),'*','fontsize',25)
%     end
%     
% end

%plot(Xs,smooth(Ys),'o','color','k','markersize',9,'linewidth',2)

%xrn1:
%[Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Yxrn1,win_size);
%plot(Xs,smooth(Ys),'-o','color','k','markerfacecolor','k','linewidth',3,'markerfacecolor','k','markersize',7)
%plot(Xs,smooth(Ys),'o','color','k','markersize',9,'linewidth',2)

%WT2014:
%[Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Y2014,win_size);
%plot(Xs,smooth(Ys),'-o','color','m','markerfacecolor','m','linewidth',3,'markerfacecolor','m','markersize',7)
%plot(Xs,smooth(Ys),'o','color','k','markersize',9,'linewidth',2)


 [r,p] = corr(X,Ygal);
 %text(mean(X),mean(Ygal),strcat('r =',{' '},num2str(r),{' '},'p =',{' '},num2str(p)),'color',Red,'fontsize',11)
 [r,p] = corr(X,Yrpl);
 %text(mean(X),mean(Yrpl),strcat('r =',{' '},num2str(r),{' '},'p =',{' '},num2str(p)),'color',Blue,'fontsize',11)
% [r,p] = corr(X,Ymut)
% text(mean(X),mean(Ymut),strcat('r =',{' '},num2str(r)),'color',[.7 .7 .7],'fontsize',18,'fontweight','bold')

%legend({'GAL1p','RPL4Ap'})
%xlim([190,400])
ylim([-0.7,1])
xlabel(Xlabel)
ylabel(Ylabel)
%title(Year)
set(gca,'fontsize',11)

%% Check if the slopes are different

figure; hold on;
colors = [Red;Blue];
for Run = 1:2
    if Run == 1; Y = Ygal; end
    if Run == 2; Y = Yrpl; end   
    for N = round(logspace(1,log10(length(Ygal)),100))
        idx = randsample(length(Ygal),N);
        mdl = GeneralizedLinearModel.fit( X(idx) , Y(idx) )    ;
        ci = coefCI(mdl) ;
        plot( [N N] , ci(2,:)  ,'-o','color',colors(Run,:) );
    end
end

set(gca,'xscale','log')
ylabel('Confidence interval of slope')
xlabel('# data points')
grid on
axis tight;


%Check if the bootstrap distributions overlap

figure; hold on;

r_gal = bootstrp(10000,@corr,X,Ygal);
r_rpl = bootstrp(10000,@corr,X,Yrpl);
[~,p] = ttest2(r_gal,r_rpl)
h = histogram(r_gal);
set(h,'facecolor',Red)
h = histogram(r_rpl);
set(h,'facecolor',Blue)

set(gca,'xscale','lin')
xlabel('Correlation')
ylabel('# bootstraps')
title(Xlabel)
grid on
axis tight;

%Calculate stats

Fraction_negative = mean( (r_gal-r_rpl) > 0)


%% Plot histograms


%idx_len = DSL.insert_l <= 400 & DSL.insert_l >= 100;
idx_len = DSL.insert_l <= 400 & DSL.insert_l >= 100 & DSL.premature_stop == 0; %for the NMD plots

DS = DSL(idx_len,:);


Xgal = (DS.Expr_corrected); 
Xrpl = (DS.RPL4_corrected); 
Xupf1 = (DS.UPF1_corrected);
Xdhh1 = (DS.DHH1_corrected);
Xxrn1 = (DS.XRN1_corrected);

%[C,X] = ksdensity(Xgal); Title = 'GALL_p_r';
[C,X] = ksdensity(Xdhh1); Title = 'DHH1';

fh = figure('units','centimeters','position',[5 5 18 10]); 
hold on;
plot(X,C,'-','color','k','linewidth',2)
xlabel('Zscore Expression')
ylabel('Fraction seqs')
title(Title)
set(gca,'fontweight','bold','fontsize',14)

