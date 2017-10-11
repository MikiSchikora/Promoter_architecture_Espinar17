%% Description
%This Script plots 3'utr length VS NMD efficiency

%% Load data

clear all
cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Data/');
%UPF1 mutant data --> cointains the information for the 3'utr length:
load('./External Data/S_cer_UPF1_mutant_RNAseq.mat');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%TATA Types:
load('Promoter_Features.mat')
%TFIID/SAGA Types:
load('Yeast_SAGAvsTFIID.mat')

%% Join datasets and add things:

All_DS = G;
%SAGA vs TFIID:
for I = 1:length(All_DS.ORF)
    ORF = All_DS.ORF{I};
    idx = find(strcmp(SAGAvsTFIID_DS.ORF,ORF)); 
    if isempty(idx) == 0
        if strcmp(SAGAvsTFIID_DS.Type(idx),'TFIID-dominated') == 1
            All_DS.TFIID(I) = 1;
        else
            All_DS.TFIID(I,:) = 0;            
        end
        if strcmp(SAGAvsTFIID_DS.Type(idx),'SAGA-dominated') == 1
            All_DS.SAGA(I) = 1;
        else
            All_DS.SAGA(I) = 0;
        end       
    else
        All_DS.TFIID(I) = 0;
        All_DS.SAGA(I) = 0;
    end
end
%TATA feats:
for I = 1:length(All_DS.ORF)
    ORF = All_DS.ORF(I);
    idx = find(strcmp(AllPromotersDS.ORF,ORF));
    if isempty(idx) == 0;
        All_DS.TATA_type(I) = AllPromotersDS.TATA_type(idx);
        All_DS.Subclassification(I) = AllPromotersDS.Subclassification(idx);
    else
        All_DS.TATA_type(I,:) = 10000;
        All_DS.Subclassification(I) = 10000;        
    end    
end

%% Plot 3'utr VS NMD efficiency

utr3threshold_down = 200;
utr3threshold_up = 400;

Red = [0.8,0,0];
Blue = [39,101,186]./255;

idx = find(All_DS.UTR3length >= utr3threshold_down & All_DS.UTR3length <= utr3threshold_up);
G = All_DS(idx,:);
win_size = 60; %for closest data --> sliding number of points
N = 1;
fh = figure('units','centimeters','position',[5 5 9 6]);
hold on;
markers = '^o';
for Class = 2 %for each of the types of data
    %Select Classification
    if Class == 1; idx_1 = find(G.SAGA == 1); idx_2 = find(G.TFIID == 1); Labels = {'SAGA','TFIID'}; w = 10;end
    if Class == 2; idx_1 = find(G.TATA_type == 1); idx_2 = find(G.TATA_type == 0); Labels = {'TATA_1','TATA_0'}; w = 10;end
    %Select Range of data to filter by log2(WT expression):
        Order_of_extreme_values = sort([max(log2(G.median_tpm(idx_1)));max(log2(G.median_tpm(idx_2)));min(log2(G.median_tpm(idx_1)));min(log2(G.median_tpm(idx_2)))]);
        Min = Order_of_extreme_values(2);
        Max = Order_of_extreme_values(3);
        idx_filt = find(log2(G.median_tpm) >= Min & log2(G.median_tpm) <= Max);
        Filt_G = G(idx_filt,:);
    %Select DS:
    DS = G; DataSet = 'All Data';
    %     DS = Filt_G; DataSet = 'Overlapping log_2(WT expression)';
    %     GoodIDX_1 = [];
    %     for I = 1:length(idx_1);
    %         idx = find(idx_filt == idx_1(I));
    %         if  isempty(idx) == 0
    %             GoodIDX_1 = [GoodIDX_1 ; idx];
    %         end
    %     end
    %     idx_1 = GoodIDX_1;
    %     
    %     GoodIDX_2 = [];    
    %     for I = 1:length(idx_2);
    %         idx = find(idx_filt == idx_2(I));
    %         if  isempty(idx) == 0
    %             GoodIDX_2 = [GoodIDX_2 ; idx];
    %         end
    %     end
    %     idx_2 = GoodIDX_2;
    %     
    % end
    colors = [Red;Blue];
    %Plot for each promoter type:
    for Ptype = 1:2
        if Ptype == 1; idx = idx_1; end
        if Ptype == 2; idx = idx_2; end
        X = (DS.UTR3length(idx)) ; Xlabel = '3''utr length';
        Y = DS.upf_o_wt(idx); Ylabel = 'log_2 Exp_U_P_F_1_k_o/Exp_W_T';
        [Xs,Ys,E] = RNAlib_Smooth_and_plot(X,Y,win_size);
        %plot(Xs,smooth(Ys),'-','color',colors(Ptype,:),'linewidth',2,'markersize',4,'markerfacecolor',colors(Ptype,:))
        %scatter(X,Y,'o','markeredgecolor',colors(Ptype,:),'markerfacecolor',colors(Ptype,:),'markerfacealpha',0.1,'markeredgealpha',0.1,'linewidth',0.1,'SizeData',10)
        s = shadedErrorBar(Xs,smooth(Ys),E,{'color',colors(Ptype,:),'linewidth',2});
        s.patch.FaceAlpha = 0.5;


        [r,p] = corr(X,Y);
        %text(mean(X),mean(Y),strcat('r =',{' '},num2str(r),{' '},'p =',{' '},num2str(p)),'color',colors(Ptype,:),'fontsize',11)
        if Ptype == 1; Y1 = Y; X1 = X; end
        if Ptype == 2; Y2 = Y; X2 = X; end
        
    end
    N = N + 1;
end
[h,p] = ttest2(Y1,Y2)

xlim([200,400])
ylim([-0.2,0.6])
%title(DataSet)
%legend(Labels);
xlabel(Xlabel)
ylabel(Ylabel)
set(gca,'fontsize',11,'xscale','log')

%% Statistics
figure; hold on;

r_tata = bootstrp(100000,@corr,X1,Y1);
r_notata = bootstrp(100000,@corr,X2,Y2);
[~,p] = ttest2(r_tata,r_notata)
h = histogram(r_tata);
set(h,'facecolor',Red)
h = histogram(r_notata);
set(h,'facecolor',Blue)

set(gca,'xscale','lin')
xlabel('Correlation')
ylabel('# bootstraps')
title(Xlabel)
grid on
axis tight;

%Calculate stats

Fraction_negative = mean( (r_tata-r_notata) > 0)

