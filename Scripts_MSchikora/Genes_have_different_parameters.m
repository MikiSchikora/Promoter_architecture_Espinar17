%% Description
%This generates the distributions of many parameters for the TATA and
%TFIID/SAGA genes

%% Load data

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%curated DEE datasets:
load('All_RNA_seq_S_cer_integrated_All_DS.mat')

%% Plot histograms

idx = DS_cor.Promoter_YFP > 0 & DS_cor.tAI < 0.5;
%idx = DS_cor.UPF1mut_WT_Expression > 0 & DS_cor.CAI < 0.8;

All_DS = DS_cor(idx,:);
clc
Red = [0.8,0,0];
fh = figure('units','centimeters','position',[5 5 20 10]);
hold on;

for Class = 2
   
   %X = log10(All_DS.thalf); Xlabel = 'log_1_0 RNA stability (Neymotin 2014)';
   %X = log10(All_DS.Keren_Expression); Xlabel = 'Expression';
   %X = (All_DS.CAI); Xlabel = 'CAI';
   X = (All_DS.tAI); Xlabel = 'tAI';
   
   
   GoodIDX = find(X>0);
   if Class == 1
       idx1 = find(All_DS.TFIID(GoodIDX) == 1); 
       idx2 = find(All_DS.SAGA(GoodIDX) == 1);
       Title = {'TFIID vs SAGA'};       
   end  
   if Class == 2
       idx1 = find(All_DS.TATA_type(GoodIDX) == 0); 
       idx2 = find(All_DS.TATA_type(GoodIDX) == 1);
       Title = {'TATA-less vs TATA'};       
   end      
   
   X = X(GoodIDX);

   %plot each type
   [C1,X1] = ksdensity(X(idx1));
   plot(X1,C1,'-','color','b','linewidth',2)
   [C2,X2] = ksdensity(X(idx2));
   plot(X2,C2,'-','color',Red,'linewidth',2)
   %plot line:
   %plot(linspace(0.5,0.5,3),linspace(min(C1),max(C1),3),'color',[.7 .7 .7],'linewidth',1.5)
   xlabel(Xlabel)
   set(gca,'xscale','lin','fontweigh','bold','fontsize',14)
   title(Title)
end

[h,p] = ttest2(X1,X2)
%% Plot correlation plots

%idx = DS_cor.Promoter_YFP > 0;
idx = DS_cor.UPF1mut_WT_Expression > 0;

All_DS = DS_cor(idx,:);

win_size = 40;
Red = [0.8,0,0];
fh = figure('units','centimeters','position',[5 5 20 10]);
for Class = 1:2
   
   X = (All_DS.CAI); Xlabel = 'CAI';
   Y = (All_DS.tAI); Ylabel = 'tAI';
   
   
   GoodIDX = find(X>0 & Y>0);
   if Class == 1
       idx1 = find(All_DS.TFIID(GoodIDX) == 1); 
       idx2 = find(All_DS.SAGA(GoodIDX) == 1);
       Title = {'TFIID vs SAGA'};       
   end  
   if Class == 2
       idx1 = find(All_DS.TATA_type(GoodIDX) == 0); 
       idx2 = find(All_DS.TATA_type(GoodIDX) == 1);
       Title = {'TATA-less vs TATA'};       
   end      
   
   X = X(GoodIDX);
   Y = Y(GoodIDX);
   subplot(1,2,Class)
   hold on;
   %plot each type
   
   %plot(X(idx1),Y(idx1),'.','color','b','linewidth',2)
   %plot(X(idx2),Y(idx2),'.','color',Red,'linewidth',2)
   %plot lines:
   plot(linspace(0.8,0.8,3),linspace(min(Y),max(Y),3),'color',[.7 .7 .7],'linewidth',1.5)
   %plot grouped:
   [Xs,Ys,E] = RNAlib_Smooth_and_plot(X(idx1),Y(idx1),win_size);
   %plot(Xs,Ys,'o','color','b','markerfacecolor','b')
   errorbar(Xs,Ys,E,'o','color','b','markerfacecolor','b','linewidth',2)
   [Xs,Ys,E] = RNAlib_Smooth_and_plot(X(idx2),Y(idx2),win_size);
   %plot(Xs,Ys,'o','color',Red,'markerfacecolor',Red)
      errorbar(Xs,Ys,E,'o','color',Red,'markerfacecolor',Red,'linewidth',2)

   xlabel(Xlabel)
   ylabel(Ylabel)
   set(gca,'xscale','lin','fontweigh','bold','fontsize',14)
   title(Title)
end

