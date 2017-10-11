%% Description
%This is for generating the FACS analysis of NMD in the different strains

%% Load data

clear all
cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Flow Cytometry/Promoters_NMD_11032017_3');
load('DataSet.mat');
load('wells.mat');

%% Define vars

idxInt = find(strcmp(DataSet.Promoter,'Nothing') == 0);
Int_DS = DataSet(idxInt,:);
Int_wells = wells(idxInt);

%% Fold changes upon putting stop for each promoter and strain
DSidx = Int_DS.NCellsAfterFiltering >= 300 & Int_DS.NCellsAfterFiltering <= 60000 ...
    & Int_DS.Innoculation_OD < 0.2;
DS = Int_DS(DSidx,:);
%Strains = unique(DS.Strain);
%Strains = {'WT','YRA1','DBP2','NAM7'}';
Strains = {'WT';'DBP2';'NAM7';'WT';'DBP2';'NAM7'};
Promoters = {'RPL4','RPL4','RPL4','GALL','GALL','GALL'};
S = 1;
Pos = 1;
ttestCount = 0;
Ymat = [];
Gmat = [];
fh = figure('units','centimeters','position',[5 5 12 8]); 
hold on;
Blue = [39,101,186]./255;

Red = [0.8,0,0];

colors = [Blue ; Blue ; Blue ; Red ; Red ;Red ];

for S = 1:length(Strains)
    Strain = Strains(S);
    Promoter = Promoters(S);
    P = 1;
        idxSTOP = find(strcmp(DS.Promoter,Promoter) & strcmp(DS.Strain,Strain) ...
            & DS.STOP == 1 & DS.mean_PE_TexasRed_A > 0 & DS.mode_PE_TexasRed_A > 0);
        idxNOSTOP = find(strcmp(DS.Promoter,Promoter) & strcmp(DS.Strain,Strain) ...
            & DS.STOP == 0 &  DS.mean_PE_TexasRed_A > 0 &  DS.mode_PE_TexasRed_A > 0);
        
        mCherrySTOP = log2(DS.mean_PE_TexasRed_A(idxSTOP)); Title = 'Using mean';
        mCherryNOSTOP = log2(DS.mean_PE_TexasRed_A(idxNOSTOP));
        
        %mCherrySTOP = DS.mode_logRFP(idxSTOP); Title = 'Using mode';
        %mCherryNOSTOP = DS.mode_logRFP(idxNOSTOP);

        Ratios = [];
        % Make all possible ratios
        for RFPstop = mCherrySTOP'
            for RFPnostop = mCherryNOSTOP'
                Rat = RFPnostop-RFPstop; %NMD effect
                Ratios = [Ratios ; Rat];
            end
        end
        Positions = linspace(Pos,Pos,length(Ratios))';
        Ymat = [Ymat;Ratios];
        Gmat = [Gmat;Positions];
        
        jitterAmount = 0.2;
        jitterValuesX = 2*(rand(size(Positions))-0.5)*jitterAmount;   % +/-jitterAmount max
        plot(Positions+jitterValuesX, Ratios,'.','color',colors(S,:),'markersize',10);                        
        
        if S==1 ; Ywt_rpl4 = Ratios; end
        if S==2 ; Ydbp2_rpl4 = Ratios; end
        if S==3 ; Yupf1_rpl4 = Ratios; end
        if S==4 ; Ywt_gall = Ratios; end
        if S==5 ; Ydbp2_gall = Ratios; end
        if S==6 ; Yupf1_gall = Ratios; end
        PreviousRats = Ratios;
        P = P + 1;
        Pos = Pos + 1;
end
p = plot(linspace(min(Gmat)-1,max(Gmat)+1,3),linspace(0,0,3),'-','color',[.7 .7 .7]);
uistack(p,'bottom')
p = plot(linspace(min(Gmat)-1,max(Gmat)+1,3),linspace(0.88,0.88,3),'-','color',[.7 .7 .7]);
uistack(p,'bottom')
p = plot(linspace(min(Gmat)-1,max(Gmat)+1,3),linspace(1.84,1.84,3),'-','color',[.7 .7 .7]);
uistack(p,'bottom')


%legend({'GALLpr','RPL4pr'})
h= boxplot(Ymat,Gmat,'notch','on','positions',unique(Gmat),'color','k');
set(h(7,:),'Visible','off')
%set(h,'linewidth',2);
set(gca,'xtick',[],'xticklabels',Strains,'fontsize',11)   
set(gca,'fontsize',11)
%ylabel('NMD effect (w/out_S_T_O_P / w/_S_T_O_P)')
ylabel('NMD')
%title(Title)
%%
%[h,p] = ttest2(Ydbp2_rpl4,Ywt_rpl4)
[h,p] = ttest2(Yupf1_gall,Yupf1_rpl4)

%% Boxplots of several things for each genotype and promoter

DSidx = Int_DS.NCellsAfterFiltering >= 300 & Int_DS.NCellsAfterFiltering <= 60000 ...
    & Int_DS.Innoculation_OD <= 0.2;
DS = Int_DS(DSidx,:);
%Strains = unique(DS.Strain);
Strains = {'DBP2','WT','NAM7'}';
fh = figure('units','centimeters','position',[5 5 20 10]); 
hold on;

S = 1;
Pos = 1;
ttestCount = 0;
Ymat = [];
Gmat = [];
colors = [[0.8,0,0];[0,0,1]];
for Strain = Strains'
    Promoters = unique(DS.Promoter);
    P = 1;
    for Promoter = Promoters'
        idxSTOP = find(strcmp(DS.Promoter,Promoter) & strcmp(DS.Strain,Strain) ...
            & DS.STOP == 1 & DS.mean_PE_TexasRed_A > 0 & DS.mode_PE_TexasRed_A > 0);
        idxNOSTOP = find(strcmp(DS.Promoter,Promoter) & strcmp(DS.Strain,Strain) ...
            & DS.STOP == 0 &  DS.mean_PE_TexasRed_A > 0 &  DS.mode_PE_TexasRed_A > 0);
        
        mCherrySTOP = log2(DS.mean_PE_TexasRed_A(idxSTOP)); Title = 'Using mean';
        mCherryNOSTOP = log2(DS.mean_PE_TexasRed_A(idxNOSTOP));
        GFPSTOP = log2(DS.mean_FITC_A(idxSTOP));
        GFPNOSTOP = log2(DS.mean_FITC_A(idxNOSTOP));
        Readthrough = GFPSTOP - mCherrySTOP;
                
%         mCherrySTOP = DS.mode_logRFP(idxSTOP); Title = 'Using mode';
%         mCherryNOSTOP = DS.mode_logRFP(idxNOSTOP);
%         GFPSTOP = DS.mode_logGFP(idxSTOP);
%         GFPNOSTOP = DS.mode_logGFP(idxNOSTOP);
%         Readthrough = GFPSTOP - mCherrySTOP;
        
        %Yvar = mCherrySTOP; Ylabel = 'log_2(mCherry_P_T_C)';
        %Yvar = mCherryNOSTOP; Ylabel = 'log_2(mCherry_n_o_P_T_C)';
        %Yvar = GFPSTOP; Ylabel = 'log_2(GFP_P_T_C)';
        %Yvar = GFPNOSTOP; Ylabel = 'log_2(GFP_n_o_P_T_C)';
        Yvar = Readthrough; Ylabel = 'log_2(GFP_s_t_o_p/mCherry_s_t_o_p)';
        
        CV = std(Yvar)/mean(Yvar);
        
        Positions = linspace(Pos,Pos,length(Yvar))';
        Ymat = [Ymat;Yvar];
        Gmat = [Gmat;Positions];
        
        jitterAmount = 0.1;
        jitterValuesX = 2*(rand(size(Positions))-0.5)*jitterAmount;   % +/-jitterAmount max
        plot(Positions+jitterValuesX, Yvar,'.','color',colors(P,:),'markersize',15);                        
        if P == 2 %perform the ttest based on every time the promoter is at 2
            [~,p] = ttest2(Yvar,PreviousRats,'tail','both');
            %text(Pos-0.5,0.3,num2str(p),'fontsize',13)
            PreviousRats = [];
        end
        Pos
        %text(Pos,mean(Yvar),num2str(CV));
        PreviousRats = Yvar;
        P = P + 1;
        Pos = Pos + 1;
    end
    S = S + 1;
end
%legend({'GALLpr','RPL4pr'})
h= boxplot(Ymat,Gmat,'notch','on','positions',unique(Gmat),'color','k');
set(h(7,:),'Visible','off')
set(h,'linewidth',2);
set(gca,'xtick',[1.5,3.5,5.5,7.5],'xticklabels',Strains,'fontweight','bold','fontsize',14)
%ylim([10,13.5]) %mCherry
%ylim([6.9,13.7]) %GFP
ylabel(Ylabel)
%title(Title)

