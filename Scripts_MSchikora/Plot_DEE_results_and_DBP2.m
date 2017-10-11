%% Description 
%This plots the result of the DEE and detailed for DBP2

%% Load Data

clear all;
cd('~/Google Drive/CareyLab/random RNA library/Miki/Data/');
%curated DEE datasets:
load('All_RNA_seq_S_cer_integrated_All_DS.mat')

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Data/');
%DEE_DS after analysis:
load ('DEE_DS.mat') 

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/DigitalExpressionExplorer/')
%SRA_S dataset:
SRAs_DS = dataset('file','./metadata.tab','Delimiter','\t');  


%% Plot DEE analysis result (Fig 4B)

idx = find(DEE_DS.r_with_YFP > 0);
DS = DEE_DS(idx,:);

%X = DS.r_with_YFP.^2; Xlabel = 'r^2 YFP vs Expression';
%X = DS.mean_Effect_TFIID; Xlabel = 'mean Impact CAI TFIID';
%X = DS.median_Effect_TFIID; Xlabel = 'Effect of CAI TFIID';
%X = DS.mean_Effect_TATAless; Xlabel = 'Effect of CAI TATA_l_e_s_s';
X = DS.median_Effect_TATAless; Xlabel = 'tAI effect TATA-less';

%Y = DS.mean_Effect_SAGA; Ylabel = 'mean Impact CAI SAGA';
%Y = DS.median_Effect_SAGA; Ylabel = 'Effect of CAI SAGA';
%Y = DS.mean_Effect_TATA; Ylabel = 'Effect of CAI TATA';
Y = DS.median_Effect_TATA; Ylabel = 'tAI effect TATA';
%Y = DS.median_Effect_TATAless; Ylabel = 'Effect of CAI TATA_l_e_s_s';


fh = figure('units','centimeters','position',[5 5 8 8]); 
hold on;
%plot all exps:
plot(X,Y,'.','color',[.7 .7 .7],'markersize',.5)
%plot identity
minimum = min([X;Y]);
maximum = max([X;Y]);

%find and plot WT experiments:
SampleAts = SRAs_DS.sample_attribute;
idxWT = find(regexpcmp(SampleAts,'WT') | regexpcmp(SampleAts,'wt') | regexpcmp(SampleAts,'wild-type')); %cointains the SRRs that are WT
SRRs = SRAs_DS.SRR_accession(idxWT);
%plot each WT SRRs:
Xwt =[]; Ywt = []; 
for SRR = SRRs'
    idxSRR = find(strcmp(DS.ERR,SRR));
    if isempty(idxSRR)==0
        Xwt = [Xwt;X(idxSRR)];
        Ywt = [Ywt;Y(idxSRR)];
    end
end

Ywt_plot = Ywt(Ywt > 0 & Xwt > 0);
Xwt_plot = Xwt(Ywt > 0 & Xwt > 0);

%plot(Xwt,Ywt,'.','color','b','markersize',5)
% COLORS = flipud(bone(100));
% COLORS = COLORS(5:50,:);
% dscatter(Xwt_plot,Ywt_plot)
% colormap(COLORS)

%Add the DBP2 experiment:
SRA = 'SRA167760'; %DBP2 experiment
idxDBP2 = find(strcmp(SRAs_DS.SRA_accession,SRA));
SRRs = SRAs_DS.SRR_accession(idxDBP2);
Atribs = SRAs_DS.sample_attribute(idxDBP2);
idxDBP2_WT = find(regexpcmp(Atribs,'WT'));
idxDBP2_mut = find(regexpcmp(Atribs,'dbp2'));
%Plot WT:
for SRR = SRRs(idxDBP2_WT)'
    idxSRR = find(strcmp(DS.ERR,SRR));
    if isempty(idx)==0
        Xwt = X(idxSRR);
        Ywt = Y(idxSRR);
        plot(Xwt,Ywt,'o','color','k','markerfacecolor','k','markersize',7) 
    end
end

%Plot DBP2:
for SRR = SRRs(idxDBP2_mut)'
    idxSRR = find(strcmp(DS.ERR,SRR));
    if isempty(idx)==0
        Xwt = X(idxSRR);
        Ywt = Y(idxSRR);
        plot(Xwt,Ywt,'o','color','k','markerfacecolor',[0 0.5 0],'markersize',7) 
    end
end

plot(linspace(minimum,maximum,3),linspace(minimum,maximum,3),'-','color','b')

xlim([minimum,maximum])
ylim([minimum,maximum])
xlabel(Xlabel)
ylabel(Ylabel)

set(gca,'fontsize',11,'xtick',[0,0.1,0.2],'ytick',[0,0.1,0.2])

%% Plot DBP2 result in boxplots

%take the genes for which we have YFP data & CAI < 0.8
All_idx = find(DS_cor.Promoter_YFP > 0 & DS_cor.CAI <= 1 & DS_cor.tAI <= 0.5);
All_DS = DS_cor(All_idx,:);
r_threshold = 0;

%Go through each experiment (WT and DBP2)
SRA = 'SRA167760'; %DBP2 experiment
idxDBP2 = find(strcmp(SRAs_DS.SRA_accession,SRA));
SRRs = SRAs_DS.SRR_accession(idxDBP2);
Atribs = SRAs_DS.sample_attribute(idxDBP2);
idxDBP2_WT = find(regexpcmp(Atribs,'WT'));
idxDBP2_mut = find(regexpcmp(Atribs,'dbp2'));
%Generate WT data:
WTexp = [];
for SRR = SRRs(idxDBP2_WT)'
    idxSRR = find(strcmp(All_DS.Properties.VarNames,SRR));
    WTexp = [WTexp;double(All_DS(:,idxSRR))'];
end
WTexp = mean(WTexp)';
%Generate DBP2 data:
DBP2exp = [];
for SRR = SRRs(idxDBP2_mut)'
    idxSRR = find(strcmp(All_DS.Properties.VarNames,SRR));
    DBP2exp = [DBP2exp;double(All_DS(:,idxSRR))'];
end
DBP2exp = mean(DBP2exp)';

Count = 1;
for Gen = 1:2
    if Gen == 1; Data = WTexp; DataType = 'WT'; end
    if Gen == 2; Data = DBP2exp; DataType = 'DBP2ko'; end    
    DataType
    %Define the expression values:
    IDX = find(Data > 0 & Data < 1000000); 
    DS = All_DS(IDX,:);
    Expression = log10(Data(IDX));
    r = corr(log10(DS.Promoter_YFP),Expression);
    %Only consider experiment if it has a good correlation with expression:
    if r >= r_threshold
       %Get the r2 ratios (Sampling on #regulated genes, subtraction and 5x CV):       
        Sampling = 1; TypeSampling = 'Sampling #regulated genes'; EndJ = 10; %this is the number of samplings
        Ymat = []; %includes the increase in r2 values of each of the 4 groups
        Gmat = []; %includes the vector of the positions: 1(TFIID), 2(SAGA), 3(TATA-less), 4(TATA)
        for Class = 1:2
           if Class == 1
               idx1 = find(DS.TFIID == 1); 
               idx2 = find(DS.SAGA == 1);
               Labels = {'TFIID';'SAGA'};       
           end  
           if Class == 2
               idx1 = find(DS.TATA_type == 0); 
               idx2 = find(DS.TATA_type == 1);
               Labels = {'TATA_0';'TATA_1'};       
           end      
           Expression_1 = Expression(idx1);
           Expression_2 = Expression(idx2);
           Predictors_1 = log10(DS.Promoter_YFP(idx1));
           Predictors_2 = log10(DS.Promoter_YFP(idx2));
           %Predictors_CAI_1 = [log10(DS.Promoter_YFP(idx1)) , DS.CAI(idx1), DS.tAI(idx1)];
           Predictors_CAI_1 = [log10(DS.Promoter_YFP(idx1)) , DS.tAI(idx1)];   
           %Predictors_CAI_2 = [log10(DS.Promoter_YFP(idx2)) , DS.CAI(idx2), DS.tAI(idx2)];   
           Predictors_CAI_2 = [log10(DS.Promoter_YFP(idx2)) , DS.tAI(idx2)];   

           %Split sequences:
           RS = 100000; %Sample RS genes each time --> if you put a very high number you get all genes into it           
           for J = 1:EndJ 
              for TypeCAI = 1    
                %define the indexes (IDX) within each of the expression and predictors:
                %splitting into CAI ranges:
                %if TypeCAI == 1; IDX_1 = DS.CAI(idx1) <= 0.75; IDX_2 = DS.CAI(idx2) <= 0.75; end
                %if TypeCAI == 2; IDX_1 = DS.CAI(idx1) > 0.75; IDX_2 = DS.CAI(idx2) > 0.75; end   
                %without splitting:
                if TypeCAI == 1; IDX_1 = find(DS.CAI(idx1) <= 10); IDX_2 = find(DS.CAI(idx2) <= 10); end        
                if isempty(IDX_1)== 0 & isempty(IDX_2)==0        
                    %define each of the expression and predictor values
                    Exp_1 = Expression_1(IDX_1);
                    Exp_2 = Expression_2(IDX_2);
                    Pred_1 = Predictors_1(IDX_1,:);
                    Pred_2 = Predictors_2(IDX_2,:);
                    Pred_CAI_1 = Predictors_CAI_1(IDX_1,:);
                    Pred_CAI_2 = Predictors_CAI_2(IDX_2,:);
                    %define the indexes in which to make the sampling:
                    N1 = length(Exp_1); 
                    N2 = length(Exp_2);   
                    %if the number of samples to run by is low then take the maximum possible:
                    %random sampling:
                    if Sampling == 1
                        if RS > N1 | RS > N2; RS = min([N1;N2]); end
                        rs_idx1 = randperm(N1,RS)';
                        rs_idx2 = randperm(N2,RS)';
                    end
                    %without sampling:
                    if Sampling == 0
                        rs_idx1 = randperm(N1,N1)';
                        rs_idx2 = randperm(N2,N2)';
                    end
                    %do K fold cross validation:
                    %K = 5; Kfold = '5';
                    K = 2; Kfold = '2';
                    %constitutive genes (1):            
                    idx_cv = crossvalind('Kfold', length(Exp_1(rs_idx1)), K);
                    for cv = 1:K
                        idx_test = idx_cv == cv; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                        glm1 = fitglm( Pred_1(rs_idx1(idx_train)), Exp_1(rs_idx1(idx_train)));
                        glmCAI1 = fitglm( Pred_CAI_1(rs_idx1(idx_train),:),Exp_1(rs_idx1(idx_train)));
                        Estimate_tAI_TATAless = glmCAI1.Coefficients.Estimate(3)
                        Measured_1 = Exp_1(rs_idx1(idx_test));                
                        Predicted_1 = predict(glm1,Pred_1(rs_idx1(idx_test)));
                        Predicted_CAI_1 = predict(glmCAI1,Pred_CAI_1(rs_idx1(idx_test),:));
                        r2_1 = rsquare(Measured_1,Predicted_1);
                        r2CAI_1 = rsquare(Measured_1,Predicted_CAI_1);
                        %Dif_1 = log2(r2CAI_1/r2_1); Type_Comp = 'log_2 Ratio';                   
                        Dif_1 = r2CAI_1 - r2_1; Type_Comp = ' Subtraction';
                        %add to Ymat the values of each CV:
                        Ymat = [Ymat ; Dif_1 ];
                        %add the index to Gmat
                        if Class == 1; Gmat = [Gmat ; 1]; end
                        if Class == 2; Gmat = [Gmat ; 3]; end
                    end
                    %inducible genes (1):            
                    idx_cv = crossvalind('Kfold', length(Exp_2(rs_idx2)), K);
                    for cv = 1:K
                        idx_test = idx_cv == cv; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                        glm2 = fitglm( Pred_2(rs_idx2(idx_train)), Exp_2(rs_idx2(idx_train)));
                        glmCAI2 = fitglm( Pred_CAI_2(rs_idx2(idx_train),:),Exp_2(rs_idx2(idx_train)));
                        Estimate_tAI_TATA = glmCAI2.Coefficients.Estimate(3)
                        
                        Measured_2 = Exp_2(rs_idx2(idx_test));                
                        Predicted_2 = predict(glm2,Pred_2(rs_idx2(idx_test)));
                        Predicted_CAI_2 = predict(glmCAI2,Pred_CAI_2(rs_idx2(idx_test),:));
                        r2_2 = rsquare(Measured_2,Predicted_2);
                        r2CAI_2 = rsquare(Measured_2,Predicted_CAI_2);
                        %Dif_2 = log2(r2CAI_2/r2_2); Type_Comp = 'log_2 Ratio';
                        Dif_2 = r2CAI_2 - r2_2; Type_Comp = 'Subtraction';
                        %add to Ymat the values of each CV:
                        Ymat = [Ymat ; Dif_2 ];
                        %add the index to Gmat
                        if Class == 1; Gmat = [Gmat ; 2]; end
                        if Class == 2; Gmat = [Gmat ; 4]; end
                    end            

                end
              end   
           end
        end   
    end
    %save data for boxplot:
    if Gen == 1; YmatWT = Ymat; GmatWT = Gmat; end
    if Gen == 2; YmatDBP2 = Ymat; GmatDBP2 = Gmat; end    
    Count = Count + 1;
end

%% %Boxplots
Red = [0.8,0,0];
Blue = [39,101,186]./255;

colors = [Blue;Blue; Red; Red ; Blue;Blue; Red; Red];
%plot for each classification of promoters:
Col = 1;
for Class = 2
    if Class == 1; IDXs = [1,3]; end
    if Class == 2; IDXs = [2,4]; end
    fh = figure('units','centimeters','position',[5 5 11 7]); 
    hold on;
    plot(linspace(0,6,3),linspace(0,0,3),'-','color',[.7 .7 .7])
    for I = IDXs
        for Gen = 1:2
            if Gen == 1; Ymat = YmatWT; Gmat = GmatWT; SumToX = 0; end
            if Gen == 2; Ymat = YmatDBP2; Gmat = GmatDBP2; SumToX = 1;end
            idx = find(Ymat > -Inf & Ymat ~= 0 & Gmat == I);
            %define the values of boxplot:
            Y = Ymat(idx);
            G = Gmat(idx);
            %plot jittering:
            jitterAmount = 0.15;
            jitterValuesX = 2*(rand(size(G))-0.5)*jitterAmount;
            X = G + jitterValuesX + SumToX;
            plot(X,Y,'o','color',colors(I,:),'markerfacecolor',colors(I,:),'markersize',5)
            %boxplot
            b = boxplot(Y,G+SumToX,'color','k','positions',[G+SumToX],'width',0.6,'notch','on');
            set(b(7,:),'Visible','off')
            
            %save values
            I,Gen
            if I==2 & Gen==1; Y1_wt = Y; end
            if I==4 & Gen==1; Y2_wt = Y; end
            if I==2 & Gen==2; Y1_dbp2 = Y; end
            if I==4 & Gen==2; Y2_dbp2 = Y; end

            Col = Col + 1;
        end
    end
    Labels = {'WT','DBP2','WT','DBP2'};
    if Class == 1; Xtick = [1,2,3,4]; end
    if Class == 2; Xtick = [2,3,4,5]; end
    
    set(gca,'xtick',Xtick,'xticklabels',Labels)

    xlim([1,6])
    ylim([-0.2,0.4])
    ylabel('Effect of tAI')
    set(gca,'fontsize',11)

end
%%
%ttests:
[h,p] = ttest2(Y1_dbp2,Y1_wt)