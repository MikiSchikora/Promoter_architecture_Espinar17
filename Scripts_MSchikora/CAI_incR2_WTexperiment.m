%% Description
% This script generate the increase in R2 upon including CAI in the model

%% Load data

clear all

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/ExternalData/');
%WT RNAseq expression. Column 4:
FitFlow_Fast_Slow = dataset('file','./vanDijk15/FitFlow__Fast_Slow.tab','delimiter','\t'); 
FitFlow_Fast_Slow.Properties.VarNames = {'ORF','Gene','Location','WTexpression','x18_2308'};
%5'end contribution to expression (use glucose):
End5_Contr = dataset('file','./Keren13/Keren13_Expression.tab','delimiter','\t'); 
End5_Contr.Properties.VarNames(3) = {'Promoter_YFP'};

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/LBC/2016_September__PromoterInfluencesRNALife__CodonBias_Predicts_Keren_Residuals/')
%DataSet of Keren (5'end contribution) with some adds:
End5_Contr_adds = dataset('file','keren_rnaseq_cai.tab','Delimiter','\t');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%ORF to seq
load('ORFtoSeq.mat');

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%TATA Types:
load('Promoter_Features.mat')
%TFIID/SAGA Types:
load('Yeast_SAGAvsTFIID.mat')


%% Join datasets:

%A dataset that cointains the information for all the native genes that are in each of the 5', 3' and WT genes.
NativeDS = join(End5_Contr,FitFlow_Fast_Slow,'Keys','ORF');
NativeDS = join(End5_Contr_adds,NativeDS,'Keys','ORF');
NativeDS = join(NativeDS,ORFseqDS,'Keys','ORF');


%% Add things:

%Correct Seq:
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq{I}; NativeDS.ORF_seq(I) = ORF;
end
%GC_content and length:
NativeDS.GC_content = 0;
NativeDS.ORF_length = 0;
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.Seq{I}; ORF = ORF{1};
    NativeDS.GC_content(I) = (length(regexp(ORF,'C')) + length(regexp(ORF,'G')))./length(ORF);
    NativeDS.ORF_length(I) = length(ORF);
end
%tAI:
NativeDS.tAI = cellfun(@calc_tAI , NativeDS.ORF_seq);
%SAGA vs TFIID:
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.ORF{I};
    idx = find(strcmp(SAGAvsTFIID_DS.ORF,ORF)); 
    if isempty(idx) == 0
        if strcmp(SAGAvsTFIID_DS.Type(idx),'TFIID-dominated') == 1
            NativeDS.TFIID(I) = 1;
        else
            NativeDS.TFIID(I,:) = 0;            
        end
        if strcmp(SAGAvsTFIID_DS.Type(idx),'SAGA-dominated') == 1
            NativeDS.SAGA(I) = 1;
        else
            NativeDS.SAGA(I) = 0;
        end       
    else
        NativeDS.TFIID(I) = 0;
        NativeDS.SAGA(I) = 0;
    end
end
%TATA feats:
for I = 1:length(NativeDS.ORF)
    ORF = NativeDS.ORF(I);
    idx = find(strcmp(AllPromotersDS.ORF,ORF));
    if isempty(idx) == 0;
        NativeDS.TATA_type(I) = AllPromotersDS.TATA_type(idx);
        NativeDS.Subclassification(I) = AllPromotersDS.Subclassification(idx);
    else
        NativeDS.TATA_type(I,:) = 10000;
        NativeDS.Subclassification(I) = 10000;        
    end    
end


%% Run the model with random sampling and splitting into >< 0.75 CAI --> Without 10x cross-validation
%Generate a matrix with the information for the  bar plot

Ymat = []; %includes the increase in r2 values of each of the 4 groups
Gmat = []; %includes the vector of the positions: 1(TFIID), 2(SAGA), 3(TATA-less), 4(TATA)

%define a dataset of overlapping CAI and expression
idx = find(NativeDS.CAI <= 0.8 & log10(NativeDS.WTexpression) > -10000);

for Class = 1:2
   %DS that cointains the data with several CAI values
   DS = NativeDS(idx,:);
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
   Expression_1 = log10(DS.WTexpression(idx1));
   Expression_2 = log10(DS.WTexpression(idx2));
   Predictors_1 = log10(DS.Promoter_YFP(idx1));
   Predictors_2 = log10(DS.Promoter_YFP(idx2));
   Predictors_CAI_1 = [log10(DS.Promoter_YFP(idx1)) , DS.CAI(idx1), DS.tAI(idx1)];
   Predictors_CAI_2 = [log10(DS.Promoter_YFP(idx2)) , DS.CAI(idx2), DS.tAI(idx2)];   
   
   %Split sequences:
   RS = 150; %Sample RS genes each time
   for J = 1:100 %Number of Random Iterations to run the model
      % Create a vector that takes predicted and measured values each time:
      Measured_1 = [];
      Measured_2 = [];       
      Predicted_1 = []; 
      Predicted_2 = []; 
      Predicted_CAI_1 = []; 
      Predicted_CAI_2 = [];      
      for TypeCAI = 1:2      
        %define the indexes (IDX) within each of the expression and predictors:
        %splitting into CAI ranges:
        %if TypeCAI == 1; IDX_1 = DS.CAI(idx1) <= 0.75; IDX_2 = DS.CAI(idx2) <= 0.75; end
        %if TypeCAI == 2; IDX_1 = DS.CAI(idx1) > 0.75; IDX_2 = DS.CAI(idx2) > 0.75; end   
        %without splitting:
        if TypeCAI == 1; IDX_1 = find(DS.CAI(idx1) <= 10); IDX_2 = find(DS.CAI(idx2) <= 10); end
        if TypeCAI == 2; IDX_1 = find(DS.CAI(idx1) > 10); IDX_2 = find(DS.CAI(idx2) > 10); end   
        
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
            if RS > N1 | RS > N2; RS = min([N1;N2]); end
            rs_idx1 = randperm(N1,RS)';
            rs_idx2 = randperm(N2,RS)';

            %constitutive genes (1):
            glm1 = fitglm( Pred_1(rs_idx1), Exp_1(rs_idx1));
            glmCAI1 = fitglm( Pred_CAI_1(rs_idx1,:),Exp_1(rs_idx1));
            Meas_1 = Exp_1(rs_idx1);
            Pred_1_s = glm1.Fitted.Response;
            Pred_CAI_1_s = glmCAI1.Fitted.Response;

            %inducible genes (2):
            glm2 = fitglm( Pred_2(rs_idx2), Exp_2(rs_idx2))  ;
            glmCAI2 = fitglm( Pred_CAI_2(rs_idx2,:),Exp_2(rs_idx2)) ; 
            Meas_2 = Exp_2(rs_idx2);
            Pred_2_s = glm2.Fitted.Response;
            Pred_CAI_2_s = glmCAI2.Fitted.Response;

            %each of the following vectors cointains the information to make the r2 calculations:
            Measured_1 = [Measured_1 ; Meas_1];
            Measured_2 = [Measured_2 ; Meas_2];       
            Predicted_1 = [Predicted_1 ; Pred_1_s];
            Predicted_2 = [Predicted_2 ; Pred_2_s];
            Predicted_CAI_1 = [Predicted_CAI_1 ; Pred_CAI_1_s];
            Predicted_CAI_2 = [Predicted_CAI_2 ; Pred_CAI_2_s];  
        end
    end
   %calculate r2 and differences:   
   r2_1 = rsquare(Measured_1,Predicted_1);
   r2_2 = rsquare(Measured_2,Predicted_2);
   r2CAI_1 = rsquare(Measured_1,Predicted_CAI_1);
   r2CAI_2 = rsquare(Measured_2,Predicted_CAI_2);
   Dif_1 = log2(r2CAI_1/r2_1);
   Dif_2 = log2(r2CAI_2/r2_2);
   
   %add to Ymat and Gmat:
   Ymat = [Ymat ; Dif_1 ; Dif_2];
   if Class == 1
       Gmat = [Gmat; 1 ; 2];
   end
   if Class == 2
       Gmat = [Gmat; 3 ; 4];
   end
   
   end
end

% Plot histogram of the 4 groups

idx = find(Ymat ~= 0);
Yval = Ymat(idx);
Gval = Gmat(idx);

Red = [0.8,0,0];

colors = [[0,0,1]; Red; [0,0,1] ; Red];
fh = figure('units','centimeters','position',[5 5 18 10]); 
hold on;
% Plot Gittered Values
for I = 1:4
    G = Gval(Gval == I);
    Y = Yval(Gval == I);
    jitterAmount = 0.1;
    jitterValuesX = 2*(rand(size(G))-0.5)*jitterAmount;
    X = G + jitterValuesX;
    plot(X,Y,'o','color','k','markerfacecolor',colors(I,:),'markersize',5)
end
%PLOT BOXPLOTS:
b = boxplot(Yval,Gval,'color','k','positions',[1,2,3,4],'width',0.6,'notch','on');
set(gca,'xtick',[1,2,3,4],'xticklabels',{'TFIID','SAGA','TATA-less','TATA'})
set(b(7,:),'Visible','off')
ylabel('r^2 increase upon CAI')
set(gca,'fontsize',14,'fontweight','bold')

%% Run the model without splitting: 10x CV 

Ymat = []; %includes the increase in r2 values of each of the 4 groups
Gmat = []; %includes the vector of the positions: 1(TFIID), 2(SAGA), 3(TATA-less), 4(TATA)

%define a dataset of overlapping CAI and expression
idx = find(NativeDS.tAI <= 0.5 & log10(NativeDS.WTexpression) > -10000);

for Class = 2
   %DS that cointains the data with several CAI values
   DS = NativeDS(idx,:);
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
   Expression_1 = log10(DS.WTexpression(idx1));
   Expression_2 = log10(DS.WTexpression(idx2));
   Predictors_1 = log10(DS.Promoter_YFP(idx1));
   Predictors_2 = log10(DS.Promoter_YFP(idx2));
   %Predictors_CAI_1 = [log10(DS.Promoter_YFP(idx1)) , DS.CAI(idx1), DS.tAI(idx1)];
   Predictors_CAI_1 = [log10(DS.Promoter_YFP(idx1)) , DS.tAI(idx1)];   
   %Predictors_CAI_2 = [log10(DS.Promoter_YFP(idx2)) , DS.CAI(idx2), DS.tAI(idx2)];   
   Predictors_CAI_2 = [log10(DS.Promoter_YFP(idx2)) , DS.tAI(idx2)];   
   
   %Split sequences:
   RS = 100000; %Sample RS genes each time --> if you put a very high number you get all genes into it
   Sampling = 1; TypeSampling = 'Sampling #regulated genes'; EndJ = 100;
   %Sampling = 0; TypeSampling = 'Taking all genes'; EndJ = 1;
   for J = 1:EndJ %Number of Random Iterations to run the model --> If we take all genes we put 1, always the same saple
      J
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
            K = 2; Kfold = '2';
            %K = 10; Kfold = '10';
            %constitutive genes (1):            
            idx_cv = crossvalind('Kfold', length(Exp_1(rs_idx1)), K);
            for I = 1:K
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                glm1 = fitglm( Pred_1(rs_idx1(idx_train)), Exp_1(rs_idx1(idx_train)));
                glmCAI1 = fitglm( Pred_CAI_1(rs_idx1(idx_train),:),Exp_1(rs_idx1(idx_train)));
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
            for I = 1:K
                idx_test = idx_cv == I; idx_train = ~idx_test; %define a training(90%) and a testing set(10%)
                glm2 = fitglm( Pred_2(rs_idx2(idx_train)), Exp_2(rs_idx2(idx_train)));
                glmCAI2 = fitglm( Pred_CAI_2(rs_idx2(idx_train),:),Exp_2(rs_idx2(idx_train)));
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

% Plot histogram of the 4 groups

idx = find(Ymat ~= 0);
%idx = find(Ymat > -Inf);

Yval = Ymat(idx);
Gval = Gmat(idx);
Blue = [39,101,186]./255;

Red = [0.8,0,0];

colors = [Blue; Red; Blue ; Red];
fh = figure('units','centimeters','position',[5 5 5 5]); 
hold on;
%line
plot(linspace(0,5,3),linspace(0,0,3),'-','color',[.7 .7 .7])
% Plot Gittered Values
for I = 3:4
    G = Gval(Gval == I);
    Y = Yval(Gval == I);
    jitterAmount = 0.15;
    jitterValuesX = 2*(rand(size(G))-0.5)*jitterAmount;
    X = G + jitterValuesX;
    hold on;
    if I == 3
        Y3 = Y;
    end
    if I == 4
        Y4 = Y;
    end    
    plot(X,Y,'o','color',colors(I,:),'markerfacecolor',colors(I,:),'markersize',5)
end

%PLOT BOXPLOTS:
b = boxplot(Yval,Gval,'color','k','positions',[3,4],'width',0.6,'notch','on');
set(gca,'xtick',[3,4],'xticklabels',{'TATA-less','TATA'})
set(b(7,:),'Visible','off')
ylabel('tAI effect')
set(gca,'fontsize',11)
%title(strcat(TypeSampling,{' '},Type_Comp,{' '},'Type CV:',{' '},Kfold))
[h,p] = ttest2(Y3,Y4)
