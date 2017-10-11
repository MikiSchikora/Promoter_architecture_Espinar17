%% Description:
%This script creates the DEE data of CAI vs Expression

%% Load data

clear all
cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Miki/Data/');
%curated DEE datasets:
load('All_RNA_seq_S_cer_integrated_All_DS.mat')
load('All_RNA_seq_S_cer_integrated_PreDataSet.mat')

%% Create the dataset

%take the genes for which we have YFP data & CAI < 0.8
All_idx = find(DS_cor.Promoter_YFP > 0 & DS_cor.tAI <= 0.5);
All_DS = DS_cor(All_idx,:);

%create a DS that cointains, for each experiment, the information of CAI:
DEE_DS = dataset();
% define a threshold for which the experiment is going to be discarted
r_threshold = 0;

Row = 1;
%Go through each experiment:
for I = 1:length(Experiments)
    I
    %Define the expression values:
    IDX = find((double(All_DS(:,I)) > 0 & double(All_DS(:,I)) < 1000000)); 
    DS = All_DS(IDX,:);
    Expression = log10(double(DS(:,I)));
    r = corr(log10(DS.Promoter_YFP),Expression);
    %Only consider experiment if it has a good correlation with expression:
    if r >= r_threshold
       %Get the r2 ratios (Sampling on #regulated genes, subtraction and 5x CV):       
        Sampling = 1; TypeSampling = 'Sampling #regulated genes'; EndJ = 20; %this is the number of samplings
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
        
        %Add things to the dataset
        DEE_DS.ERR(Row,:) = All_DS.Properties.VarNames(I);
        DEE_DS.r_with_YFP(Row) = r;
        %For each of the types of promoters
        for G = 1:4
            Gidx = find(Gmat == G & Ymat ~= 0 & Ymat > -Inf);
            Y = Ymat(Gidx);
            %TFIID:
            if G == 1
                DEE_DS.mean_Effect_TFIID(Row) = mean(Y);
                DEE_DS.median_Effect_TFIID(Row) = median(Y);
                DEE_DS.sem_Effect_TFIID(Row) = sem(Y);
                DEE_DS.std_Effect_TFIID(Row) = std(Y);
            end            
            %SAGA:
            if G == 2
                DEE_DS.mean_Effect_SAGA(Row) = mean(Y);
                DEE_DS.median_Effect_SAGA(Row) = median(Y);
                DEE_DS.sem_Effect_SAGA(Row) = sem(Y);
                DEE_DS.std_Effect_SAGA(Row) = std(Y);
            end
            %TATA-less:
            if G == 3
                DEE_DS.mean_Effect_TATAless(Row) = mean(Y);
                DEE_DS.median_Effect_TATAless(Row) = median(Y);
                DEE_DS.sem_Effect_TATAless(Row) = sem(Y);
                DEE_DS.std_Effect_TATAless(Row) = std(Y);
            end    
            %TATA:
            if G == 4
                DEE_DS.mean_Effect_TATA(Row) = mean(Y);
                DEE_DS.median_Effect_TATA(Row) = median(Y);
                DEE_DS.sem_Effect_TATA(Row) = sem(Y);
                DEE_DS.std_Effect_TATA(Row) = std(Y);
            end            
        end
    Row = Row + 1;  
    end
end
    
% Save

cd('/Users/miki_schikora_tamarit/Google Drive/CareyLab/random RNA library/Promoter effect Manuscript/Data/');
save DEE_DS.mat DEE_DS