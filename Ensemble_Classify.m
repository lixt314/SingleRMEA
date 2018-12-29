function [predicted_celltype] = Ensemble_Classify(data,training_celltype)
%addpath('F:\x86_64\alg')

New_Data2=data;
%New_Data2=New_Data1.New_Data;
%New_Data2 = normalizeData(New_Data2);
id = training_celltype~=0;
New_labels = training_celltype;%(id);
New_labels=New_labels+1;
New_labels1=New_labels';
% New_Singcell = New_Data2;%(id,:);
% New_Singcell_python = matpy.mat2nparray(New_Singcell);
% New_labels_python = matpy.mat2nparray(New_labels1);
% X_r=py.sklearn.discriminant_analysis.LinearDiscriminantAnalysis(pyargs('n_components',int32(300))).fit(New_Singcell_python,New_labels_python).transform(New_Singcell_python)
% X_r  = matpy.nparray2mat(X_r)
% New_Singcell=X_r;
    New_Singcell = New_Data2;%(id,:);
    New_Singcell_python = matpy.mat2nparray(New_Singcell);
    New_labels_python = matpy.mat2nparray(New_labels1);
    X_r=py.sklearn.discriminant_analysis.LinearDiscriminantAnalysis(pyargs('n_components',int32(300))).fit(New_Singcell_python,New_labels_python).transform(New_Singcell_python)
    X_r  = matpy.nparray2mat(X_r)
    New_Singcell_1=X_r; 
    New_Singcell = New_Data2;%(id,:);
    New_Singcell_python = matpy.mat2nparray(New_Singcell);
    New_labels_python = matpy.mat2nparray(New_labels1);
    X_r=py.sklearn.feature_selection.SelectKBest(pyargs('k',int32(100))).fit_transform(New_Singcell_python,New_labels_python)
    X_r  = matpy.nparray2mat(X_r)
    New_Singcell_2=X_r;
    New_Singcell = New_Data2;%(id,:);
    [features,weights] = MI(New_Singcell,New_labels,12);
    New_Singcell_3=New_Singcell(:,features(1:200));
    New_Singcell = New_Data2;%(id,:);
    [ranks,weights] = relieff(New_Singcell,New_labels',10);
    New_Singcell_4=New_Singcell(:,ranks(1:200));
    
for k=1:1:100
 cc=randi(4);
 if cc==1
    New_Singcell=New_Singcell_1;
  elseif cc==2  
   New_Singcell=New_Singcell_2;
  elseif cc==3 
   New_Singcell=New_Singcell_3;
 elseif cc==4 
   New_Singcell=New_Singcell_4;
 end
    
dd=randi(5);
if dd==1
    eclf=py.sklearn.ensemble.RandomForestClassifier(pyargs('n_estimators',int32(10)));
elseif dd==2
    eclf = py.sklearn.svm.SVC()
elseif dd==3
    eclf=py.sklearn.tree.DecisionTreeClassifier();    
elseif dd==4
    eclf=py.sklearn.naive_bayes.BernoulliNB()
elseif dd==5
    eclf=py.sklearn.neighbors.KNeighborsClassifier(pyargs('n_neighbors',int32(2))) 
end
indices = crossvalind('Kfold',New_labels,10);

    for i = 1:10 
        test = (indices == i); 
        train = ~test;
        NewTrain = New_Singcell(train,:);
        NewTest  = New_Singcell(test,:);
        y_train  = New_labels(train)';
        y_test   = New_labels(test)';

        NewTrain = matpy.mat2nparray(NewTrain);
        y_train  = matpy.mat2nparray(y_train);
        NewTest  = matpy.mat2nparray(NewTest);
    
        subLabels_ = eclf.fit(NewTrain, y_train).predict(NewTest);
        subLabels  = matpy.nparray2mat(subLabels_);
        celltype(test) = subLabels';
    end
    Fuction(k,:)=size(find(celltype==New_labels),2)./size(celltype,2);
    predicted_celltype(k,:)=celltype;
end


   D=100;
   M=2;
   p1 = [99 13  7  5  4  3  3  2  3];
   p2 = [ 0  0  0  0  1  2  2  2  2];
   p1 = p1(M-1);
   p2 = p2(M-1);
   [N,W] = F_weight(p1,p2,M)
    W(W==0) = 0.000001;
    T = 4;
    B = zeros(N);
    for i = 1 : N
        for j = i : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    MaxFES = 1000;
    Population = rand(N,D);
    Population1 = Population>0.5
    for k=1:1:N
        [~,Index]=find(Population1(k,:)==1);

        ensemble_labels=predicted_celltype(Index,:);
        for j=1:1:size(ensemble_labels,2)
            [ out ] = majorityvote(ensemble_labels(:,j));
            New_label(j)=out;
        end
        RRR(k,:)=New_label;
         FunctionValue(k,1) = -size(find(New_labels==New_label),2)./size(New_labels,2);
         FunctionValue(k,2) = sum(Population1(k,:));
    end

    delta = 0.9;
    Z = min(FunctionValue);
    Coding='Binary';
    Boundary=[ones(1,D); zeros(1,D)];
    Generations=10;
     nr = 2;
     F=0.5;
     CR=0.9
    for Gene = 1 : Generations
        %对每个个体执行操作
        for i = 1 : N
            %归一化
            Fmax = max(FunctionValue);
            Fmin = Z;
            FunctionValue = (FunctionValue-repmat(Fmin,N,1))./repmat(Fmax-Fmin,N,1);

                P = 1:N;
 
            k = randperm(length(P));

             Offspring = F_generator(Population(P(k(1)),:),Population(P(k(2)),:),Population(P(k(3)),:),Boundary, F, CR);
            Offspring1=Offspring;
       
            
              [~,Index]=find(Offspring1==1);

            ensemble_labels=predicted_celltype(Index,:);
        for j=1:1:size(ensemble_labels,2)
             [ out ] = majorityvote(ensemble_labels(:,j));
            New_label(j)=out;
        end
         OffFunValue(1) = -size(find(New_labels==New_label),2)./size(New_labels,2);
         OffFunValue(2) = sum(Offspring1);
         RX_off=New_label;
          OffFunValue = (OffFunValue-Fmin)./(Fmax-Fmin);
            
            %更新最优理想点
            Z = min(Z,OffFunValue);

            %更新邻居个体
            for j = 1 : T
                 g_old = sum(FunctionValue(B(i,j),:).*W(B(i,j),:));
                    g_new = sum(OffFunValue.*W(B(i,j),:));
                    %g_old = max(abs(FunctionValue(B(i,j),:)-Z).*W(B(i,j),:));
                    %g_new = max(abs(OffFunValue-Z).*W(B(i,j),:));
                if g_new < g_old
                    %更新当前向量的个体
                    Population(B(i,j),:) = Offspring;
                    FunctionValue(B(i,j),:) = OffFunValue;
                    RRR(B(i,j),:)=RX_off;
                end
            end

            %反归一化
            FunctionValue = FunctionValue.*repmat(Fmax-Fmin,N,1)+repmat(Fmin,N,1);

        end
        Gene=Gene+1;
      min(FunctionValue)  
    end

[final,Indexii]=min(FunctionValue)
predicted_celltype=RRR(Indexii(1),:)-1;
end





