function [ShapeData TrainingData]= ASM_MakeShapeModel2D(TrainingData, ...
                                                  percent_var, xy_ref_type)

if nargin < 2
    percent_var = 0.98;
end

if nargin < 3
    xy_ref_type = 'mean';
end

% Number of datasets (= total number of frames)
s=length(TrainingData);

% Number of landmarks
nl = size(TrainingData(1).Vertices,1);


%% Shape model

% Remove rotation and translation 
for i=1:s
      lengthaverage = 0;
      for j=1:size(TrainingData(i).Vertices(:,1),1)-1
        lengthaverage = lengthaverage + sqrt((TrainingData(i).Vertices(j+1,1)-TrainingData(i).Vertices(j,1))^2 + (TrainingData(i).Vertices(j+1,2)-TrainingData(i).Vertices(j,2))^2);     
      end
      
      if strcmpi(xy_ref_type, 'mean')
          TrainingData(i).Vertices(:,2) = (TrainingData(i).Vertices(:,2) - mean(TrainingData(i).Vertices(:,2)))./(lengthaverage);
          TrainingData(i).Vertices(:,1) = ...
              (TrainingData(i).Vertices(:,1)-mean(TrainingData(i).Vertices(:,1)))./(lengthaverage);
      elseif strcmpi(xy_ref_type, 'first')
          TrainingData(i).Vertices(:,2) = (TrainingData(i).Vertices(:,2) - TrainingData(i).Vertices(1,2))./(lengthaverage);
          TrainingData(i).Vertices(:,1) = ...
              (TrainingData(i).Vertices(:,1)-TrainingData(i).Vertices(1,1))./(lengthaverage);
      end
end

% Construct a matrix with all contour point data of the training data set
x=zeros(nl*2,s);
for i=1:length(TrainingData)
    x(:,i)=[TrainingData(i).Vertices(:,1)' TrainingData(i).Vertices(:,2)']';
end

[Evalues, Evectors, x_mean]=PCA(x);

% % Keep required % of all eigen vectors, (remove contour noise)
 i=find(cumsum(Evalues)>sum(Evalues)*percent_var,1,'first'); 
 Evectors=Evectors(:,1:i);
 Evalues=Evalues(1:i);

% Store the Eigen Vectors and Eigen Values
ShapeData.Evectors=Evectors;
ShapeData.Evalues=Evalues;
ShapeData.x_mean=x_mean;
ShapeData.x=x;


