function make_motion_model(training_file_list, output_name, percent_var)    
       
% Creates shape and motion models from segmented examples.  The
% list of examples to use is stored in the training_file_list.
% This is a text file where the examples are listed in the
% following format:
% segmentation_file1.mat movie_file1.avi
% segmentation_file2.mat movie_file2.avi
% ...
% segmentation_filen.mat movie_filen.avi
% 
% The function creates shape model, which will be stored in
% output_name_ShapeModel.mat, and a motion model, which will be
% stored in the file output_name_MotionModel.mat
%
% The percent var parameter should take a value between 0 and 1
% (typically much closer to 1, e.g., 0.98) and represents the fraction of the
% variation present in the data that must be represented by the
% shape model.
    
    addpath('intersections');
    addpath('linlinintersect');
    % get manual segmentation results from input .mat file list
    file_list = fopen(training_file_list, 'r');
    data = textscan(file_list, '%s\t%s');
    seg_files = data{1};
    avi_files = data{2};
    fclose(file_list);
    
    ix = 0;
    
    for ii = 1:length(seg_files)
        [LeftEdgePoint1(ii,:), LeftEdgePoint2(ii,:)] = make_mask(avi_files{ii});
       seg_data{ii} = load(seg_files{ii});
       for f = 1:seg_data{ii}.state.NFRAME
           ix = ix + 1;           
           TrainingData(ix).Vertices = seg_data{ii}.state.XY(:,:,f);
            
           % Modify training data to force 1st point to lie on left edge of
           % ultrasound image
           xy_line = [LeftEdgePoint1(ii,:); LeftEdgePoint2(ii,:)];
          
           nsnaxels = length(TrainingData(ix).Vertices);
           [x_intersect, y_intersect, iout, jout] = intersections(xy_line(:,1), ...
                                                             xy_line(:,2),...
                                                             TrainingData(ix).Vertices(:,1),...
                                                             TrainingData(ix).Vertices(:,2),1);    
    
           if length(x_intersect) > 0        
               % cut snake at intersection point (remove snaxels that lie
               % before the edge)
               TrainingData(ix).Vertices = ...
                   TrainingData(ix).Vertices(floor(jout):end,:);
        
               % make first vertex of the snake be the intersection point
               TrainingData(ix).Vertices(1,:) = [x_intersect, ...
                                   y_intersect];
        
               % reinterpolate required number of snaxels along rest of
               % snake
               arclength = [0; ...
                            cumsum(sqrt(sum(diff(TrainingData(ix).Vertices).^2,2)))];
               inc = arclength(end)/(nsnaxels-1);
               TrainingData(ix).Vertices = interp1(arclength, ...
                                                     TrainingData(ix).Vertices, ...
                                                     0:inc:arclength(end), 'spline');
    
           else
               pt = linlinintersect([LeftEdgePoint1(ii,1), LeftEdgePoint1(ii,2);...
                                   LeftEdgePoint2(ii,1), LeftEdgePoint2(ii,2);...
                                   TrainingData(ix).Vertices(1,1),...
                                   TrainingData(ix).Vertices(1,2);...
                                   TrainingData(ix).Vertices(2,1),...
                                   TrainingData(ix).Vertices(2,2)]);
        
               % make first vertex of the snake be the intersection point
               TrainingData(ix).Vertices(1,:) = pt;
        
               % reinterpolate required number of snaxels along rest of
               % snake
               arclength = [0; ...
                            cumsum(sqrt(sum(diff(TrainingData(ix).Vertices).^2,2)))];
               inc = arclength(end)/(nsnaxels-1);
               TrainingData(ix).Vertices = interp1(arclength, ...
                                                   TrainingData(ix).Vertices, ...
                                                   0:inc:arclength(end), 'spline');
               
           end
       end
    end

    
    % Make the Shape model, which finds the variations between contours
    % in the training data sets and makes a PCA model describing typical
    % contours
            
    [ShapeData, TrainingData]= ASM_MakeShapeModel2D(TrainingData, ...
                                                    percent_var, 'first');
   
    ShapeModel_filename = [output_name,'_ShapeModel.mat'];
    save(ShapeModel_filename, 'ShapeData');

    % Make the Motion model, which looks at shape differences
    % between consecutive contours in the training data sets and
    % computes the covariance matrix of the shape parameters.
    
    state_diff = [];
    % get state variables for each frame in the database
    
    for ii=1:length(seg_files)      
        state_vec = zeros(3+size(ShapeData.Evectors,2), ...
                          seg_data{ii}.state.NFRAME);
        % state vector contents:
        % 1 & 2 - x, y position of leftmost part of curve
        % (constrained to lie along US
        % image edge)
        % 3 - length of curve
        % 4:end - shape parameters in ASM
                
        tongue_length = zeros(1,seg_data{ii}.state.NFRAME);
        normalized_pts = zeros(length(ShapeData.x_mean),seg_data{ii}.state.NFRAME);

        xy = [];
        xy = seg_data{ii}.state.XY;
                
        % fix xy data to constrain first point to lie on
        % US image edge
        xy_line = [LeftEdgePoint1(ii,:); LeftEdgePoint2(ii,:)];
        
        for f = 1:seg_data{ii}.state.NFRAME
            nsnaxels = length(xy(:,:,1));
            [x_intersect, y_intersect, iout, jout] = intersections(xy_line(:,1), ...
                                                              xy_line(:,2),...
                                                              xy(:,1,f),xy(:,2,f),1);    
            if length(x_intersect) > 0
                
                % cut snake at intersection point (remove snaxels that lie
                % before the edge)
                xy_new = xy(floor(jout):end,:,f);
        
                % make first vertex of the snake be the intersection point
                xy_new(1,:) = [x_intersect, y_intersect];
        
                % reinterpolate required number of snaxels along rest of
                % snake
                arclength = [0; ...
                             cumsum(sqrt(sum(diff(xy_new).^2,2)))];
                inc = arclength(end)/(nsnaxels-1);
                xy(:,:,f) = interp1(arclength, xy_new, 0:inc:arclength(end), 'spline');
    
            else
                pt = linlinintersect([LeftEdgePoint1(ii,1), LeftEdgePoint1(ii,2);...
                                    LeftEdgePoint2(ii,1), LeftEdgePoint2(ii,2);...
                                    xy(1,1,f),xy(1,2,f); xy(2,1,f),xy(2,2,f)]);
                
                % make first vertex of the snake be the intersection point
                xy(1,:,f) = pt;
        
                % reinterpolate required number of snaxels along rest of
                % snake
                arclength = [0; ...
                             cumsum(sqrt(sum(diff(xy(:,:,f)).^2,2)))];
                inc = arclength(end)/(nsnaxels-1);
                xy(:,:,f) = interp1(arclength,xy(:,:,f), 0:inc:arclength(end), 'spline');        
            end
        end               
                
        %state_vec(1:2,:) = mean(xy,1);
        state_vec(1:2,:) = xy(1,:,:);
                
        % convert xy position of first point to single
        % coordinate along US image edge                
        for f = 1:seg_data{ii}.state.NFRAME                    
            tongue_length(f) = sum(sqrt(sum(diff(xy(:,:, f)).^2,2)));
            normalized_pts(:,f) = reshape(bsxfun(@minus,xy(:,:,f),xy(1,:,f)),[],1);
        end

        state_vec(3,:) = tongue_length./tongue_length(1);

        normalized_pts = bsxfun(@rdivide, normalized_pts, tongue_length);
                
        normalized_pts = bsxfun(@minus,normalized_pts, ...
                                ShapeData.x_mean);
                
        for(f = 1:seg_data{ii}.state.NFRAME)
            state_vec(4:end,f) = normalized_pts(:,f)'*ShapeData.Evectors;
        end
                
        motion = diff(state_vec,1,2);
        motion(1:2,:) = bsxfun(@rdivide, motion(1:2,:), tongue_length(1:end-1));
        state_diff= [state_diff,motion];
    end 
             
    motion_model_var = std(state_diff,0,2)
    motion_model_corr_coef = corrcoef(state_diff')
    edge_point1 = LeftEdgePoint1;
    edge_point2 = LeftEdgePoint2;
            
    
    MotionModel_filename = [output_name,'_MotionModel.mat'];
    save(MotionModel_filename, 'motion_model_var', ...
         'motion_model_corr_coef', 'edge_point1', 'edge_point2');
                           

                


        
