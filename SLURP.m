function [xy,anchors] = SLURP(varargin)
%GETCONTOUR  - extract contour from movie frame
%
%	usage:  GetContour  % to initialize
%	        
%
% Open a movie via the menu
% Choose the reference frame via the slider
% click in order along contour to place anchor points
% click and drag on a point to reposition it
% modified click on a point deletes it
% double-click clears all points
%
% use the menu items to adjust anchor points and how the image is displayed
%
% mkt 02/12



if nargin == 0
    varargin{1} ='';
    varargin{2}= '';
end

% strip src,evt on callbacks
if ishandle(varargin{1}), varargin(1:2) = []; end;

% branch by action

switch varargin{1},

%-----------------------------------------------------------------------------
% ANCHORS:  redistribute anchor points

	case 'ANCHORS',
		state = get(gcbf,'userData');
                state.ANCHORS_VISIBLE = 1;
		nAnchors = size(state.ANCHORS,1);
		if nAnchors < 3, return; end;	% need at least three anchors
		k = [0 ; cumsum(sqrt(sum(diff(state.ANCHORS).^2,2)))];
		nAnchors = nAnchors + strmatch(varargin{2},{'DELETE','REDISTRIBUTE','INSERT'}) - 2;
		a = interp1(k,state.ANCHORS,linspace(0,k(end),nAnchors),'linear');
		delete(state.ALH);
		state.ALH = [];
		state.ANCHORS = a;                
		for k = 1 : nAnchors,
			state.ALH(k) = line(a(k,1),a(k,2),'marker','.','color','r','tag','CONTOUR','buttonDownFcn',{@SLURP,'DOWN','POINT'});
		end;
		set(gcbf,'userData',UpdateContour(state));
		

%-----------------------------------------------------------------------------
% CLOSE:  shutdown handler

	case 'CLOSE',
		state = get(gcbf,'userData');
		if ishandle(state.RH), delete(state.RH); end;		% close contrast adjustment (if open)
		delete(gcbf);										% close window		

%-----------------------------------------------------------------------------
% DOWN:  mouseDown handler
%
% click in image creates new anchor point
% click on an anchor repositions it
% modified click on a point deletes it
% double-click clears all points

	case 'DOWN',
		gotPoint = (length(varargin) ==2);	% nonzero for click on existing point       
        state = get(gcbf,'userData');
        state.ANCHORS_VISIBLE = 1;
		switch get(gcbf,'selectionType'),		

			case 'normal',			% unmodified

% move existing point
				if (gotPoint && state.CORREC) ,		
			
					set(gcbf, 'windowButtonMotionFcn',{@SLURP,'MOVE',gcbo}, ...
								'windowButtonUpFcn',{@SLURP,'UP',gcbo}, ...
								'pointer','crosshair');
% add new point								
				else, 
                   if (state.CORREC)                        
                    cp = get(gca, 'currentPoint');
                    cp = cp(1,1:2);
                    lh = line(cp(1),cp(2),'marker','.','color','r','tag','CONTOUR','buttonDownFcn',{@SLURP,'DOWN','POINT'});
                    state = get(gcbf,'userData');
                    state.ALH(end+1) = lh;
                    state.ANCHORS(end+1,:) = cp;
                    set(gcbf, 'windowButtonMotionFcn',{@SLURP,'MOVE',lh}, ...
                                'windowButtonUpFcn',{@SLURP,'UP',lh}, ...
                                'pointer','crosshair', ...
                                'userData',UpdateContour(state));
                            
                   end
				end;

% delete all points
			case 'open',			% double-click
				delete(findobj(gca,'tag','CONTOUR'));

% delete existing point
			otherwise,				% (ctl/shift) modified
				if gotPoint,
					state = get(gcbf,'userData');
					k = find(gcbo == state.ALH);
					state.ANCHORS(k,:) = [];
					state.ALH(k) = [];
					delete(gcbo);
					set(gcbf,'userData',UpdateContour(state));
				end;
		end;
		
	
%-----------------------------------------------------------------------------
% EMIT:  emit current contour to current frame

	case 'EMIT',
		state = get(gcbf,'userData');
		v = evalin('base',state.VNAME);
                %size(state.XY)
		v(:,:,state.CURFRAME) = state.XY(:,:,state.CURFRAME);
		assignin('base',state.VNAME,v);
		fprintf('current contour written to frame %d in %s\n',state.CURFRAME,state.VNAME);
						
		
%-----------------------------------------------------------------------------
% FLIP:  invert image

	case 'FLIP',
		if strcmp(varargin{2},'HORIZONTAL'),
			if strcmp(get(gca,'xdir'),'normal'), set(gca,'xdir','reverse'); else, set(gca,'xdir','normal'); end;
		else,
			if strcmp(get(gca,'ydir'),'normal'), set(gca,'ydir','reverse'); else, set(gca,'ydir','normal'); end;
		end;
		
		
%-----------------------------------------------------------------------------
% FRAME:  set image frame

	case 'FRAME',
		state = get(gcbf,'userData');
		f = state.CURFRAME;
		switch varargin{2},
			case 'PREV', if f > 1, f = f - 1; end;
			case 'NEXT', if f < state.NFRAME, f = f + 1; end;
			case 'TARG', f = state.TARGFRAME;
		end;
		if state.CURFRAME == f, return; end;
		set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
		set(state.TH,'string',sprintf('%04d',f));
		state.CURFRAME = f;
                set(state.CLH,'xdata',state.XY(:,1,state.CURFRAME),...
                              'ydata',state.XY(:,2,state.CURFRAME));
		set(gcbf,'userData',state);
		
	
%-----------------------------------------------------------------------------
% MAP:  set colormap

	case 'MAP',
		set(get(get(gcbo,'parent'),'children'),'checked','off');
		set(gcbo,'checked','on');
		mapName = varargin{2};
		if strcmp(mapName,'Inv Gray'),
			map = 1-gray;
		else,
			map = eval(lower(get(gcbo,'label')));
		end;
		set(gcbf,'colormap',map);
		

%-----------------------------------------------------------------------------
% MOVE:  mouseMvt handler

	case 'MOVE',
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);
		lh = varargin{2};
		x = get(lh,'xdata');
		y = get(lh,'ydata');
		if x ~= cp(1) || y ~= cp(2), 
			set(lh, 'xdata',cp(1), 'ydata',cp(2));
			state = get(gcbf,'userData');
			k = find(lh == state.ALH);
			state.ANCHORS(k,:) = [x,y];
			set(gcbf,'userData',UpdateContour(state));
		end;


%-----------------------------------------------------------------------------
% RANGE:  adjust range

	case 'RANGE',
		if strcmp(varargin{2},'FULL'),
			set(gca,'clim',[0 255]);
			return;
		end;
		state = get(gcbf,'userData');
		if ishandle(state.RH),
			figure(state.RH);
			return;
		end;
		state.RH = imcontrast(state.IH);
		set(gcbf,'userData',state);

%-----------------------------------------------------------------------------
% SETPARAM:  Set snake and tracking parameters
                
  case 'SETPARAM',
 
    
    state = get(gcbf, 'userData');
  prompt={'Image gradient smoothness (\sigma)', ...
         'Snake local search margin', ...
         'Snake internal energy weight (\alpha)', ...
         'Snake external energy weight (\beta)', ...
         'Snake internal smoothness weight (\lambda_1)', ...
         'Snake internal segment evenness weight (\lambda_2)', ...
         'Snake band penalty value', ...
         'Nparticles (use -1 for adaptive particle filtering)',...
         'Minimum number of particles (for adaptive particle filter tracking only)', ...
         ['Maximum number of particles (for adaptive particle filter ' ...
          'tracking only)']};
  
 
  name = 'Snake and particle filter tracking parameters';
  
  state.OPT
  
  defaultans = {num2str(state.OPT.Sigma), num2str(state.OPT.Delta), ...
                num2str(state.OPT.Alpha), num2str(1-state.OPT.Alpha), ...
                num2str(state.OPT.Lambda1), num2str(1-state.OPT.Lambda1), ...
                num2str(state.OPT.BandPenalty), num2str(state.OPT.Nparticles),...
                num2str(state.OPT.MinParticles), num2str(state.OPT.MaxParticles)};
  options.Interpreter = 'tex';
  answer = inputdlg(prompt,name,[1 40],defaultans, options);

  state.OPT.Sigma = str2double(answer{1});
  state.OPT.Delta = str2num(answer{2});
  state.OPT.Alpha = str2double(answer{3})/(str2double(answer{3})+ ...
                                        str2double(answer{4}));
  state.OPT.Lambda1 = str2double(answer{5})/(str2double(answer{5})+ ...
                                          str2double(answer{6}));
  state.OPT.BandPenalty = str2double(answer{7});
  state.OPT.Nparticles = str2num(answer{8});
  state.OPT.MinParticles = str2num(answer{9});
  state.OPT.MaxParticles = str2num(answer{10});
          
  set(gcbf, 'userData', state);

%-----------------------------------------------------------------------------
% SETMODEL:  set particle filter shape and motion models

  case 'SETMODEL',
    state = get(gcbf, 'userData');
    [filename,pathname]= uigetfile('*.mat',...
                        ['Select shape model ' ...
                        'from file...']);
    load(fullfile(pathname,filename));
        
    [filename,pathname] = uigetfile('*.mat',...
                         ['Select motion model ' ...
                        'from file...']);
    
    load(fullfile(pathname,filename));
    
    if exist('ShapeData') && exist('motion_model_var')
        if isempty(state.MODEL)
            state.MODEL = struct('Evectors', ShapeData.Evectors,...
                                 'Evalues', ShapeData.Evalues,...
                                 'x_mean', ShapeData.x_mean,...
                                 'motion_cov', kron(motion_model_var, ...
                                                    motion_model_var').* ...
                                 corr_coef);                                            
        else
            state.MODEL.Evectors = ShapeData.Evectors;
            state.MODEL.Evalues = ShapeData.Evalues;
            state.MODEL.x_mean = ShapeData.x_mean;
            state.MODEL.motion_cov = kron(motion_model_var,motion_model_var').* motion_model_corr_coef;    
        end
    else
        errordlg(['Selected files do not contain valid shape and/or ' ...
                  'motion models']);
    end
    
    set(gcbf, 'userData', state);
	
%-----------------------------------------------------------------------------
% SNAKE:  fit the snake

  case 'SNAKE',
    state = get(gcbf,'userData');
    
    state.ANCHORS_VISIBLE = 0;
    
        %Get Npoints
        hnpoints = findobj('tag', 'npoint');
        nPoints = get(hnpoints, 'String');
        state.NPOINTS = str2double(nPoints);
        state.XY = zeros(state.NPOINTS, 2, state.NFRAME);      
        set(gcbf, 'userData');
        
        %Get scale
        hnscale = findobj('tag', 'scale');
        scale = get(hnscale, 'string');
        state.SCALE = str2double(scale);
        set(gcbf, 'userData');
                               
        k = find(~ishandle(state.ALH));
        state.ANCHORS(k,:) = [];
        state.ALH(k) = [];
        anchors = state.ANCHORS;
        if size(anchors,1) < 4, return; end;
        set(gcbf,'pointer','watch'); drawnow;

        % compute image external energy (based on gradient) here to
        % save time
        Im = im2double(uint8(get(state.IH, 'cdata')));                
        Egradient = Egrad(Im, 5.0, [1, 1, size(Im,2), size(Im,1)]); 
        
        % interpolate anchor points for initial snake approximation (use
        % MATLAB's splines here rather than reprogram them in C)
        arclength = [0; cumsum(sqrt(sum(diff(anchors).^2,2)))];
        inc = arclength(end)/(state.NPOINTS-1);
        init_pts = interp1(arclength, anchors, 0:inc:arclength(end), ...
                           'spline');
        arclength = [0; cumsum(sqrt(sum(diff(init_pts).^2,2)))];
        state.LENGTHAVERAGE = arclength(end);
        set(gcbf, 'UserData',state);
          
        [xy, state.START_ENERGY] = make_snake(Im', ...
                                          Egradient',...
                                          init_pts, state.OPT.Delta*ones(state.NPOINTS,1), ...
                                          state.OPT.BandPenalty, state.OPT.Alpha, state.OPT.Lambda1, 1);
                                                    
        set(gcbf,'pointer','arrow');
        state.XY(:,:,state.CURFRAME) = xy;
        set(state.CLH,'xdata',xy(:,1),'ydata',xy(:,2));
        %set(state.ALH,'xdata',[],'ydata',[]);
        state.ALH=[];
        set(gcbf,'userData',state);
                
%-----------------------------------------------------------------------------
% TRACK: track the snake across the entire sequence

 case 'TRACK',
   state = get(gcbf,'userData');  
    state.ANCHORS_VISIBLE = 0;
    
    % Initialize energy matrice
    energy = zeros(state.NPOINTS, state.NFRAME);
    energy(:,state.CURFRAME) = state.START_ENERGY;
   
   % start is the reference frame 
    start = state.CURFRAME;
   % The reference frame's snake is used for the next frame 
    xy = state.XY(:,:,state.CURFRAME);
  
    % delta determines the research space
    delta = state.OPT.Delta*ones(state.NPOINTS,1);
    
    for frameInc = -1:2:1
        if frameInc < 0
            endFrame = 1;
        else
            endFrame = state.NFRAME-1;
        end
        for f = start:frameInc:endFrame
            frame = f;           
            set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
            set(state.TH,'string',sprintf('%04d',f));
            state.CURFRAME = f;
            set(state.CLH,'xdata',state.XY(:,1,f+1),...
                          'ydata',state.XY(:,2,f+1));
            
            % compute image external energy (based on gradient) here to
            % save time
            Im = im2double(uint8(get(state.IH,'cdata')));
            Egradient = Egrad(Im, 5.0, ...
                              [1, 1, size(Im,2), size(Im,1)]);
            
            % Make the snake    
            [xy_new, energy(:,f)] = make_snake(Im', ...
                                               Egradient', ...
                                               xy, delta, state.OPT.BandPenalty, ...
                                               state.OPT.Alpha, state.OPT.Lambda1, 1);
            
            xy = xy_new;
            

            % Save the snake
            state.XY(:,:,f) = xy;
            
    
            % Set the contour
            set(state.CLH,'xdata',state.XY(:,1,f),...
                          'ydata',state.XY(:,2,f));
            
            %update drawing
            set(gcbf,'pointer','watch'); drawnow;   
            
        end
    end
   
  set(gcbf,'userData',state);
    
  %Display energy map
  fig =  figure('name', 'energy');
  
  %Normalize energy
  for i=1:state.NFRAME
    energy_visu(1:state.NPOINTS-1,i) =100*(abs((energy(1:state.NPOINTS-1,i)-energy(1:state.NPOINTS-1,start)))./energy(1:state.NPOINTS-1,start));
  end 

   state.ENERGY(:,:)= energy_visu(:,:);
   state.ENERGY(state.NPOINTS,:)= energy_visu(state.NPOINTS-1,:);
   
  % Save the reference energy at the end of the array
  state.ENERGY(:,state.NFRAME+1) = energy(:,start);
   
  set(gcbf,'userData',state); 
  % Display energy map
  imagesc(energy_visu(1:state.NPOINTS-1,1:state.NFRAME), [0 400]);
  colormap('gray');
  xlabel('frame no.');
  ylabel('snaxel no.');
  
  % Attached main figure to energy map's axes.
  setappdata(gca,'fig',gcbf);
  %get cursor
  dcmObj = datacursormode(fig);
  %Modify  UpdateFcn of the cursor (When a user clicks on the map)
  set(dcmObj,'UpdateFcn',@myupdatefcn,'Enable','on');
  

%-----------------------------------------------------------------------------
% PARTICLE:  Track using particle filter

  case 'PARTICLE',
    state = get(gcbf,'userData');
    
    if isempty(state.MODEL)
        errordlg(['No shape and motion models loaded.  Please load ' ...
                  'suitable shape and motion models before using the ' ...
                  'particle filter, or use tracking without particle ' ...
                  'filtering.']);
    else
    
    frameInc = 1;
    endFrame = state.NFRAME-1;
    
    s = RandStream('mcg16807','Seed',0);
    RandStream.setGlobalStream(s);
    state.ANCHORS_VISIBLE = 0;   
  
    energy = zeros(state.NPOINTS, state.NFRAME);
    energy(:,state.CURFRAME) = state.START_ENERGY;
   
    Evectors = state.MODEL.Evectors; 
    Evalues = state.MODEL.Evalues; 
    x_mean = state.MODEL.x_mean;
    motion_cov = state.MODEL.motion_cov;
        
    % start is the reference frame 
    start = state.CURFRAME;
    % The reference frame's snake is used for the next frame 
    xy = state.XY(:,:,state.CURFRAME);
    start_length = sum(sqrt(sum(diff(xy).^2,2)));
    
    
    
    % get state vector for initial frame
    for frameInc = -1:2:1        
        if state.OPT.Nparticles < 0
            max_particles = state.OPT.MaxParticles;
        else
            max_particles = state.OPT.Nparticles;
        end
        min_particles = state.OPT.MinParticles;
        Nparticles = max_particles;    
        pfstate = zeros(3+size(Evectors,2),Nparticles);
        if frameInc < 0 
            endFrame = 1;
        else
            endFrame = state.NFRAME-1;
        end
        xy = state.XY(:,:,start);
        pfstate(1:2,:) = repmat(xy(1,:)',1,Nparticles);    
        normalized_pts = reshape((xy - repmat(pfstate(1:2,1)',length(xy),1)),[],1)./ ...
            start_length;
        pfstate(3,:) = ones(1,Nparticles); % length scale wrt original
                                           % tongue length
        
        pfstate(4:end,:) = repmat(((normalized_pts - x_mean)'*Evectors)',1,Nparticles);
        
        % get asm contour from projected particle states
        for pp=1:Nparticles
            pt_vec = (x_mean + sum(repmat(pfstate(4:end,pp)', ...
                                          length(Evectors),1).*Evectors,2));
            pf_xy(:,:,pp) = reshape(pt_vec, length(xy), 2);
            pf_xy(:,:,pp) = pf_xy(:,:,pp).* pfstate(3,pp).*start_length + ...
                repmat(pfstate(1:2,pp)',length(xy),1);
        end                               

        for f = start+frameInc:frameInc:endFrame
            set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
            set(state.TH,'string',sprintf('%04d',f));
            state.CURFRAME = f;
            set(state.CLH,'xdata',state.XY(:,1,f-frameInc),...
                          'ydata',state.XY(:,2,f-frameInc));
            
            % sample new particles using state transition model
            % for the moment, assume uncorrelated gaussian noise for
            % the 5 variables                
            old_state = pfstate;
            
            rv = mvnrnd(zeros(1,size(pfstate,1)), motion_cov, ...
                        max_particles);
            
            pfstate(1:2,:) = pfstate(1:2,:) + start_length* ...
                repmat(pfstate(3,:),[2 1]).*...
                rv(:,1:2)';
            
            pfstate(3:end,:) = pfstate(3:end,:) + rv(:,3:end)';                
            
            % evaluate new particles using external energy terms in
            % snake model
            % pre-compute image gradient
            Im = im2double(uint8(get(state.IH,'cdata')));
            Egradient = Egrad(Im, 5.0, [1, 1, size(Im,2), size(Im,1)], state.MASK);
            
            ImTrans = Im';
            EgradientTrans = Egradient';
            
            pf_xy = zeros(length(xy), 2, max_particles);
            pxy = zeros(size(pf_xy));
            p_energy = zeros(length(xy),max_particles);
            
            pp = 0;
            cumlike = 0;
            minlike = 7*exp(-sum(state.START_ENERGY));
            like_thresh = minlike;        
            
            while (pp < max_particles-1 && (cumlike < like_thresh || ...
                                            state.OPT.Nparticles > 0)) || pp < min_particles-1
                pp  = pp + 1;
                
                % get snaxel positions from particle state vector            
                pt_vec = (x_mean + sum(repmat(pfstate(4:end,pp)',length(Evectors),1).*Evectors,2));
                pf_xy(:,:,pp) = reshape(pt_vec, length(xy), 2);
                pf_xy(:,:,pp) = pf_xy(:,:,pp).* pfstate(3,pp).*start_length + ...
                    repmat(pfstate(1:2,pp)',length(xy),1);           
                
                [pxy(:,:,pp), p_energy(:,pp)] = ...
                    make_snake(ImTrans, ...
                               EgradientTrans, ...
                               pf_xy(:,:,pp), state.OPT.Delta*ones(state.NPOINTS,1), ...
                               state.OPT.BandPenalty, state.OPT.Alpha, state.OPT.Lambda1, 0);
                
                Eext(pp) = sum(p_energy(:,pp));
                
                len = sum(sqrt(sum(diff(pxy(:,:,pp)).^2,2)));
                lratio(pp) = max(len/start_length, start_length/len);
                %lratio(pp) = 1;
                
                cumlike = cumlike + exp(-Eext(pp)*lratio(pp));
            end               
            Nparticles = pp + 1;
            
            % compute particle weights from external energy
            like(1:Nparticles-1) = exp(-Eext(1:Nparticles-1).*lratio(1:Nparticles-1));
            
            % show curve corresponding to most likely particle
            [max_w(f), index] = max(like(1:Nparticles-1));
            
            % get snaxel positions from processed particle state vector            
            xy = pxy(:,:,index);                              
            [xy, energy(:,f)] = ...
                make_snake(ImTrans, ...
                           EgradientTrans, xy, state.OPT.Delta*ones(state.NPOINTS,1), ...
                           state.OPT.BandPenalty, state.OPT.Alpha, state.OPT.Lambda1, 1);         
            
            
            len = sum(sqrt(sum(diff(xy).^2,2)));
            lr = max(len/start_length, start_length/len);
            
            % save updated best state to new particle
            pfstate(:,Nparticles) = zeros(size(pfstate,1),1);
            pfstate(1:2,Nparticles) = xy(1,:);
            pfstate(3,Nparticles) = sum(sqrt(sum(diff(xy).^2,2)))./start_length;
            normalized_pt = reshape((xy(:,:)-...
                                     repmat(pfstate(1:2,Nparticles)',length(xy(:,:)),1)),...
                                    [],1)./(pfstate(3,Nparticles)*start_length);        
            pfstate(4:end,Nparticles) = (normalized_pt - x_mean)'*Evectors;              
            
            saved_state(:,f-frameInc) = pfstate(:,Nparticles);
            
            % get asm contour from projected particle state
            best_pt_vec = (x_mean + sum(repmat(pfstate(4:end,Nparticles)', ...
                                               length(Evectors),1).*Evectors,2));
            best_pf_xy = reshape(best_pt_vec, length(xy(:,:)), 2);
            best_pf_xy = best_pf_xy.* pfstate(3, Nparticles).*start_length + ...
                repmat(pfstate(1:2,Nparticles)',length(xy),1);  
            
            
            state.XY(:,:,f) = xy;
            set(state.CLH,'xdata',state.XY(:,1,f),'ydata',state.XY(:,2, ...
                                                              f));                 
            set(gcbf,'pointer','watch'); drawnow;
            set(gcbf,'userData',state);
            like(Nparticles) = exp(-sum(energy(:,f)*lr));
            weight = like(1:Nparticles)./sum(like(1:Nparticles)); 
            
            cdf = cumsum(weight);
            cdf_prev = circshift(cdf,[0 1]);
            cdf_prev(1) = 0.0;
            samples = rand(max_particles,1);
            for pp = 1:max_particles        
                try,
                    ix(pp) = find(cdf >= samples(pp) & cdf_prev <= ...
                                  samples(pp));
                catch,
                    samples(pp)
                    cdf
                    pause;
                end;
            end
            
            new_pfstate = zeros(size(pfstate,1),max_particles);
            new_pfstate = pfstate(:,ix(1:max_particles));
            
            pfstate = new_pfstate;        
            
        end
    end
    
    fig =  figure('name', 'energy');
    %Normalize energy
  for i=1:state.NFRAME
      energy_visu(1:state.NPOINTS-1,i) = 100*((energy(1: ...
                                                      state.NPOINTS-1,i)-energy(1:state.NPOINTS-1,start))./energy(1:state.NPOINTS-1,start));
      energy_visu = energy_visu .* (energy_visu >= 0);
  end
  
  state.RAW_ENERGY = energy;
  state.ENERGY(:,:)= energy_visu(:,:);
  state.ENERGY(state.NPOINTS,:)= energy_visu(state.NPOINTS-1,:);
   
  % Save the reference energy at the end of the array
  state.ENERGY(:,state.NFRAME+1) = energy(:,start);
   
  set(gcbf,'userData',state); 
  % Display energy map
  imagesc(energy_visu(1:state.NPOINTS-1,1:state.NFRAME), [0 400]);
  colormap('gray');
  xlabel('frame no.');
  ylabel('snaxel no.');
  
  % Attached main figure to energy map's axes.
  setappdata(gca,'fig',gcbf);  
  dcmObj = datacursormode(fig);
  %Modify  UpdateFcn of the cursor (When an user click on the map)
  set(dcmObj,'UpdateFcn',@myupdatefcn,'Enable','on');
    end
      
%-----------------------------------------------------------------------------
% UP:  mouseUp handler

	case 'UP',
		set(gcbf, 'windowButtonMotionFcn','', 'windowButtonUpFcn','', 'pointer','arrow');
		cp = get(gca, 'currentPoint');
		cp = cp(1,1:2);
		lh = varargin{2};
		x = get(lh,'xdata');
		y = get(lh,'ydata');
		if x ~= cp(1) || y ~= cp(2), 
			set(lh, 'xdata',cp(1), 'ydata',cp(2));
			state = get(gcbf,'userData');
                        state.ANCHORS_VISIBLE = 1;
			k = find(lh == state.ALH);
			state.ANCHORS(k,:) = [x,y];
			set(gcbf,'userData',UpdateContour(state));
		end;
        

    
%-----------------------------------------------------------------------------
% UNDO:  Undo correction


  case 'UNDO'
    %t = toc;
    %state.USER_TIME = state.USER_TIME + t;
        state = get(gcbf,'userData');
        state.ANCHORS_VISIBLE = 0;
        
        for f=1:state.NFRAME
            state.CURFRAME = f;
            state.XY(:,:,f) = state.OLDANCHORS(:,:,f);
            set(state.CLH,'xdata',state.XY(:,1, ...
                 state.CURFRAME),'ydata',state.XY(:,2,state.CURFRAME));
            
            set(gcbf,'pointer','watch'); drawnow;
            
        end
        set(gcbf,'userData',state)
		state.OLDANCHORS = [];
                %tic;
        

%-----------------------------------------------------------------------------
% CORRECT:  Correct frames as edgetrack do
  case 'CORRECTASEDGETRACK',   % initialize a frame and track the
                               % snake from this frame
                               %t = toc;    
   % If gotClicked is true, means the user clicked on the image
    gotClicked = false;   
    state = get(gcbf,'userData');   
    %state.USER_TIME = state.USER_TIME + t;
    
    if(length(varargin) >1)
        % Means the user clicked on the image
        gotClicked = true;        
    end
      
    
    if(state.CORRECTED) % Means the user has choosed 3 points --> The correction is over, we want to do a new correction 
       % Initialize the points used for correction 
       state.PTS_ANCRAGE = [0 0];
       % Don't run the 'DOWN' case --> Don't add new points to the snake
       state.CORRECTED = false;
       findobj('Type', 'figure', 'name', 'CorrecParameters')
       % If the parameters's figure doesn't exist
       if (isempty(findobj('Type', 'figure', 'name', 'CorrecParameters')));
          % Create the new figure
           paramFig = figure('position',...
                   [520 380 300 75],...
                  'units', 'pixels',...
                   'name', 'CorrecParameters');
            edit1 = uicontrol('style', 'edit', ...
                'units', 'pixels', ...
                'position',[2 50 150 20],...
                'tag', 'frame',...
                'String', num2str(state.CURFRAME),...
                'Enable', 'Inactive',...
                'ButtonDownFcn', @editfunction);
           uicontrol('style', 'text', ...
                 'units', 'pixels', ...
                 'position',[160 50 150 20],...
                 'String', 'start');
           edit2 = uicontrol('style', 'edit', ...
                'units', 'pixels', ...
                'position',[2 5 150 20],...
                'tag', 'frame',...
                'String', num2str(state.CURFRAME),...
                'Enable', 'Inactive',...
                'ButtonDownFcn', @editfunction);  
            uicontrol('style', 'text', ...
                'units', 'pixels', ...
                'position',[160 5 150 20],...
                'String', 'end');
            
            set(edit1, 'UserData', struct('otherEdit', edit2, 'selected',0));
            set(edit2, 'UserData', struct('otherEdit', edit1, 'selected',0))
       
        % Initalize the parameters for correction with the CURFRAME
        state.PARAM_CORREC.FBEF = str2double(get(edit1, 'String'));
        state.PARAM_CORREC.FAFT = str2double(get(edit2, 'String'));

       else
       % If the window exists, put it on foreground
       handle = findobj('Type', 'figure', 'name', 'CorrecParameters');  
       figure(handle);
           
       end
       set(gcbf,'userData',state)
    end
    
    if(gotClicked) % If the user clicked on the image
        state.CORREC = false;  
        cp = get(gca, 'currentPoint'); %get point
        cp = cp(1,1:2);
        line(cp(1),cp(2),'marker','+','color','g'); 
        state.PTS_ANCRAGE(end+1,:) = [cp(1) cp(2)];  %save the point 
        set(gcbf,'userData',state);
        
        if(size(state.PTS_ANCRAGE,1) == 4) %If 3 points has been choosed 
            state.CORRECTED = true; 
            endup = false;
            enddown = false; 
            cp =  state.PTS_ANCRAGE(end,:);
            f = state.CURFRAME;     
            % Calculate the distance between the point 1 and 2 and the snake
            % gives minup/down = distance between the point and the snake
            % gives departup/down = the closest snaxel from the correction point
            [mindown,departdown] = min(sqrt((state.XY(:,1,f)-state.PTS_ANCRAGE(2,1)).^2+(state.XY(:,2,f)-state.PTS_ANCRAGE(2,2)).^2));
            [minup,departup] = min(sqrt((state.XY(:,1,f)-state.PTS_ANCRAGE(3,1)).^2 +(state.XY(:,2,f)-state.PTS_ANCRAGE(3,2)).^2));                   
            if( minup > 10 )
                endup = true; % The selected point is not on the snake
            end
            
            if (mindown> 10)
                enddown = true; % The selected point is not on the snake
            end
            
            line(cp(1),cp(2),'marker','.','color','b');           
            delta = state.OPT.Delta*ones(state.NPOINTS,1);
            state.OLDANCHORS(:,:,:) = state.XY(:,:,:); % Save points before correction
            state.ANCHORS = [];
            state.ALH = [];                      
            npoints = size(state.XY,1);



            fdown = state.PARAM_CORREC.FBEF;  %Start frame for the correction
            fup = state.PARAM_CORREC.FAFT;    % End frame for the correction
            energy(:,:) = zeros(state.NPOINTS, state.NFRAME);
            
            % Initialize the snake of the curframe
            start = state.CURFRAME;
            f = state.CURFRAME;   
            set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));  
            set(state.TH,'string',sprintf('%04d',f));
            state.OLDANCHORS(:,:,f) = state.XY(:,:,f); % Save points before correction
            if (endup) % If endup, cut the snake from departdown
                state.ANCHORS(:,:) = cat(1,state.XY(1:departdown,:,f) ,[cp(1) cp(2)], [state.PTS_ANCRAGE(3,1) state.PTS_ANCRAGE(3,2)]);
            elseif (enddown) % If endown, cut the snake until departup
                state.ANCHORS(:,:) = cat(1,[state.PTS_ANCRAGE(2,1) state.PTS_ANCRAGE(2,2)], [cp(1) cp(2)],state.XY(departup:npoints,:,f));
            else    % else, cut the snake between departdown and departup
                state.ANCHORS(:,:) = cat(1,state.XY(1:departdown,:,f) ,[cp(1) cp(2)],state.XY(departup:npoints,:,f));                            
            end  
            k = [0 ; cumsum(sqrt(sum(diff(state.ANCHORS).^2,2)))];
            xy = interp1(k,state.ANCHORS,linspace(0,k(end),state.NPOINTS),'pchip');
            
            Im = im2double(uint8(get(state.IH,'cdata')));
            Egradient = Egrad(Im, 5.0, [1, 1, size(Im,2), size(Im,1)]);
            [xy_new, energy(:,f)] = make_snake(Im', ...
                                               Egradient', ...
                                               xy, delta, ...
                                               state.OPT.BandPenalty, state.OPT.Alpha, ...
                                               state.OPT.Lambda1, 1);                              
            state.XY(:,:,f) = xy_new(:,:); % Save  the new snake

            set(state.CLH,'xdata',state.XY(:,1, ...
                                                   f),'ydata',state.XY(:,2,f));


            set(gcbf,'pointer','watch'); drawnow;
            state.CORREC = false;   
            set(gcbf,'userData',state);
            
            % Start the tracking from the reference frame to fup
            for f=start+1:fup
                        state.ANCHORS = [];
                        state.ALH = [];
                        state.CURFRAME = f;
                        xy = state.XY(:,:,f-1);
                        set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));  
                        set(state.TH,'string',sprintf('%04d',f)); ...
                            
                        Im = im2double(uint8(get(state.IH,'cdata')));
                        Egradient = Egrad(Im, 5.0, [1, 1, size(Im,2), size(Im,1)]);
                        [xy_new, energy(:,f)] = make_snake(Im', ...
                                                           Egradient', ...
                                                           xy, delta, ...
                                                           state.OPT.BandPenalty, ...
                                                           state.OPT.Alpha, ...
                                                           state.OPT.Lambda1, 1);                   
                        state.XY(:,:,f) = xy_new(:,:);    

                        set(state.CLH,'xdata',state.XY(:,1, ...
                                                               state.CURFRAME),'ydata',state.XY(:,2,state.CURFRAME));
                        set(gcbf,'pointer','watch'); drawnow;
                        state.CORREC = false;   
                        set(gcbf,'userData',state);

            end
           
           % Start the tracking from the reference frame to fdown
           for f=start-1:-1:fdown
                        state.ANCHORS = [];
                        state.ALH = [];
                        state.CURFRAME = f;
                        xy = state.XY(:,:,f+1);
                        set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));  
                        set(state.TH,'string',sprintf('%04d',f));                            
                        Im = im2double(uint8(get(state.IH,'cdata')));
                        Egradient = Egrad(Im, 5.0, [1, 1, size(Im,2), size(Im,1)]);
                        [xy_new, energy(:,f)] = make_snake(Im', ...
                                                           Egradient', ...
                                                           xy, delta, ...
                                                           state.OPT.BandPenalty, ...
                                                           state.OPT.Alpha, ...
                                                           state.OPT.Lambda1, 1);                   
                        state.XY(:,:,f) = xy_new(:,:);    

                        set(state.CLH,'xdata',state.XY(:,1, ...
                                                              f),'ydata',state.XY(:,2,f));


                        set(gcbf,'pointer','watch'); drawnow;
                        state.CORREC = false;   
                        set(gcbf,'userData',state);

            end
            
            % Save the old energy in the matrice energy_visu
            energy_visu(:,:) = state.ENERGY(:,1:state.NFRAME);           
            % Modify the old energy with the energy between fdown and fup
            for i=fdown:fup
                % Normalize the energy with the energy from the reference
                % frame (saved in state.ENERGY(:,state.NFRAME+1))
                energy_visu(1:state.NPOINTS-1,i) = 100*(abs((energy(1:state.NPOINTS-1,i)-state.ENERGY(1:state.NPOINTS-1,state.NFRAME+1)))./state.ENERGY(1:state.NPOINTS-1,state.NFRAME+1));
            end
            
            % Save the new energy in the variable state
            state.ENERGY(:,1:state.NFRAME) = energy_visu(:,1:state.NFRAME);
            set(gcbf,'userData',state);
            
            %Find the current map of energy
            figure_energy = findobj('Type', 'figure', 'name', 'energy'); 
            ax1 = get(figure_energy, 'CurrentAxes');
            
            %Create a new figure with the new energy
            figure_energy_new= figure;
            imagesc(energy_visu(1:state.NPOINTS-1,:),[0 400]); 
            colormap('gray');
            xlabel('frame no.');
            ylabel('snaxel no.');
            ax2 = gca;
            
            %Create a new figure to show the 2 maps
            h3 = figure('name', 'energy');
 
            s1 = subplot(1,2,1); %create and get handle to the subplot axes
            s2 = subplot(1,2,2);
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            fig2 = get(ax2,'children');
            copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
            copyobj(fig2,s2);
            colormap(gray)
            
            close(figure_energy);
            close(figure_energy_new);
            colormap('gray');
            xlabel(s1,{'frame no.'});
            ylabel(s1,{'snaxel no.'});
            title(s1 ,{'old energy'});
            xlabel(s2,{'frame no.'});
            ylabel(s2,{'snaxel no.'});
            title(s2 ,{'new energy'});
            set(s1, 'ydir', 'reverse');
            set(s2,'ydir', 'reverse');
            set(s1, 'Clim', [0 400]);
            set(s2, 'Clim', [0 400]);           
            % Attached main figure to energy map's axes.
            setappdata(s1,'fig',gcbf);
            % Attached main figure to energy map's axes.
            setappdata(s2,'fig',gcbf);
            %get cursor
            dcmObj = datacursormode(h3);
            %Modify  UpdateFcn of the cursor (When an user click on the map)
            set(dcmObj,'UpdateFcn',@myupdatefcn,'Enable','on');
            
            
        end
        
    
    else  
        % When the user click on the image, the function
        % @SLURP,'CORRECTASEDGETRACK','CLICK' is called
        set(gcbf,'WindowButtonDownFcn', {@SLURP,'CORRECTASEDGETRACK','CLICK'}); 
        % state.CORREC is false --> When the user click on the image, the
        % function @GetContour,'CORRECTASEDGETRACK','CLICK' is called, and
        % the function @GetContour,'DOWN','POINT' can't be run
        state.CORREC = false; 
       
    end
    %tic
%-----------------------------------------------------------------------------
% REV:  Reverse Tracking for correction

    case 'REV'
       
        StartTracking = false;     %StartTracking is true when the user push "Track" 
        
        if(length(varargin) == 2)
            StartTracking = true;        
        end    
        
        
        
        if (StartTracking)
            %t = toc;        
            main_fig=get(gcbf, 'UserData')
            close(gcbf);
            state = get(main_fig,'userData');
            delta = state.OPT.Delta*ones(state.NPOINTS,1);
            xy = state.XY(:,:,state.PARAM_TRACK.FBEF);
            energy(:,:) = zeros(state.NPOINTS, state.NFRAME);
            state.OLDANCHORS(:,:,:) = state.XY(:,:,:); % Save points before correction
            for f = state.PARAM_TRACK.FBEF-1:-1:state.PARAM_TRACK.FAFT
             
                set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
                set(state.TH,'string',sprintf('%04d',f));
                state.CURFRAME = f;
                set(state.CLH,'xdata',state.XY(:,1,f+1),...
                              'ydata',state.XY(:,2,f+1));




                % compute image external energy (based on gradient) here to
                % save time
                Im = im2double(uint8(get(state.IH, 'cdata')));
                Egradient = Egrad(Im, 5.0, ...
                                  [1, 1, size(Im,2), size(Im,1)]);

                [xy_new, energy(:,f)] = make_snake(Im', ...
                                                   Egradient', ...
                                                   xy, delta, ...
                                                   state.OPT.BandPenalty, state.OPT.Alpha, ...
                                                   state.OPT.Lambda1, 1);

  

                xy = xy_new;
  
               state.XY(:,:,f) = xy;
                    
               set(state.CLH,'xdata',state.XY(:,1,f),...
                  'ydata',state.XY(:,2,f));
               
              set(main_fig,'pointer','watch'); drawnow;  
         
            end

            set(main_fig,'userData',state);
            
            energy_visu(:,:) = state.ENERGY(:,1:state.NFRAME);           
            for i=state.PARAM_TRACK.FAFT:state.PARAM_TRACK.FBEF-1            
                energy_visu(1:state.NPOINTS-1,i) = 100*(abs((energy(1:state.NPOINTS-1,i)-state.ENERGY(1:state.NPOINTS-1,state.NFRAME+1)))./state.ENERGY(1:state.NPOINTS-1,state.NFRAME+1));
            end
            
            state.ENERGY(:,1:state.NFRAME) = energy_visu(:,1:state.NFRAME);
            set(main_fig,'userData',state);
            
            
            %interesting_figures = findobj('Type', 'figure', 'name', ...
            %                             'energy');
            
            % figure_energy_no = max(interesting_figures.Numbers);
            
            %figure_energy = findobj('Type', 'figure', 'number', figure_energy_no);
            figure_energy = findobj('Type', 'figure', 'name', 'energy'); 
            ax1 = get(figure_energy, 'CurrentAxes');
            
            
            figure_energy_new= figure;
            imagesc(energy_visu(1:state.NPOINTS-1,:),[0 400]); 
            colormap('gray');
            xlabel('frame no.');
            ylabel('snaxel no.');
            ax2 = gca;
            
            h3 = figure('name', 'energy');
 
            s1 = subplot(1,2,1); %create and get handle to the subplot axes
            s2 = subplot(1,2,2);
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            fig2 = get(ax2,'children');
            copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
            copyobj(fig2,s2);
            
            close(figure_energy);
            close(figure_energy_new);
            colormap('gray');
            xlabel(s1,{'frame no.'});
            ylabel(s1,{'snaxel no.'});
            title(s1 ,{'old energy'});
            xlabel(s2,{'frame no.'});
            ylabel(s2,{'snaxel no.'});
            title(s2 ,{'new energy'});
            set(s1, 'ydir', 'reverse');
            set(s2,'ydir', 'reverse');
              set(s1, 'Clim', [0 400]);
              set(s2, 'Clim', [0 400]);     
            % Attached main figure to energy map's axes.
            setappdata(s1,'fig',main_fig);
            % Attached main figure to energy map's axes.
            setappdata(s2,'fig',main_fig);
            %get cursor
            dcmObj = datacursormode(h3);
            %Modify  UpdateFcn of the cursor (When an user click on the map)
            set(dcmObj,'UpdateFcn',@myupdatefcn,'Enable','on');
            
            %tic;
            
        else
            % GUI to choose parameters
            state = get(gcbf,'userData');
            main_fig = gcbf; 
            if (isempty(findobj('Type', 'figure', 'name', 'ReverseTracking')));
                figure('position',...
                       [520 380 300 100],...
                      'units', 'pixels',...
                       'name', 'ReverseTracking',...
                       'UserData', main_fig);
                edit1 = uicontrol('style', 'edit', ...
                    'units', 'pixels', ...
                    'position',[2 75 150 20],...
                    'tag', 'trackFrame',...
                    'String', num2str(state.CURFRAME),...
                    'Enable', 'Inactive',...
                    'ButtonDownFcn', @editfunction);
               uicontrol('style', 'text', ...
                     'units', 'pixels', ...
                     'position',[160 75 150 20],...
                     'String', 'start');
               edit2 = uicontrol('style', 'edit', ...
                    'units', 'pixels', ...
                    'position',[2 50 150 20],...
                    'tag', 'trackFrame',...
                    'String', num2str(state.CURFRAME),...
                    'Enable', 'Inactive',...
                    'ButtonDownFcn', @editfunction);  
                uicontrol('style', 'text', ...
                    'units', 'pixels', ...
                    'position',[160 50 150 20],...
                    'String', 'end');
                uicontrol('style', 'pushbutton', ...
                    'units', 'pixels', ...
                    'position',[2 5 100 20],...
                    'tag', 'done',...
                    'String', 'Track ', ...
                    'callback', {@SLURP, 'REV','Track'}); 

                %To know which edit is selected 
                set(edit1, 'UserData', struct('otherEdit', edit2, 'selected',0));
                set(edit2, 'UserData', struct('otherEdit', edit1, 'selected',0));
                
               %FBEF = Reference Frame
               state.PARAM_TRACK.FBEF = str2double(get(edit1, 'String'));
               %FATF = Track to that frame
               state.PARAM_TRACK.FAFT = str2double(get(edit2, 'String'));
               
               
                          
               set(gcbf,'userData',state)
                
            end
                    

        end        
        
%-----------------------------------------------------------------------------
% SAVE:  Save initialization
   case 'SAVE'     
   state = get(gcbf, 'userData');
   initialization = struct('anchors', state.ANCHORS(:,:), ...
       'frame', state.CURFRAME);
   
   [filename, pathname] = uiputfile('.mat', 'Save Initialization As','Ini');

    if ~(isequal(filename,0) || isequal(pathname,0))
       save(fullfile(pathname, filename), 'initialization')
    end
    
    
%-----------------------------------------------------------------------------
% SAVESTATE:  Save state
  case 'SAVESTATE'     
   state = get(gcbf, 'userData');
   fName = state.VIDEO_NAME;
   [filename, pathname] = uiputfile('.mat', 'Save State As',fName);

    if ~(isequal(filename,0) || isequal(pathname,0))
       save(fullfile(pathname, filename), 'state')
    end    
   
    %tic;
%-----------------------------------------------------------------------------
% LOADSTATE:  Load state
  case 'LOADSTATE'     
    state = get(gcbf, 'userData');
   [filename, pathname] = uigetfile('.mat', 'Load state...');
   if ~(isequal(filename,0) || isequal(pathname,0))
       old_state = load(fullfile(pathname, filename), 'state');
       if isfield(old_state.state, 'NPOINTS')
           state.NPOINTS = old_state.state.NPOINTS;    
       end
       %state.ANCHORS = old_state.state.ANCHORS;
       state.XY = old_state.state.XY;
       if ~isfield(old_state.state, 'NPOINTS')
           state.NPOINTS = length(state.XY(:,1,1));
       end
       state.ANCHORS = state.XY(:,:,state.CURFRAME);
       size(state.ANCHORS)
       if isfield(old_state.state, 'OPT')
           state.OPT = old_state.state.OPT;  
       end
       if isfield(old_state.state, 'SCALE')
           state.SCALE = old_state.state.SCALE;
       end
       
       state.ENERGY = old_state.state.ENERGY;
       
       if isfield(old_state.state, 'START_ENERGY')
           state.START_ENERGY = old_state.state.START_ENERGY;
       else
           state.START_ENERGY = state.ENERGY(:,1);
       end
       
       clear old_state;
       set(state.CLH,'xdata',state.XY(:,1, state.CURFRAME),...
                     'ydata',state.XY(:,2, state.CURFRAME));       

       set(gcbf, 'userData', UpdateContour(state));
       
   end   

%-----------------------------------------------------------------------------
% ENERGYMAP:  View energy map
%Display energy map
  case 'ENERGYMAP'
  state = get(gcbf,'userData');
  fig =  figure('name', 'energy');  
  % Display energy map
  size(state.ENERGY)
  state.NPOINTS
  state.NFRAME
  
  imagesc(state.ENERGY(1:state.NPOINTS-1,1:state.NFRAME), [0 400]);
  colormap('gray');
  xlabel('frame no.');
  ylabel('snaxel no.');
  % Attached main figure to energy map's axes.
  setappdata(gca,'fig',gcbf);
  %get cursor
  dcmObj = datacursormode(fig);
  %Modify  UpdateFcn of the cursor (When an user click on the map)
  set(dcmObj,'UpdateFcn',@myupdatefcn,'Enable','on');    
  
%-----------------------------------------------------------------------------
% SAVEPTS:  Save points in a text file    
   case 'SAVEPTS'     
   state = get(gcbf, 'userData');
   
   [filename, pathname] = uiputfile('.con', 'Save Points As','SEG');

    if ~(isequal(filename,0) || isequal(pathname,0))
       xy = state.XY;
       distr = computer
       if (strcmp(distr,'PCWIN64') || strcmp(distr,'PCWIN'))
           l = 'pc';
       else
           l = 'unix';
       end
        dlmwrite(fullfile(pathname, filename), xy, 'delimiter','\t', ...
            'precision', 5, 'newline', l);
    end            
    
%-----------------------------------------------------------------------------
% LOAD:  Load Anchors
    case 'LOAD' 
      
    [filename, pathname] = uigetfile('*.mat', 'Pick a mat file');
    if ~(isequal(filename,0) || isequal(pathname,0))
       load(fullfile(pathname, filename))
       state = get(gcbf,'userData');       
       state.CURFRAME = initialization.frame;
       set(state.IH,'cdata',uint8(mean(read(state.MH,state.CURFRAME),3)));
       set(state.TH,'string',sprintf('%04d',state.CURFRAME));
       state.ANCHORS =  initialization.anchors;
       for k = 1 : length(initialization.anchors),
           state.ALH(k) = line(initialization.anchors(k,1),initialization.anchors(k,2),...
                               'marker','.','color','r','tag','CONTOUR','buttonDownFcn',{@SLURP,'DOWN','POINT'});
       end;              
       state.ANCHORS_VISIBLE = 1;
       set(gcbf,'userData',UpdateContour(state));
      
    end       

%-----------------------------------------------------------------------------
% INTER:  Interpolate the curve to have more points
    case 'INTER' 
        StartInter = false;     %StartInter is true when the user push "Inter" 
        
        if(length(varargin) == 2)
            StartInter = true;        
        end    
        
        
        
        if (StartInter)
            
            main_fig=get(gcbf, 'UserData');
            hnpoints = findobj('tag', 'Npoints');
            Npoints = str2double(get(hnpoints, 'String'));
            close(gcbf);
            state = get(main_fig,'userData');            
            XYnew = zeros(Npoints, 2, state.NFRAME);
            for f = 1:state.NFRAME
                state.CURFRAME = f;
                arclength = [0; cumsum(sqrt(sum(diff(state.XY(:,:,state.CURFRAME)).^2,2)))];
                inc = arclength(end)/(Npoints-1);
                xy = state.XY(:,:,state.CURFRAME);
                XYnew(:,:,f) = interp1(arclength, xy...
                    , 0:inc:arclength(end), 'spline');
                set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
                set(state.TH,'string',sprintf('%04d',f));
                set(state.CLH,'xdata',XYnew(:,1,f),...
                                  'ydata',XYnew(:,2,f));
            end
            state.XY = XYnew;
            state.NPOINTS = Npoints;
            set(main_fig,'userData',state);
            
            
            
        else                  
            % GUI to choose parameters
            state = get(gcbf,'userData');
            main_fig = gcbf; 
            if (isempty(findobj('Type', 'figure', 'name', 'InterpolateSnake')));
                figure('position',...
                       [520 380 300 100],...
                      'units', 'pixels',...
                       'name', 'InterpolateSnake',...
                       'UserData', main_fig);
                edit1 = uicontrol('style', 'edit', ...
                    'units', 'pixels', ...
                    'position',[2 75 150 20],...
                    'tag', 'Npoints',...
                    'String', num2str(state.NPOINTS));
               uicontrol('style', 'text', ...
                     'units', 'pixels', ...
                     'position',[160 75 150 20],...
                     'String', 'Npoints');
                uicontrol('style', 'pushbutton', ...
                    'units', 'pixels', ...
                    'position',[2 5 100 20],...
                    'tag', 'done',...
                    'String', 'Interpolate ', ...
                    'callback', {@SLURP, 'INTER','Done'}); 
                       
               set(gcbf,'userData',state);
            end
        end

    
    
    
    
%-----------------------------------------------------------------------------
% OPEN:  Open a video


    case 'OPEN'
        
        

        close(gcbf)
        currentFigure = Initialize();
        
        [p,f,e] = uigetfile('*.avi', 'Choose your video');
        
        
        
        % defaults
        anchors = [];
        nPointsDef = 39;
        scaleDef = 98/321;
        vName = '';
        opt = struct('Sigma', 5.0, 'Delta', 2, 'BandPenalty', 2.0, ...
                     'Alpha', 0.8, 'Lambda1', 0.95, 'MinParticles', ...
                     10, 'MaxParticles', 1000, 'Nparticles', -1); 

        
        % look for default shape and motion models.
        model = [];
        default_shape_filename = fullfile(pwd,'ShapeModel.mat');
        default_motion_filename = fullfile(pwd,'MotionModel.mat');
        if exist(default_shape_filename) == 2 && ...
                exist(default_motion_filename) ==2
            load(default_shape_filename);
            load(default_motion_filename);
            
            if exist('ShapeData') && exist('motion_model_var')
                model = struct('Evectors', ShapeData.Evectors,...
                               'Evalues', ShapeData.Evalues,...
                               'x_mean', ShapeData.x_mean,...
                               'motion_cov', kron(motion_model_var, ...
                                                  motion_model_var').*motion_model_corr_coef);               
                disp('Default shape and motion models loaded');
            end
        end
        
        
        if ~(isequal(p,0) || isequal(f,0))
            fName = fullfile(f,p);
            [pathn, fn, en] = fileparts(fName);
            try,                
                mh = VideoReader(fName);
            catch,
                error('unable to open %s', fName);
            end;
            try,
                                                
                frame = 1;
                img = read(mh,frame);
                nFrames = get(mh,'numberOfFrames');
                for ii = 1:min(nFrames,200)
                    img_data(:,:,ii) = double(rgb2gray(read(mh,ii)));
                end                               
                
                % find image mask by looking at parts where there is variation from one 
                % image to the next (the background should be 99% the same between all
                % images since they come from the same machine
                var_img = var(img_data, 0, 3);
                mask = var_img > 2;
                
                
                CC = bwconncomp(mask);
                numOfPixels = cellfun(@numel,CC.PixelIdxList);
                [unused,indexOfMax] = max(numOfPixels);
                unmask = zeros(size(mask));
                unmask(CC.PixelIdxList{indexOfMax}) = 1;                
                
                mask = unmask;                                
                
        catch,
            error('unable to load frame %d in %s', frame, fName);
        end;

            
            % display image
            img = uint8(mean(img,3));
            ih = imshow(img);
            th = text(50,80,sprintf('%04d',frame),'fontsize',18,'color','w');                            
        
        % create emit array
            if ~isempty(vName),
                uimenu(menu,'label','Emit Contour','separator','on','accelerator','E','callback',{@SLURP,'EMIT',1});
                assignin('base',vName,NaN(nPointsDef,2,nFrames,'single'));
                fprintf('created %s [%d contour points x X,Y x %d movie frames] in base ws\n',vName,nPointsDef,nFrames);
            end;

            


            % init parameters for frames correction
            param_correc = struct('FBEF', 0,...
                                 'FAFT', 0,...
                                 'TRESHENER', 0,...
                                 'AUTO', true);

            % init parameters for reverse tracking
            param_tracking = struct('FBEF', 0,...
                                 'FAFT', 0);

            % init internal state
            state = struct('IH', ih, ...				% image handle
                            'MH', mh, ...				% movie handle
                            'TH', th, ...				% text handle
                            'TARGFRAME', frame, ...		% target frame
                            'CURFRAME', frame, ...		% currently displayed frame
                            'NFRAME', nFrames, ...		% number of available frames
                           'MASK', mask, ...
                            'RH', [], ...				% contrast adjustment handle
                            'NPOINTS', nPointsDef, ...		% number of contour points
                            'ANCHORS', anchors, ...		% current anchor points
                            'ALH', [], ... % and their line handles
                            'ANCHORS_VISIBLE', 1,... 
                            'XY', [], ... % current contour points                     
                            'CLH', [], ...  % and their line handle
            'FWD', [],...
                'FWD_CLH', [],...
                'BWD', [],...
                'BWD_CLH',[],...
                            'OLDANCHORS', [],...        % contour points before correction
                            'SCALE', scaleDef,...          % Image scale
                            'VNAME', vName,...          % emit array name
                            'CORREC', true,...          % if false, "Corrected frames" has been called 
                            'ENERGY', [],...            
                            'START_ENERGY', [], ...
                            'LENGTHAVERAGE', 1, ...    % Length of the first snake
                            'PTS_ANCRAGE', [0 0],...   % Corrected points
                            'CORRECTED', true,...      % if true, "3 corrected points has been chosen 
                            'PARAM_CORREC', param_correc,...
                            'PARAM_TRACK' , param_tracking, ...
                           'USER_TIME', 0, ...
                           'OPT', opt,...
                           'MODEL', model,...
                           'VIDEO_NAME', [pathn,'/',fn]);


            state.XY = zeros(state.NPOINTS, 2, state.NFRAME);
            
           
            % Slide to choose the reference frame
            position = get(currentFigure, 'Position');
            width = position(3);            
            maxframe = state.NFRAME;
            hSlider = uicontrol('Style', 'slider',...
                'userData', currentFigure,...
                'Max', maxframe, ...
                'Min', 1, ...
                'SliderStep', [1 1]/(maxframe-1), ...
                'Value', state.CURFRAME,...
                'Position', [5 5 width/2 20],...
                'Callback', @slider_targframe);
            
            % Edit with the number of point
            uicontrol('Style', 'text', ...
                'String', 'Npoint', ...
                'Position',  [width/2+20 25 50 20])
            hnpoint = uicontrol('Style', 'edit',...
                'Value', nPointsDef,...
                'tag', 'npoint', ...
                'String', num2str(nPointsDef),...
                'Position', [width/2+20 5 50 20]);
            
           % Edit with the scale
            uicontrol('Style', 'text', ...
                'String', 'Scale', ...
                'Position',  [width/2+100 25 50 20])
            hnscale = uicontrol('Style', 'edit',...
                'Value', scaleDef,...
                'tag', 'scale', ...
                'String', num2str(scaleDef),...
                'Position', [width/2+100 5 50 20]);           
                                    
            set(currentFigure,'name',sprintf('%s  frame %03d',f,frame), ...
                    'busyAction','cancel', ...
                    'closeRequestFcn',{@SLURP,'CLOSE'}, ...
                    'tag','SLURP', ...
                    'userData',state);
            set(ih,'buttonDownFcn', {@SLURP,'DOWN'});
            

            % initialize contour from passed-in anchor points
            if ~isempty(anchors),
                for k = 1 : size(anchors,1),
                    state.ALH(k) = line(anchors(k,1),anchors(k,2),'marker','.','color','r','tag','CONTOUR','buttonDownFcn',{@SLURP,'DOWN','POINT'});
                end;
                set(currentFigure,'userData',UpdateContour(state)); 
            end;
            if isunix, [s,r] = unix('osascript -e ''tell application "MATLAB" to activate'''); end;
        
             
        end
        

        
        

	
%-----------------------------------------------------------------------------
% INITialize

    otherwise
              
        Initialize(varargin{:});
end;


%===============================================================================
% INITIALIZE

%function Initialize(fName, frame, scale, varargin)
function fh = Initialize(varargin)


fh = figure;

% init menu
menu = uimenu(fh, 'label','SLURP');

uimenu(menu,'label','Open a video', ...
			'accelerator', 'B', ...
			'callback',{@SLURP,'OPEN'});
                
uimenu(menu,'label','Previous Frame', ...
			'accelerator', 'B', ...
			'callback',{@SLURP,'FRAME','PREV'});
uimenu(menu,'label','Next Frame', ...
			'accelerator', 'F', ...
			'callback',{@SLURP,'FRAME','NEXT'});
uimenu(menu,'label','Target Frame', ...
			'accelerator', 'T', ...
			'callback',{@SLURP,'FRAME','TARG'});
uimenu(menu,'label','Insert Anchor', ...
			'separator','on', ...
			'accelerator','I', ...
			'callback',{@SLURP,'ANCHORS','INSERT'});
uimenu(menu,'label','Delete Anchor', ...
			'accelerator','D',...
			'callback',{@SLURP,'ANCHORS','DELETE'});
uimenu(menu,'label','Redistribute Anchors', ...
			'accelerator','R',...
			'callback',{@SLURP,'ANCHORS','REDISTRIBUTE'});
uimenu(menu,'label','Fit Snake', 'callback',{@SLURP,'SNAKE', ...
                    1});
uimenu(menu,'label','Track snake without particle filter', 'callback',{@SLURP,'TRACK'});

h = uimenu(menu,'label','Track snake using particle filter', 'callback', {@SLURP,'PARTICLE'});

uimenu(menu, 'label', 'Set snake and tracking parameters', 'callback', ...
       {@SLURP,'SETPARAM'});

uimenu(menu, 'label', 'Set particle filter models', 'callback', ...
       {@SLURP, 'SETMODEL'});

uimenu(menu,'label','Interpolate snake ', 'callback',{@SLURP,'INTER'});

h = uimenu(menu,'label','Contrast','separator','on');
uimenu(h,'label','Full Range','callback',{@SLURP,'RANGE','FULL'});
uimenu(h,'label','Adjust Range','callback',{@SLURP,'RANGE','ADJUST'});

h = uimenu(menu,'label','Colormap');
uimenu(h,'label','Gray','checked','on','callback',{@SLURP,'MAP','Gray'});
uimenu(h,'label','Inv Gray','callback',{@SLURP,'MAP','Inv Gray'});
cMaps = {'Autumn','Bone','Colorcube','Cool','Copper','Flag','Hot','HSV','Jet','Lines','Pink','Prism','Spring','Summer','Winter'};
for k = 1 : length(cMaps),
	if k == 1, sep = 'on'; else, sep = 'off'; end;
	uimenu(h,'label',cMaps{k},'separator',sep,'callback',{@SLURP,'MAP',cMaps{k}});
end;

h = uimenu(menu,'label','Flip');
uimenu(h,'label','Horizontal', ...
			'callback',{@SLURP,'FLIP','HORIZONTAL'});
uimenu(h,'label','Vertical', ...
			'callback',{@SLURP,'FLIP','VERTICAL'});
        
uimenu(menu,'label','Save initialization', ...
			'accelerator', 'H', ...
			'callback',{@SLURP,'SAVE'});

 uimenu(menu,'label','Save State', ...
			'accelerator', 'I', ...
			'callback',{@SLURP,'SAVESTATE'});    
 uimenu(menu,'label','Save pts in a text file', ...
			'accelerator', 'I', ...
			'callback',{@SLURP,'SAVEPTS'});         
        
uimenu(menu,'label','Load anchors', ...
			'accelerator', 'C', ...
			'callback',{@SLURP,'LOAD'});        

uimenu(menu, 'label','Load segmentation','accelerator', 'L', 'callback',...
       {@SLURP, 'LOADSTATE'});

uimenu(menu, 'label', 'View energy map', 'accelerator', 'E', ...
       'callback', {@SLURP, 'ENERGYMAP'});

uimenu(menu,'label','Correct frames', ...
			'accelerator', 'C', ...
			'callback',{@SLURP,'CORRECTASEDGETRACK'});          
        
uimenu(menu,'label','Reverse Tracking', ...
			'accelerator', 'T', ...
			'callback',{@SLURP,'REV'}); 	

uimenu(menu,'label','Undo correction', ...
			'accelerator', 'Z', ...
			'callback',{@SLURP,'UNDO'});    
%tic;
        
%===============================================================================
% UPDATECONTOUR  - update displayed contour

function state = UpdateContour(state)

% reset (following clear by double-click)
if ~ishandle(state.CLH),
	state.ANCHORS(1:end-1,:) = [];
	state.ALH(1:end-1) = [];
	state.XY(:,:,state.CURFRAME) = 0;
	state.CLH = [];
end;

% purge anchor points following deletion
k = find(~ishandle(state.ALH));
state.ANCHORS(k,:) = [];
state.ALH(k) = [];

if state.ANCHORS_VISIBLE == 0
    state.ALH = [];
end

% create contour line
if isempty(state.CLH),
	state.CLH = line(state.ANCHORS(1,1),state.ANCHORS(1,2), ...
                         'color','y','LineWidth', 1, 'tag','CONTOUR','hitTest','off');
end;

if isempty(state.FWD_CLH),
    state.FWD_CLH = line(state.ANCHORS(1,1),state.ANCHORS(1,2), ...
                         'color','b','LineWidth', 1, 'tag','FWD_CONTOUR','hitTest','off');
end

if isempty(state.BWD_CLH),
    state.BWD_CLH = line(state.ANCHORS(1,1),state.ANCHORS(1,2),...
                         'color','g','LineWidth',1,'tag', ...
                         'BWD_CONTOUR','hitTest', 'off');
end

% update contour line from anchor points
switch size(state.ANCHORS,1),
	case 1,
		%state.XY(:,:,state.CURFRAME) = state.ANCHORS;
		set(state.CLH,'xdata',state.ANCHORS(1,1),'ydata',state.ANCHORS(1,2));
	case 2,
		k = [0 ; sqrt(sum(diff(state.ANCHORS).^2,2))];
		state.XY(:,:,state.CURFRAME) = interp1(k,state.ANCHORS,linspace(0,k(end),state.NPOINTS),'linear');
		set(state.CLH,'xdata',state.XY(:,1,state.CURFRAME),'ydata',state.XY(:,2,state.CURFRAME));
	otherwise,
		k = [0 ; cumsum(sqrt(sum(diff(state.ANCHORS).^2,2)))];
		state.XY(:,:,state.CURFRAME) = interp1(k,state.ANCHORS,linspace(0,k(end),state.NPOINTS),'pchip');
		set(state.CLH,'xdata',state.XY(:,1, ...
                                               state.CURFRAME), ...
                              'ydata',state.XY(:,2, ...
                                               state.CURFRAME));
                
                
                if ~isempty(state.FWD)
                    set(state.FWD_CLH, 'xdata', state.FWD(:,1,state.CURFRAME),...
                                      'ydata', state.FWD(:,2, ...
                                                         state.CURFRAME));
                end
                
                if ~isempty(state.BWD)
                    set(state.BWD_CLH, 'xdata', state.BWD(:,1,state.CURFRAME),...
                                      'ydata', state.BWD(:,2, ...
                                                         state.CURFRAME));
                end
end;
drawnow;



