%This function is called when ths user clicked on the energy map
% It used to refresh the current frame and to save the correction's
% parameters


function txt = myupdatefcn(empty, event_obj)


% Get cursor's position
pos = get(event_obj, 'Position');

%Get figure intensity
energy = get(event_obj.Target);
intensity = energy.CData;
%Get figure's axes
hAxes = get(get(event_obj,'Target'),'Parent');
%Get the attached figure and the associated state 
figure = getappdata(hAxes, 'fig');
state = get(figure,'userData');

f= pos(1);
%Set main figure's property
set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
set(state.TH,'string',sprintf('%04d',f));
 set(state.CLH,'xdata',state.XY(:,1,f),...
                              'ydata',state.XY(:,2,f));
 %set(state.ASM_CLH, 'xdata', state.ASM(:,1,f),'ydata',state.ASM(:,2:,f));
 
state.CURFRAME = f;                           

 % If the window for correction parameters is open
 if ~isempty(findobj('Type', 'figure', 'name', 'CorrecParameters'));
     correcParameters = findobj('style', 'edit', 'tag', 'frame');
     data1 = get(correcParameters(1), 'UserData');
     selected1 = data1.selected;
     if(selected1);  
         set(correcParameters(1), 'String', num2str(f));            
     else
         set(correcParameters(2), 'String', num2str(f));                
     end
     frame1 = str2double(get(correcParameters(1), 'String'));
     frame2 = str2double(get(correcParameters(2), 'String'));
     state.PARAM_CORREC.FBEF = min(frame1, frame2);
     state.PARAM_CORREC.FAFT = max(frame1, frame2);
 end


 if ~isempty(findobj('Type', 'figure', 'name', 'ReverseTracking')) ||...
         ~isempty(findobj('Type', 'figure', 'name', 'Reverse PF Tracking')) ...
     correcParameters = findobj('style', 'edit', 'tag', 'trackFrame');
     data1 = get(correcParameters(1), 'UserData');
     selected1 = data1.selected;
     if(selected1);  
         set(correcParameters(1), 'String', num2str(f));            
     else
         set(correcParameters(2), 'String', num2str(f));                
     end
     frame1 = str2double(get(correcParameters(1), 'String'));
     frame2 = str2double(get(correcParameters(2), 'String'));
     state.PARAM_TRACK.FBEF = max(frame1, frame2);
     state.PARAM_TRACK.FAFT = min(frame1, frame2);
     
     
 end
 
 
set(figure,'userData',state);




txt = {['Frame: ',num2str(pos(1))],...
       ['Snaxel: ',num2str(pos(2))], ...
       ['Value:', num2str(intensity(pos(2),pos(1)))]};
   
   
