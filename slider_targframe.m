function slider_targframe(hobject, event)
%t = toc;
    
    %Get the main Figure (saved as User Data of the slider)
    currentFigure = get(hobject, 'userData'); 
    state = get(currentFigure, 'UserData');
    
    %state.USER_TIME = state.USER_TIME + t;
    
    % Get the slider's vaue
    f = get(hobject, 'Value');
    
    %Round the value
    f = round(f);
    
    %Set main figure's property
    img = im2double(uint8(mean(read(state.MH,f),3)));
% $$$     img = phasesymmono(img, ...
% $$$                        6, 3, 2.1, 0.55, 2.0, 1, -2).*state.MASK;
% $$$     img = (img-min(img(:)))./range(img(:));
% $$$                     
% $$$     thresh = graythresh(img);
% $$$     bwimg = im2bw(img, thresh);                                        
% $$$    
% $$$     s = strel('square', 5);
% $$$     bwimg = imdilate(bwimg, s);
% $$$     %    bwimg = imerode(bwimg,s);
% $$$     
% $$$     bwimg1=bwpropfilt(bwimg,'area',20);
% $$$     bwimg2=bwpropfilt(bwimg,'eccentricity',20);
% $$$     bwimg = bwimg1 & bwimg2;
% $$$     %    bwimg=bwpropfilt(bwimg,'area',5);
% $$$     stats = regionprops(bwimg,'all')
%img = im2uint8(bwimg);
    img = im2uint8(img);
    %set(state.IH,'cdata',uint8(mean(read(state.MH,f),3)));
    set(state.IH,'cdata',img);
    set(state.TH,'string',sprintf('%04d',f));
    set(state.CLH,'xdata',state.XY(:,1,f),...
                              'ydata',state.XY(:,2,f));
    state.CURFRAME = f;     
    state.TARGFRAME = f;
    set(currentFigure,'userData',state);
    
    tic
end