
% This function is called when the user select an edit. 
% It used to know which edit is selected

function editfunction(hObject, eventdata, handles)
data = get(hObject,'UserData');
% This edit is selected
data.selected = 1;
% data.OtherEdit is the handle of the other Edit in the Parameters windows
otherEdit = data.otherEdit;
set(hObject, 'UserData', data);
dataOtherEdit = get(otherEdit, 'UserData');
% The other Edit is not selected
dataOtherEdit.selected = 0;
set(otherEdit, 'UserData', dataOtherEdit);
end

