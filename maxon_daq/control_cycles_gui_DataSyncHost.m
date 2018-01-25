function control_cycles_gui_DataSyncHost(hObject, structSpecifier, valueInput, handles)
% prserve handles data by default
% handles.rpm = handles.rpm;
% handles.currentlimit = handles.currentlimit;
% handles.angle = handles.angle;
global control
	switch structSpecifier
		case 'rpm'
			if (-40000 <= valueInput) && (valueInput <= 40000)
				set(handles.rpmSlider,'Value', valueInput)
				set(handles.rpmField, 'String', num2str(valueInput))
				control.rpm = valueInput;
				display(['New RPM Accepted: ', num2str(control.rpm)])
				refresh_motor(handles)
			else
				display('Enter Valid RPM')
			end
		case 'currentlimit'
			if (0 <= valueInput) && (valueInput <= 12)		
				set(handles.currentlimitSlider,'Value', valueInput)
				set(handles.currentlimitField, 'String', num2str(valueInput))
				control.currentlimit = valueInput;
				display('New Current Limit Accepted')
				refresh_motor(handles)
			else
				display('Enter Valid Current Limit')
			end
		case 'deg'
			control.angle = valueInput/360*(2*pi);
			set(handles.degField, 'String', num2str(control.angle/(2*pi)*360))
			set(handles.radField, 'String', num2str(control.angle))
			set(handles.revField, 'String', num2str(control.angle/(2*pi)))
		case 'rad'
			control.angle = valueInput;
			set(handles.degField, 'String', num2str(control.angle/(2*pi)*360))
			set(handles.radField, 'String', num2str(control.angle))
			set(handles.revField, 'String', num2str(control.angle/(2*pi)))
		case 'rev'
			control.angle = valueInput*(2*pi);
			set(handles.degField, 'String', num2str(control.angle/(2*pi)*360))
			set(handles.radField, 'String', num2str(control.angle))
			set(handles.revField, 'String', num2str(control.angle/(2*pi)))
	end
end