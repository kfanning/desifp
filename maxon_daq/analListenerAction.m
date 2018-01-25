function analListenerAction(source, event)
% 			plot_actual_speed(event, handles);
%			handles.ylimit = handles_new.ylimit;
global control
	if control.rotationIntegratedControl == 1
		rotate_motor_integrated(event,handles);
	end
end