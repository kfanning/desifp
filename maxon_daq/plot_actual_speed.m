function plot_actual_speed(event, handles)
global control
	newlimit = max(abs(event.Data))*control.rpmVscale*1.13;
	if mod(event.TimeStamps(1), control.PlotInterval) == 0
		% reset plot, restart from left
		hold(handles.axesSpeed,'off')
		control.ylimit = newlimit;
		plot(handles.axesSpeed, ...
				event.TimeStamps, event.Data*control.rpmVscale, '-r')
		xlim(handles.axesSpeed,	...
			[event.TimeStamps(1), event.TimeStamps(1)+control.PlotInterval])
		ylim(handles.axesSpeed, ...
			[-control.ylimit, control.ylimit])
		xlabel(handles.axesSpeed, ...
			'Time/s', 'FontSize', 7.5)
		ylabel(handles.axesSpeed, ...
			'Speed/rpm', 'FontSize', 7.5)
% 				title(handles.axesSpeed, ...
% 					'Actual Speed Monitoring from Controller Output')
% 				legend(handles.axesSpeed, ...
% 					'Actual Speed at 20kHz')
	else
		% add data to graph
		if newlimit>control.ylimit
			control.ylimit = newlimit;
		end
		hold(handles.axesSpeed,'on')
		plot(handles.axesSpeed,...
			event.TimeStamps, event.Data*control.rpmVscale, '-r')
		ylim(handles.axesSpeed, ...
			[-control.ylimit, control.ylimit])
	end
end