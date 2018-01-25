function excess_data = reproducibility
	global control repeatIndex					% everything non-GUI related
	control.devices = daq.getDevices;			% scan DAQ devices
	control.digi_2_value = 0;
	control.gearRatio = 1296/1;
	control.motorEnabled = 0;
	control.rotationIntegratedControl = 0;
	control.rotating = 0;
	control.angleIntegrated = 0;
	control.rpmVscale = 40000/4;

	handles.sessionDigi = daq.createSession('ni');
	digi_1 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line0', 'OutputOnly');
	digi_2 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line1', 'OutputOnly');
	digi_3 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line2', 'InputOnly');
	digi_4 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line3', 'InputOnly');

	handles.sessionAnalOut = daq.createSession('ni');
	handles.sessionAnalOut.Rate = 5000;
	% Analogue Output:	AO.0 - J6 1 Analogue In 1+
	analOut_1 = addAnalogOutputChannel(handles.sessionAnalOut, 'dev1', 'ao0', 'Voltage');
	% Analogue Output:	AO.1 - J6 3 Analogue In 2+
	analOut_2 = addAnalogOutputChannel(handles.sessionAnalOut, 'dev1', 'ao1', 'Voltage');

	% Analogue Input:	AI.0 - J6 5 Analogue Out 1
	% Actual Speed at 20kHz
	handles.sessionAnalIn_0 = daq.createSession('ni');
	handles.sessionAnalIn_0.IsContinuous = 1;
	% 	handles.sessionAnalIn_0.DurationInSeconds = 10; % for debugging
	handles.sessionAnalIn_0.Rate = 20000;
	handles.sessionAnalIn_0.NotifyWhenDataAvailableExceeds = 2000;
	control.PlotInterval = 5;
	control.ylimit = 100;

	analIn_1 = addAnalogInputChannel(handles.sessionAnalIn_0, 'dev1', 'ai0', 'Voltage');
	handles.analListener = addlistener(handles.sessionAnalIn_0,'DataAvailable',...
								@analListenerAction1);	
		function analListenerAction1(source, event)
		% 			plot_actual_speed(event, handles);
		%			handles.ylimit = handles_new.ylimit;
			if control.rotationIntegratedControl == 1
				rotate_motor_integrated(event,handles);
			end
		end
	% Analogue Input:	AI.1 - J6 6 Analogue Out 2
	% Actual Speed Averaged at 10kHz
	% handles.sessionAnalIn_1 = daq.createSession('ni');
	% handles.sessionAnalIn_1.IsContinuous = 1;
	% handles.sessionAnalIn_1.Rate = 10000;
	% handles.sessionAnalIn_1.NotifyWhenDataAvailableExceeds = 1000;
	% analIn_2 = addAnalogInputChannel(handles.sessionAnalIn_1, 'dev1', 'ai1', 'Voltage');
	% 
	% analListener = addlistener(handles.sessionAnalIn,'DataAvailable', ...
	% 							@acquireActualSpeedAverage);
	% startBackground(handles.sessionAnalIn_1)

	control.rpm = 2500;
	control.currentlimit = 12;
	control.angle = pi;
	control.angleControlMode = 'Integrated';
	display(['Initial RPM: ', num2str(control.rpm)])
	startBackground(handles.sessionAnalIn_0)

	for repeatIndex = 1:3000
		rotate_motor(handles)
		pause(24)
	end
	excess_data = control.excess;
	save('excess_data', 'excess_data')
end