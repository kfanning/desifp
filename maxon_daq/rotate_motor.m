function rotate_motor(handles)
global control
	% motion parametres
	gearRatio		= control.gearRatio;
	digi_2_value	= control.digi_2_value;
	rpm				= control.rpm;
	currentlimit	= control.currentlimit;
	% convert to analogue voltages
	rpmV			= rpm / 4000;
	currentlimitV	= currentlimit/1.2;
	% calculate samples for timed motion control
	rate			= handles.sessionAnalOut.Rate;		% in samples/second
	angle			= control.angle;					% in radians
	% before motion, set Speed and Current Limit
	analCommand = [rpmV, currentlimitV];
	outputSingleScan(handles.sessionAnalOut, analCommand)
	switch control.angleControlMode
		case 'Sampled'
			runtime			= angle / rpm2rad(rpm/gearRatio);	% in seconds
			nSamples		= round(rate * runtime);	% number of samples needed
			% generate command sequence
			analCommandSequence = repmat(analCommand, nSamples-1, 1);
			outputSingleScan(handles.sessionAnalOut, analCommand)
			% Start rotation
			queueOutputData(handles.sessionAnalOut, analCommandSequence)
			display('Command Sequence Queued')
			outputSingleScan(handles.sessionDigi, [1, digi_2_value])
			display('Command Sequence Initiated')
			startForeground(handles.sessionAnalOut);
			% Stop rotation
			display('Rotating, wait...')
			outputSingleScan(handles.sessionDigi, [0, digi_2_value])
			display('Command Sequence Halted')
		case 'Timed'
			runtime	= angle / rpm2rad(rpm/gearRatio);	% in seconds
			% before motion, set Speed and Current Limit
			% Start rotation
			outputSingleScan(handles.sessionDigi, [1, digi_2_value])
			display('Waiting for Motion to Complete')
			pause(runtime-1/rate)			
			% Stop roration
			outputSingleScan(handles.sessionDigi, [0, digi_2_value])
			display('Motion Complete ')
			display(['Total Runtime ' num2str(runtime) ' seconds'])
		case 'Integrated'
			control.rotationIntegratedControl = 1;
			control.angleIntegrated = 0;
			control.angleExcess = NaN;
	end
end