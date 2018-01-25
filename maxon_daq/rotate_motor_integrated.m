function rotate_motor_integrated(event, handles)
% preserve handle values by default
	global control repeatIndex
	digi_2_value	= control.digi_2_value;
	switch control.rotating
		case 0
			% motor is stationary; start integration-controlled rotation
			disp('Motor Stationary. Starting Rotation')
			rpm				= control.rpm;
			currentlimit	= control.currentlimit;
			% convert to analogue voltages
			rpmV			= rpm / 4000;
			currentlimitV	= currentlimit/1.2;
			% before motion, set Speed and Current Limit
			outputSingleScan(handles.sessionAnalOut, [rpmV, currentlimitV])
			control.motorEnabled = 1;
			control.rotating = 1;
			outputSingleScan(handles.sessionDigi, [1, digi_2_value])
		case 1
% 			disp('Still Rotating')
			speeddata = rpm2rad(event.Data * control.rpmVscale) / control.gearRatio;
			angleIntegrated = trapz(event.TimeStamps, speeddata);
			control.angleIntegrated = control.angleIntegrated + angleIntegrated;
			if control.angleIntegrated >= control.angle
				disp('Reached Desired Angle')
				switch control.motorEnabled
					case 1
						% motor still enabled; disable it
						% will keep going due to moment of inertia
						disp('Disabling Motor')
						control.motorEnabled = 0;
						outputSingleScan(handles.sessionDigi, [0, digi_2_value])
					case 0
						% motor already disabled
						disp('Motor State is Disabled')
						if min(event.Data * control.rpmVscale) > 0
							% still moving at greater than 10 rpm
							% do nothing, keep integrating and wait for slowdown
							display('Momentum Still Going, Wait...')
						else
							% motor had just stopped rotating
							% finalise integrated control and clear tags to quit
							control.rotating = 0;
							control.rotationIntegratedControl = 0;
							control.angleExcess = control.angleIntegrated - control.angle;
							disp('Full Stop Reached. Rotation Completed.')
							disp(['Defined Angle:     ', num2str(control.angle)])
							disp(['Actual Angle:      ', num2str(control.angleIntegrated)])
							disp(['Excess Angle/rad:  ', num2str(control.angleExcess)])
							disp(['Excess Angle/deg:  ', num2str(control.angleExcess/2/pi*360)])
							if exist('repeatIndex','var')
								control.excess(repeatIndex) = control.angleExcess/2/pi*360;
							end
						end
				end
			end
	end
end