function acquireActualSpeed(source, event)
	global speedData
	
	if mod(event.TimeStamps(1), speedData.PlotInterval) == 0
		% reset speeddata, erase previous data
		speedData.SpeedT = NaN(speedData.Rate*speedData.PlotInterval,1);
		speedData.SpeedV = NaN(speedData.Rate*speedData.PlotInterval,1);
	end
	% Indices for Data Batch when Notified
	% size of batch is defined by NotifyWhenDataAvailableExceeds
	startIndex = uint64(	...
						mod(event.TimeStamps(1),speedData.PlotInterval)...
						*speedData.Rate...
						)...
				 +1;
	endIndex = startIndex + speedData.NotifyWhenDataAvailableExceeds -1;
	speedData.SpeedT(startIndex:endIndex) = event.TimeStamps;
	speedData.SpeedV(startIndex:endIndex) = event.Data;
	
	display(event.TimeStamps(1))