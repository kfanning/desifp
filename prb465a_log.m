%% PRB465A clean room log
% load data
load('particleCount.mat')
%% Prompt for new data
i=1;


particleCount(i).Cumulative05 = A(:,1);
particleCount(i).Cumulative10 = A(:,2);
particleCount(i).Cumulative30 = A(:,3);
particleCount(i).Cumulative50 = A(:,4);
particleCount(i).Cumulative100 = A(:,5);
particleCount(i).Cumulative05Avg	=...
				mean(particleCount(i).Cumulative05);
particleCount(i).Cumulative10Avg	=...
				mean(particleCount(i).Cumulative10);
particleCount(i).Cumulative30Avg	=...
				mean(particleCount(i).Cumulative30);
particleCount(i).Cumulative50Avg	=...
				mean(particleCount(i).Cumulative50);
particleCount(i).Cumulative100Avg	=...
				mean(particleCount(i).Cumulative100);
i=i+1;
%% Store new data
particleCount=orderfields(particleCount, {'Date', 'FanSpeedSettings', ...
	'ContaminationSources','NumberOfPeople','Cumulative05','Cumulative10','Cumulative30','Cumulative50','Cumulative100'})
save('particle_count_data.mat','particleCount')
%% Display confirmation