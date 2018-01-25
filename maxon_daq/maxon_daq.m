%% Initialise Session

devices = daq.getDevices;				% scan DAQ devices
session = daq.createSession('ni');		% create session

%% Create Channels

% Digital Output:	P0.0 - J5 1 Digital In 1 - Enable
Digi_1 = addDigitalChannel(session, 'dev1', 'Port0/Line0', 'OutputOnly');
% Digital Output:	P0.1 - J5 2 Digital In 2 - 
Digi_2 = addDigitalChannel(session, 'dev1', 'Port0/Line1', 'OutputOnly');
% Digital I/O:		P0.2 - J5 3 Digital I/O 1
Digi_3 = addDigitalChannel(session, 'dev1', 'Port0/Line2', 'InputOnly');
% Digital I/O:		P0.3 - J5 4 Digital I/O 2
Digi_4 = addDigitalChannel(session, 'dev1', 'Port0/Line3', 'InputOnly');

% Analogue Output:	AO.0 - J6 1 Analogue In 1+
AnalOut_1 = addAnalogOutputChannel(session, 'dev1', 'ao0', 'Voltage');
% Analogue Output:	AO.1 - J6 3 Analogue In 2+
AnalOut_2 = addAnalogOutputChannel(session, 'dev1', 'ao1', 'Voltage');
% Analogue Input:	AI.0 - J6 5 Analogue Out 1
AnalIn_1 = addAnalogInputChannel(session, 'dev1', 'ai0', 'Voltage');
% Analogue Input:	AI.1 - J6 6 Analogue Out 2
AnalIn_2 = addAnalogInputChannel(session, 'dev1', 'ai1', 'Voltage');

%% Output Control Commands

% Direction and rpm: -10 to +10 V for -40k (CW) to +40k (CCW) rpm
Digi_2_value = 0;
RPM = 40000;
CurrentLimit = 12; % amp
RPM_V = RPM / 4000;
CurrentLimit_V = CurrentLimit/3;

% enable sequence
enable = 1;
outputSingleScan(session, [enable, Digi_2_value, RPM_V, CurrentLimit_V])

% disable sequence
enable = 0;
outputSingleScan(session, [enable, Digi_2_value, RPM_V, CurrentLimit_V])

%% Input Monitoring

