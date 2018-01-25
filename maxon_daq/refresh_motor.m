function refresh_motor(handles)
global control
% To send command to motor controller with given parametres stored in handles
motorEnabled	= control.motorEnabled;
digi_2_value	= control.digi_2_value;
rpm				= control.rpm;
currentlimit	= control.currentlimit;

rpmV			= rpm / 4000;
currentlimitV	= currentlimit/1.2;

outputSingleScan(handles.sessionAnalOut, [rpmV, currentlimitV])
outputSingleScan(handles.sessionDigi, [motorEnabled, digi_2_value])
display('Control Command Sent')