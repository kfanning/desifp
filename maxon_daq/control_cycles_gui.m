% ------------------- GUI Initialisation --------------------------------

function varargout = control_cycles_gui(varargin)
	% CONTROL_CYCLES_GUI MATLAB code for control_cycles_gui.fig
	%      CONTROL_CYCLES_GUI, by itself, creates a new CONTROL_CYCLES_GUI or raises the existing
	%      singleton*.
	%
	%      H = CONTROL_CYCLES_GUI returns the handle to a new CONTROL_CYCLES_GUI or the handle to
	%      the existing singleton*.
	%
	%      CONTROL_CYCLES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
	%      function named CALLBACK in CONTROL_CYCLES_GUI.M with the given input arguments.
	%
	%      CONTROL_CYCLES_GUI('Property','Value',...) creates a new CONTROL_CYCLES_GUI or raises the
	%      existing singleton*.  Starting from the left, property value pairs are
	%      applied to the GUI before control_cycles_gui_OpeningFcn gets called.  An
	%      unrecognized property name or invalid value makes property application
	%      stop.  All inputs are passed to control_cycles_gui_OpeningFcn via varargin.
	%
	%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
	%      instance to run (singleton)".
	%
	% See also: GUIDE, GUIDATA, GUIHANDLES

	% Edit the above text to modify the response to help control_cycles_gui

	% Last Modified by GUIDE v2.5 15-Jul-2015 21:06:30

	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
					   'gui_Singleton',  gui_Singleton, ...
					   'gui_OpeningFcn', @control_cycles_gui_OpeningFcn, ...
					   'gui_OutputFcn',  @control_cycles_gui_OutputFcn, ...
					   'gui_LayoutFcn',  [] , ...
					   'gui_Callback',   []);
	if nargin && ischar(varargin{1})
		gui_State.gui_Callback = str2func(varargin{1});
	end

	if nargout
		[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
	else
		gui_mainfcn(gui_State, varargin{:});
	end
	% End initialization code - DO NOT EDIT
end

% --- Executes just before control_cycles_gui is made visible.
function control_cycles_gui_OpeningFcn(hObject, eventdata, handles, varargin)
	% This function has no output args, see OutputFcn.
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	% varargin   command line arguments to control_cycles_gui (see VARARGIN)

	%% Initialise Session
	global control					% everything non-GUI related
	control.devices = daq.getDevices;			% scan DAQ devices
	control.digi_2_value = 0;
	control.gearRatio = 1296/1;
	control.motorEnabled = 0;
	control.rotationIntegratedControl = 0;
	control.rotating = 0;
	control.angleIntegrated = 0;
	control.rpmVscale = 40000/4;
	%% Create Channels

	handles.sessionDigi = daq.createSession('ni');
	% Digital Output:	P0.0 - J5 1 Digital In 1 - rpmField
	digi_1 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line0', 'OutputOnly');
	% Digital Output:	P0.1 - J5 2 Digital In 2 - 
	digi_2 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line1', 'OutputOnly');
	% Digital I/O:		P0.2 - J5 3 Digital I/O 1
	digi_3 = addDigitalChannel(handles.sessionDigi, 'dev1', 'Port0/Line2', 'InputOnly');
	% Digital I/O:		P0.3 - J5 4 Digital I/O 2
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
								@analListenerAction);	
		function analListenerAction(source, event)
			plot_actual_speed(event, handles);
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

	control.rpm = str2double(get(handles.rpmField, 'String'));
	control.currentlimit = str2double(get(handles.currentlimitField, 'String'));
	control.angle = str2double(get(handles.radField, 'String'));
	control.angleControlMode = handles.modeBG.SelectedObject.String;
	display(['Initial RPM: ', num2str(control.rpm)])
 	startBackground(handles.sessionAnalIn_0)

	% Choose default command line output for control_cycles_gui
	handles.output = hObject;

	% Update handles structure
	guidata(hObject, handles);
end
% UIWAIT makes control_cycles_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = control_cycles_gui_OutputFcn(hObject, eventdata, handles) 
	% varargout  cell array for returning output args (see VARARGOUT);
	% hObject    handle to figure
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)

	% Get default command line output from handles structure
	varargout{1} = handles.output;
end
% ----------------------- Control Panel -------------------------------------

function rpmField_Callback(hObject, eventdata, handles)
	rpmInput = str2double(get(handles.rpmField, 'String'));
	control_cycles_gui_DataSyncHost(hObject, 'rpm', rpmInput, handles)
% 	handles.rpm = handles_new.rpm;
% 	guidata(hObject, handles)
end

% --- Executes on slider movement.
function rpmSlider_Callback(hObject, eventdata, handles)
	rpmInput = round(get(handles.rpmSlider, 'Value'), -2, 'Decimal');
	control_cycles_gui_DataSyncHost(hObject, 'rpm', rpmInput, handles)
% 	handles.rpm = handles_new.rpm;
% 	guidata(hObject, handles)
end

function currentlimitField_Callback(hObject, eventdata, handles)
	currentlimitInput = str2double(get(handles.currentlimitField, 'String'));
	control_cycles_gui_DataSyncHost(hObject, 'currentlimit', currentlimitInput, handles)
% 	handles.currentlimit = handles_new.currentlimit;
% 	guidata(hObject, handles)
end
% --- Executes on slider movement.
function currentlimitSlider_Callback(hObject, eventdata, handles)
	currentlimitInput = round(get(handles.currentlimitSlider, 'Value'), 1, 'Decimal');
	control_cycles_gui_DataSyncHost(hObject, 'currentlimit', currentlimitInput, handles)
% 	handles.currentlimit = handles_new.currentlimit;
% 	guidata(hObject, handles)
end

function degField_Callback(hObject, eventdata, handles)
	degInput = str2double(get(handles.degField,'String'));
	control_cycles_gui_DataSyncHost(hObject, 'deg', degInput, handles)
% 	handles.angle = handles_new.angle;
% 	guidata(hObject, handles)
end

function radField_Callback(hObject, eventdata, handles)
	radInput = str2double(get(handles.radField,'String'));
	control_cycles_gui_DataSyncHost(hObject, 'rad', radInput, handles)
% 	handles.angle = handles_new.angle;
% 	guidata(hObject, handles)
end

function revField_Callback(hObject, eventdata, handles)
	revInput = str2double(get(handles.revField,'String'));
	control_cycles_gui_DataSyncHost(hObject, 'rev', revInput, handles)
% 	handles.angle = handles_new.angle;
% 	guidata(hObject, handles)
end
% --- Executes when selected object is changed in radiobutton1.
function modeBG_SelectionChangedFcn(hObject, eventdata, handles)
global control
	% hObject    handle to the selected object in radiobutton1 
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	control.angleControlMode = handles.modeBG.SelectedObject.String;
	display(['Angle Control Mode is set to: ', control.angleControlMode])
end

function enableButton_Callback(hObject, eventdata, handles)
global control
	control.motorEnabled = get(handles.enableButton,'Value');
	refresh_motor(handles)
end

% --- Executes on button press in startButton.
function startButton_Callback(hObject, eventdata, handles)
% hObject    handle to startButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	rotate_motor(handles)
% 	guidata(hObject, handles)
end


% --- Executes on button press in terminateListener.
function terminateListener_Callback(hObject, eventdata, handles)
	% hObject    handle to terminateListener (see GCBO)
	% eventdata  reserved - to be defined in a future version of MATLAB
	% handles    structure with handles and user data (see GUIDATA)
	delete(handles.analListener)
	% Hint: get(hObject,'Value') returns toggle state of terminateListener
end
