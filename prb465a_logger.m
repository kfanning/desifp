function varargout = prb465a_logger(varargin)
% PRB465A_LOGGER MATLAB code for prb465a_logger.fig
%      PRB465A_LOGGER, by itself, creates a new PRB465A_LOGGER or raises the existing
%      singleton*.
%
%      H = PRB465A_LOGGER returns the handle to a new PRB465A_LOGGER or the handle to
%      the existing singleton*.
%
%      PRB465A_LOGGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRB465A_LOGGER.M with the given input arguments.
%
%      PRB465A_LOGGER('Property','Value',...) creates a new PRB465A_LOGGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before prb465a_logger_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to prb465a_logger_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help prb465a_logger

% Last Modified by GUIDE v2.5 29-Jul-2015 13:22:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prb465a_logger_OpeningFcn, ...
                   'gui_OutputFcn',  @prb465a_logger_OutputFcn, ...
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


% --- Executes just before prb465a_logger is made visible.
function prb465a_logger_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prb465a_logger (see VARARGIN)
[handles.FileName,handles.PathName] = uigetfile('*.mat','Select Particle Count Data');
load(fullfile(handles.PathName, handles.FileName));
handles.particleCount = particleCount;
handles.RecordIndex = length(particleCount)+1;
handles.fanSpeedSettings.SelectedObject=[];
handles.ContaminationSourcesList = {'Ahlen', 'Duan', 'Cassese', handles.fieldOtherSources.String};
handles.OtherSourcesExistence = 0;
% Choose default command line output for prb465a_logger
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prb465a_logger wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prb465a_logger_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function fieldDate_Callback(hObject, eventdata, handles)
% hObject    handle to fieldDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of fieldDate as text
%        str2double(get(hObject,'String')) returns contents of fieldDate as a double
handles.particleCount(handles.RecordIndex).Date = int32(str2double(get(hObject,'String')));
handles.ChangeMade = 1;
guidata(hObject, handles);

% % --- Executes on selection change in fanSpeedSettings.
% function fanSpeedSettings_Callback(hObject, eventdata, handles)
% % hObject    handle to fanSpeedSettings (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: contents = cellstr(get(hObject,'String')) returns fanSpeedSettings contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from fanSpeedSettings
% handles.particleCount(handles.RecordIndex).FanSpeedSettings = handles.fanSpeedSettings.SelectedObject.String;
% guidata(hObject, handles);

% --- Executes when selected object is changed in fanSpeedSettings.
function fanSpeedSettings_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in fanSpeedSettings 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleCount(handles.RecordIndex).FanSpeedSettings = handles.fanSpeedSettings.SelectedObject.String;
handles.ChangeMade = 1;
guidata(hObject, handles);

% --- Executes on button press in checkAhlen.
function checkAhlen_Callback(hObject, eventdata, handles)
% hObject    handle to checkAhlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkAhlen
handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
handles.ChangeMade = 1;
guidata(hObject, handles);

% --- Executes on button press in checkDuan.
function checkDuan_Callback(hObject, eventdata, handles)
% hObject    handle to checkDuan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkDuan
handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
handles.ChangeMade = 1;
guidata(hObject, handles);

% --- Executes on button press in checkCassese.
function checkCassese_Callback(hObject, eventdata, handles)
% hObject    handle to checkCassese (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkCassese
handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
handles.ChangeMade = 1;
guidata(hObject, handles);

function fieldOtherSources_Callback(hObject, eventdata, handles)
% hObject    handle to fieldOtherSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of fieldOtherSources as text
%        str2double(get(hObject,'String')) returns contents of fieldOtherSources as a double
handles.OtherSourcesExistence=1;
handles.ContaminationSourcesList = {'Ahlen', 'Duan', 'Cassese', handles.fieldOtherSources.String};
handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
handles.ChangeMade = 1;
guidata(hObject, handles);

function fieldNumberOfPeople_Callback(hObject, eventdata, handles)
% hObject    handle to fieldNumberOfPeople (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of fieldNumberOfPeople as text
%        str2double(get(hObject,'String')) returns contents of fieldNumberOfPeople as a double
handles.particleCount(handles.RecordIndex).NumberOfPeople = ...
						int32(str2double(handles.fieldNumberOfPeople.String));
handles.ChangeMade = 1;
guidata(hObject, handles);

% --- Executes on button press in buttonSave.
function buttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.ChangeMade
	case 0
		display('No Change Made')
	case 1
		handles.particleCount(handles.RecordIndex).Date = ...
								int32(str2double(handles.fieldDate.String));
		if ~isempty(handles.fanSpeedSettings.SelectedObject)
			handles.particleCount(handles.RecordIndex).FanSpeedSettings = ...
								handles.fanSpeedSettings.SelectedObject.String;
		else
			handles.particleCount(handles.RecordIndex).FanSpeedSettings = [];
		end
		handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
		handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
		if ~isempty(handles.fieldNumberOfPeople.String)
			handles.particleCount(handles.RecordIndex).NumberOfPeople = ...
						int32(str2double(handles.fieldNumberOfPeople.String));
		else
			handles.particleCount(handles.RecordIndex).NumberOfPeople = [];
		end
		datatable=NaN(size(handles.table.Data));
		for i=1:numel(handles.table.Data)
			if ~isempty(handles.table.Data{i})
				datatable(i) = handles.table.Data{i};
			end
		end
		handles.particleCount(handles.RecordIndex).Cumulative05	= datatable(:,1);
		handles.particleCount(handles.RecordIndex).Cumulative10	= datatable(:,2);
		handles.particleCount(handles.RecordIndex).Cumulative30	= datatable(:,3);
		handles.particleCount(handles.RecordIndex).Cumulative50	= datatable(:,4);
		handles.particleCount(handles.RecordIndex).Cumulative100= datatable(:,5);
		handles.particleCount(handles.RecordIndex).Cumulative05Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative05);
		handles.particleCount(handles.RecordIndex).Cumulative10Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative10);
		handles.particleCount(handles.RecordIndex).Cumulative30Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative30);
		handles.particleCount(handles.RecordIndex).Cumulative50Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative50);
		handles.particleCount(handles.RecordIndex).Cumulative100Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative100);

		particleCount = handles.particleCount;
		save(fullfile(handles.PathName, handles.FileName),'particleCount')
		handles.RecordIndex = length(handles.particleCount) + 1;
		handles.ChangeMade = 0;
		guidata(hObject, handles);
end

% --- Executes on button press in buttonSaveAs.
function buttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to buttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.ChangeMade
	case 0
		display('No Change Made')
	case 1
		handles.particleCount(handles.RecordIndex).Date = ...
								int32(str2double(handles.fieldDate.String));
		if ~isempty(handles.fanSpeedSettings.SelectedObject)
			handles.particleCount(handles.RecordIndex).FanSpeedSettings = ...
								handles.fanSpeedSettings.SelectedObject.String;
		else
			handles.particleCount(handles.RecordIndex).FanSpeedSettings = [];
		end
		handles.SourcesLogical = logical([	handles.checkAhlen.Value,...
									handles.checkDuan.Value,...
									handles.checkCassese.Value,...
									handles.OtherSourcesExistence]);
		handles.particleCount(handles.RecordIndex).ContaminationSources = ...
					handles.ContaminationSourcesList(handles.SourcesLogical);
		if ~isempty(handles.fieldNumberOfPeople.String)
			handles.particleCount(handles.RecordIndex).NumberOfPeople = ...
						int32(str2double(handles.fieldNumberOfPeople.String));
		else
			handles.particleCount(handles.RecordIndex).NumberOfPeople = [];
		end
		datatable=NaN(size(handles.table.Data));
		for i=1:numel(handles.table.Data)
			if ~isempty(handles.table.Data{i})
				datatable(i) = handles.table.Data{i};
			end
		end
		handles.particleCount(handles.RecordIndex).Cumulative05		= datatable(:,1);
		handles.particleCount(handles.RecordIndex).Cumulative10		= datatable(:,2);
		handles.particleCount(handles.RecordIndex).Cumulative30		= datatable(:,3);
		handles.particleCount(handles.RecordIndex).Cumulative50		= datatable(:,4);
		handles.particleCount(handles.RecordIndex).Cumulative100	= datatable(:,5);
		handles.particleCount(handles.RecordIndex).Cumulative05Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative05);
		handles.particleCount(handles.RecordIndex).Cumulative10Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative10);
		handles.particleCount(handles.RecordIndex).Cumulative30Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative30);
		handles.particleCount(handles.RecordIndex).Cumulative50Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative50);
		handles.particleCount(handles.RecordIndex).Cumulative100Avg	=...
						mean(handles.particleCount(handles.RecordIndex).Cumulative100);
		particleCount = handles.particleCount;
		uisave('particleCount',handles.FileName);
		handles.RecordIndex = length(handles.particleCount)+1;
		handles.ChangeMade = 0;
		guidata(hObject, handles);
end

% --- Executes on button press in buttonViewDataTable.
function buttonViewDataTable_Callback(hObject, eventdata, handles)
% hObject    handle to buttonViewDataTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global particleCount
command = fullfile(handles.PathName, handles.FileName);
dos(command)
% load(fullfile(handles.PathName, handles.FileName))
openvar('particleCount')


% --- Executes when entered data in editable cell(s) in table.
function table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.ChangeMade = 1;
guidata(hObject, handles)


% --- Executes on button press in buttonClearTable.
function buttonClearTable_Callback(hObject, eventdata, handles)
% hObject    handle to buttonClearTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.table.Data=cell(10,5);