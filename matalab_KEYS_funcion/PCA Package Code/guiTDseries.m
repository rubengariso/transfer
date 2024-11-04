function varargout = guiTDseries(varargin)
% GUITDSERIES MATLAB code for guiTDseries.fig
%      GUITDSERIES, by itself, creates a new GUITDSERIES or raises the existing
%      singleton*.
%
%      H = GUITDSERIES returns the handle to a new GUITDSERIES or the handle to
%      the existing singleton*.
%
%      GUITDSERIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUITDSERIES.M with the given input arguments.
%
%      GUITDSERIES('Property','Value',...) creates a new GUITDSERIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiTDseries_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiTDseries_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiTDseries

% Last Modified by GUIDE v2.5 08-Dec-2015 14:12:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiTDseries_OpeningFcn, ...
                   'gui_OutputFcn',  @guiTDseries_OutputFcn, ...
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


% --- Executes just before guiTDseries is made visible.
function guiTDseries_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiTDseries (see VARARGIN)

% Choose default command line output for guiTDseries
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiTDseries wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiTDseries_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_process.
function popupmenu_process_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_process contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_process
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
scenario=contents{get(hObject,'Value')};% returns selected item from popupmenu1


switch scenario
    case 'static'
        set(handles.text_mu_hObject,'Visible','Off'); 
        set(handles.edit_mu_hObject,'Visible','Off'); 
        set(handles.text_phi_hObject,'Visible','Off'); 
        set(handles.edit_phi_hObject,'Visible','Off'); 
        set(handles.text_freq_hObject,'Visible','Off'); 
        set(handles.edit_freq_hObject,'Visible','Off'); 
        set(handles.popupmenu_model_hObject,'Visible','Off');
        set(handles.uipanel_train_hObject,'Visible','On');
        set(handles.refSigma_hObject,'Visible','On');
    case {'AR','MA','ARI','IMA'}
        set(handles.text_mu_hObject,'Visible','On'); 
        set(handles.edit_mu_hObject,'Visible','On'); 
        set(handles.text_phi_hObject,'Visible','On'); 
        set(handles.edit_phi_hObject,'Visible','On'); 
        set(handles.text_freq_hObject,'Visible','Off'); 
        set(handles.edit_freq_hObject,'Visible','Off'); 
        set(handles.popupmenu_model_hObject,'Visible','Off');
        set(handles.uipanel_train_hObject,'Visible','On');
        set(handles.refSigma_hObject,'Visible','On');
    case 'NSS'
        set(handles.text_mu_hObject,'Visible','Off'); 
        set(handles.edit_mu_hObject,'Visible','Off'); 
        set(handles.text_phi_hObject,'Visible','Off'); 
        set(handles.edit_phi_hObject,'Visible','Off'); 
        set(handles.text_freq_hObject,'Visible','On'); 
        set(handles.edit_freq_hObject,'Visible','On');
        set(handles.popupmenu_model_hObject,'Visible','Off');
        set(handles.uipanel_train_hObject,'Visible','On');
        set(handles.refSigma_hObject,'Visible','On');
    case 'Saved model'
        set(handles.text_mu_hObject,'Visible','Off'); 
        set(handles.edit_mu_hObject,'Visible','Off'); 
        set(handles.text_phi_hObject,'Visible','Off'); 
        set(handles.edit_phi_hObject,'Visible','Off'); 
        set(handles.text_freq_hObject,'Visible','Off'); 
        set(handles.edit_freq_hObject,'Visible','Off');
        set(handles.popupmenu_model_hObject,'Visible','On'); 
        set(handles.uipanel_train_hObject,'Visible','Off');
        set(handles.refSigma_hObject,'Visible','Off');
        
        % update popupmenu
        vars=evalin('base','who');
        set(handles.popupmenu_model_hObject,'String',['-- Model --';vars])
        index=find(strcmp(handles.nameModel, vars))+1;
        if isempty(index)==1,
            index=1;
        end
        set(handles.popupmenu_model_hObject,'Value',index)
end

handles.scenarioProcess=scenario;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_process_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
scenario={'static'; 'AR'; 'MA';'ARI';'IMA';'NSS';'Saved model'};
set(hObject,'String',['-- Process --'; scenario])
handles.process_hObject=hObject;
guidata(hObject, handles);


function edit_phi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phi as text
%        str2double(get(hObject,'String')) returns contents of edit_phi as a double
handles.phi=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Visible','Off'); 
handles.edit_phi_hObject=hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function text_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.text_phi_hObject=hObject;
guidata(hObject, handles);


function edit_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_mu as a double
handles.mu=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Visible','Off'); 
handles.edit_mu_hObject=hObject;
guidata(hObject, handles);


function edit_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freq as text
%        str2double(get(hObject,'String')) returns contents of edit_freq as a double
handles.freq=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Visible','Off'); 
handles.edit_freq_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.text_freq_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.text_mu_hObject=hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status_hObject,'String','running')
pause (1e-3)

try
    
    % get process parameters
    if get(handles.process_hObject,'Value')==1,
       error('Chose type of process.'); 
    end
    if get(handles.popupmenu_model_hObject,'Value')==1 && strcmp(handles.scenarioProcess,'Saved model')==1,
       error('Load model.'); 
    end
    if get(handles.faultType_hObject,'Value')==1,
       error('Chose fault type.'); 
    end
    if get(handles.faultScenario_hObject,'Value')==1,
       error('Chose fault scenario.'); 
    end
    if strcmp(handles.scenarioProcess,'Saved model')==1,
        handles.model=evalin('base',handles.nameModel);
    end
    if get(handles.refSigma_hObject,'Value')==1,
       error('Chose reference variance.'); 
    end

    n0=5000;
    
    handles.phi=str2num(get(handles.edit_phi_hObject,'String'));
    handles.mu=str2num(get(handles.edit_mu_hObject,'String'));
    handles.freq=str2num(get(handles.edit_freq_hObject,'String'));
    
    handles.m=str2double(get(handles.edit_m_hObject,'String'));
    handles.n=str2double(get(handles.edit_n_hObject,'String'));
    handles.p=str2double(get(handles.edit_p_hObject,'String'));
    
    handles.n_noc=str2double(get(handles.edit_n_noc_hObject,'String'));
    handles.n_fault=str2double(get(handles.edit_n_fault_hObject,'String'));
    handles.d=str2double(get(handles.edit_d_hObject,'String'));
    
    switch handles.refSigma
        case 'Signal variance'
            fault.refSigma='signal';
        case 'Errors variance'
            fault.refSigma='error';
    end
    
    switch handles.scenarioProcess
        case 'static'
            T_params=[];
            A_params=[];
            scenarioProcess=handles.scenarioProcess;
        case {'AR','MA','ARI','IMA'}
            T_params.phi=handles.phi*ones(1,handles.p);
            T_params.mu=handles.mu*ones(1,handles.p);
            A_params=[];
            scenarioProcess=handles.scenarioProcess;
        case 'NSS'
            T_params=[];
            A_params=handles.freq*ones(1,handles.m-1);
            scenarioProcess=handles.scenarioProcess;
        case 'Saved model'
            T_params=handles.model.T_params;
            A_params=handles.model.A_params;
            scenarioProcess=handles.model.scenario;
            fault.refSigma=handles.model.refSigma;
            handles.n=0;
            handles.m=handles.model.m;
            handles.p=handles.model.p;   
    end
    
    fault.scenario=handles.faultScenario;
    fault.type=handles.faultType;
    
    fault.loc=n0+handles.n+handles.n_noc+1;%The the observation corresponding to the location where the fault begins.
    switch fault.scenario
        case 'step',
            fault.shift=handles.d;% The magnitude of the final shift in the step and ramp in standard deviations
        case 'ramp'
            fault.shift=handles.d;% The magnitude of the final shift in the step and ramp in standard deviations
            fault.speed=handles.n_fault;
        case 'scale'
            fault.rescale=handles.d;
    end
    
    % % train
    % output_train = TDseries(n0+handles.n,handles.m,handles.p,T_params,A_params,handles.scenarioProcess);
    % A_params.base=output_train.A_params.base;
    % T_params.base=output_train.T_params.base;
    
    % % test
    output = TDseries(n0+handles.n+handles.n_noc+handles.n_fault,handles.m,handles.p,T_params,A_params,scenarioProcess, fault);
    
    handles.Xtrain=output.X(n0+1:n0+handles.n,:);
    handles.Xtest=output.X(n0+handles.n+1:end,:);
    handles.Ttrain=output.T(n0+1:n0+handles.n,:);
    handles.Ttest=output.T(n0+handles.n+1:end,:);
    
    handles.model.T_params=T_params;
    %handles.model.T_params.base=output.T_params.base;
    handles.model.A_params=A_params;
    handles.model.A_params.base=output.A_params.base;
    handles.model.scenario=scenarioProcess;
    handles.model.refSigma=fault.refSigma;
    handles.model.m=handles.m;
    handles.model.p=handles.p;
    
    set(handles.status_hObject,'String','done')
    guidata(hObject, handles);
    
catch errorObj
        errordlg(errorObj.message,'Error');
%     errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
    
    set(handles.status_hObject,'String','error')
    guidata(hObject, handles);
end


function edit_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_m as text
%        str2double(get(hObject,'String')) returns contents of edit_m as a double
handles.m=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_m_hObject=hObject;
guidata(hObject, handles);


function edit_p_Callback(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_p as text
%        str2double(get(hObject,'String')) returns contents of edit_p as a double
handles.p=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_p_hObject=hObject;
guidata(hObject, handles);


function edit_n_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n as text
%        str2double(get(hObject,'String')) returns contents of edit_n as a double
handles.n=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_n_hObject=hObject;
guidata(hObject, handles);



function edit_n_noc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_noc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_noc as text
%        str2double(get(hObject,'String')) returns contents of edit_n_noc as a double
handles.n_noc=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_n_noc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_noc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_n_noc_hObject=hObject;
guidata(hObject, handles);


function edit_n_fault_Callback(hObject, eventdata, handles)
% hObject    handle to edit_n_fault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_n_fault as text
%        str2double(get(hObject,'String')) returns contents of edit_n_fault as a double
handles.n_fault=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_n_fault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_n_fault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_n_fault_hObject=hObject;
guidata(hObject, handles);


function edit_d_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d as text
%        str2double(get(hObject,'String')) returns contents of edit_d as a double
handles.d=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit_d_hObject=hObject;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_faultScenario.
function popupmenu_faultType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_faultScenario (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_faultScenario contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_faultScenario
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
handles.faultType=contents{get(hObject,'Value')};% returns selected item from popupmenu1
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_faultType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_faultScenario (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
type={'score'; 'sensor'};
set(hObject,'String',['-- Fault type --'; type])
handles.faultType_hObject=hObject;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_faultType.
function popupmenu_faultScenario_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_faultType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_faultType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_faultType
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
handles.faultScenario=contents{get(hObject,'Value')};% returns selected item from popupmenu1
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_faultScenario_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_faultType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
scenario={'step'; 'ramp'};
set(hObject,'String',['-- Fault scenario --'; scenario])
handles.faultScenario_hObject=hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if strcmp(handles.scenarioProcess,'Saved model')==0,
        assignin('base','Xtrain',handles.Xtrain);
        assignin('base','Ttrain',handles.Ttrain);
        assignin('base','model',handles.model);
    end
    assignin('base','Xtest',handles.Xtest);
    assignin('base','Ttest',handles.Ttest);
catch
    errordlg('Need to run the function first.','Error');
end



% --- Executes on selection change in popupmenu_model.
function popupmenu_model_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_model
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
index=get(hObject,'Value');
handles.nameModel=contents{index};% returns selected item from popupmenu1
if index~=1,
    handles.model=evalin('base',handles.nameModel);
end
vars=evalin('base','who');
set(hObject,'String',['-- Model --';vars])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
vars=evalin('base','who');
set(hObject,'String',['-- Model --';vars])
set(hObject,'Visible','Off'); 
handles.popupmenu_model_hObject=hObject;
handles.nameModel='';
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.status_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanel_train_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uipanel_train_hObject=hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    figure('name', 'Plot','NumberTitle','off')
    subplot(2,2,1)
    plot(handles.Ttrain)
    xlabel('Sample')
    ylabel('T_{train}')
    subplot(2,2,2)
    plot(handles.Xtrain)
    xlabel('Sample')
    ylabel('X_{train}')
    subplot(2,2,3)
    plot(handles.Ttest)
    xlabel('Sample')
    ylabel('T_{test}')
    subplot(2,2,4)
    plot(handles.Xtest)
    xlabel('Sample')
    ylabel('X_{test}')
catch
    close 'Plot'
    errordlg('Need to run the function first.','Error');
end


% --- Executes on selection change in popupmenuSTD.
function popupmenuSTD_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSTD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSTD
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
handles.refSigma=contents{get(hObject,'Value')};% returns selected item from popupmenu1
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuSTD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
options={'Signal variance'; 'Errors variance'};
set(hObject,'String',['-- Reference variance --'; options])
handles.refSigma_hObject=hObject;
guidata(hObject, handles);
