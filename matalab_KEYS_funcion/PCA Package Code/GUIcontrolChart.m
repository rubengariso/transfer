function varargout = GUIcontrolChart(varargin)
% GUICONTROLCHART MATLAB code for GUIcontrolChart.fig
%      GUICONTROLCHART, by itself, creates a new GUICONTROLCHART or raises the existing
%      singleton*.
%
%      H = GUICONTROLCHART returns the handle to a new GUICONTROLCHART or the handle to
%      the existing singleton*.
%
%      GUICONTROLCHART('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUICONTROLCHART.M with the given input arguments.
%
%      GUICONTROLCHART('Property','Value',...) creates a new GUICONTROLCHART or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIcontrolChart_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIcontrolChart_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIcontrolChart

% Last Modified by GUIDE v2.5 16-Nov-2018 14:54:47

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIcontrolChart_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIcontrolChart_OutputFcn, ...
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


% --- Executes just before GUIcontrolChart is made visible.
function GUIcontrolChart_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIcontrolChart (see VARARGIN)

% Choose default command line output for GUIcontrolChart
handles.output = hObject;

handles.aG=1;
handles.monType='off-line';
handles.maxLag=10;
handles.eta=0.999;
handles.h=500;
handles.lagMethod='-- Type of selection --';
handles.finalPlot=0;
handles.plotMode=0;
handles.ff=0.02;
handles.kRef=0.5;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIcontrolChart wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIcontrolChart_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_Xtrain.
function popupmenu_Xtrain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Xtrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Xtrain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Xtrain

% if isfield(handles,'aG')==0,% default values
%     handles.aG=1;
%     handles.monType='off-line';
%     handles.maxLag=10;
%     handles.eta=0.999;
%     handles.h=500;
%     handles.lagMethod='-- Type of selection --';
% end

contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
index=get(hObject,'Value');
handles.nameXtrain=contents{index};% returns selected item from popupmenu1
if index~=1,
    handles.Xtrain=evalin('base',handles.nameXtrain);
end

% update popupmenus
vars=evalin('base','who');

set(handles.Xtrain_hObject,'String',['-- Train data set --';vars])
index=find(strcmp(handles.nameXtrain, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtrain_hObject,'Value',index)

set(handles.Xtest_hObject,'String',['-- Test data set --';vars])
index=find(strcmp(handles.nameXtest, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtest_hObject,'Value',index)

set(handles.savedModel_hObject,'String',['-- Saved model --';vars])
index=find(strcmp(handles.savedModel_nameModel, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.savedModel_hObject,'Value',index)


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_Xtrain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Xtrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
vars=evalin('base','who');
set(hObject,'String',['-- Train data set --';vars])
handles.Xtrain_hObject=hObject;
handles.nameXtrain='';
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_Xtest.
function popupmenu_Xtest_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Xtest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Xtest contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Xtest

contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
index=get(hObject,'Value');
handles.nameXtest=contents{index};% returns selected item from popupmenu1
if index~=1,
    handles.Xtest=evalin('base',handles.nameXtest);
end

% update popupmenus
vars=evalin('base','who');

set(handles.Xtrain_hObject,'String',['-- Train data set --';vars])
index=find(strcmp(handles.nameXtrain, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtrain_hObject,'Value',index)

set(handles.Xtest_hObject,'String',['-- Test data set --';vars])
index=find(strcmp(handles.nameXtest, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtest_hObject,'Value',index)

set(handles.savedModel_hObject,'String',['-- Saved model --';vars])
index=find(strcmp(handles.savedModel_nameModel, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.savedModel_hObject,'Value',index)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_Xtest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Xtest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
vars=evalin('base','who');
set(hObject,'String',['-- Test data set --';vars])
handles.Xtest_hObject=hObject;
handles.nameXtest='';
guidata(hObject, handles);

% --- Executes on button press in buttonRun.
function buttonRun_Callback(hObject, eventdata, handles)
% hObject    handle to buttonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status_hObject,'String','running')
pause (1e-3)

try
    % get updated values
    if get(handles.method_hObject,'Value')==1,
       error('Chose monitoring method.'); 
    end
    if get(handles.kSelect_hObject,'Value')==1 &&  (strcmp(handles.method,'PCA')==1 || strcmp(handles.method,'DPCA')==1 || strcmp(handles.method,'DPCADR')==1 || strcmp(handles.method,'RPCA')==1 || strcmp(handles.method,'MWPCA')==1),
       error('Chose method to select number of latent variables.'); 
    end
    if get(handles.preproc_hObject,'Value')==1 && (strcmp(handles.method,'PCA')==1 || strcmp(handles.method,'DPCA')==1 || strcmp(handles.method,'DPCADR')==1),
       error('Chose scaling method.'); 
    end
    if (get(handles.lagMethod_DPCA_hObject,'Value')==1 && strcmp(handles.method,'DPCA')==1) || (get(handles.lagMethod_DPCADR_hObject,'Value')==1 && strcmp(handles.method,'DPCADR')==1),
        error('Chose method to select the number of lags.');
    end
    if get(handles.typeCL_hObject,'Value')==1,
       error('Chose type of control limits.'); 
    end
    
    handles.h=str2double(get(handles.h_hObject,'String'));
    if isnan(handles.h)==1,
        handles.h=[];
    end
    handles.eta=str2num(get(handles.ff_hObject,'String'));
    if isnan(handles.eta)==1,
        handles.eta=[];
    end
    
    handles.trainValRatio=str2double(get(handles.trainValRatio_hObject,'String'));

    if strcmp(handles.method,'DPCA')==1,
        contents = cellstr(get(handles.lagMethod_DPCA_hObject,'String'));
        handles.lagMethod=contents{get(handles.lagMethod_DPCA_hObject,'Value')};
    elseif strcmp(handles.method,'DPCADR')==1,
        contents = cellstr(get(handles.lagMethod_DPCADR_hObject,'String'));
        handles.lagMethod=contents{get(handles.lagMethod_DPCADR_hObject,'Value')};
    else
        handles.lagMethod=[];
    end
    handles.maxLag=str2double(get(handles.maxLag_DPCA_hObject,'String'));% DPCA and DPCADR have the same maxLag
    handles.aG=str2double(get(handles.FDR_hObject,'String'));
    handles.cpvThresh=str2double(get(handles.cpvThresh_hObject,'String'));
    
    handles.ff=str2num(get(handles.ffcov_hObject,'String'));
    handles.kRef=str2num(get(handles.kRef_hObject,'String'));
    
    if strcmp(handles.method,'Saved model')==0,
        handles.Xtrain=evalin('base',handles.nameXtrain);
    else
        handles.Xtrain=handles.savedModel.preproc.loc;
    end
    handles.Xtest=evalin('base',handles.nameXtest);

    % pre-processing
    if strcmp(handles.method,'PCA')==1 || strcmp(handles.method,'DPCA')==1 || strcmp(handles.method,'DPCADR')==1,
        switch handles.preproc
            case 'None'
                loc=zeros(1,size(handles.Xtrain,2));
                sd=ones(1,size(handles.Xtrain,2));
            case 'Mean centering'
                loc=mean(handles.Xtrain);
                sd=ones(1,size(handles.Xtrain,2));
            case 'Auto-scaling'
                loc=mean(handles.Xtrain);
                sd=std(handles.Xtrain);
        end
    elseif strcmp(handles.method,'RPCA')==1 || strcmp(handles.method,'MWPCA')==1,
        loc=zeros(1,size(handles.Xtrain,2));
        sd=ones(1,size(handles.Xtrain,2));
    elseif strcmp(handles.method,'Saved model')==1
        loc=handles.savedModel.preproc.loc;
        sd=handles.savedModel.preproc.sd;
    else
        loc=zeros(1,size(handles.Xtrain,2));
        sd=ones(1,size(handles.Xtrain,2));
    end
    Ztrain=(handles.Xtrain-repmat(loc,size(handles.Xtrain,1),1))./repmat(sd,size(handles.Xtrain,1),1);
    Ztest=(handles.Xtest-repmat(loc,size(handles.Xtest,1),1))./repmat(sd,size(handles.Xtest,1),1);
    
    % run control chart
    switch handles.method
        case {'PCA','DPCA','DPCADR','RPCA','MWPCA'}
            [ handles.stat, handles.ucl, handles.model, handles.model_new ] = controlChart( Ztrain, Ztest, handles.method, handles.kSelect, handles.aG/100, 'monType', handles.monType, 'lagMode',handles.lagMethod,'maxLag',handles.maxLag,'eta',handles.eta,'h', handles.h, 'trainValRatio', handles.trainValRatio, 'cpvThresh', handles.cpvThresh/100, 'typeLimit', handles.typeCL );
        case 'Saved model'
            [ handles.stat, handles.ucl, handles.model, handles.model_new ] = controlChart( Ztrain, Ztest, [], [], [], 'model', handles.savedModel, 'monType', handles.monType );
        case {'MCUSUM','MEWMA','M1Z2','MAXD','LRT'}
            [ handles.stat, handles.ucl, handles.model, handles.model_new ] = controlChart( Ztrain, Ztest, handles.method, [], handles.aG/100, 'monType', handles.monType, 'eta',handles.ff,'kRef',handles.kRef );
    end
    handles.model.preproc.loc=loc;
    handles.model.preproc.sd=sd;
        
    % final plot
    if handles.finalPlot==1,
        statnames=fieldnames(handles.stat);
        sn=length(statnames);
        
        figure('name','Plot stat','NumberTitle','off')
        for i=1:sn,
            subplot(sn,1,i)
            plotChart( [], handles.stat.(statnames{i}), handles.ucl.(statnames{i}) );
            ylabel(statnames{i})
            if handles.plotMode==1,
                set(gca,'YScale','log');
            end
        end
    end
    
    set(handles.status_hObject,'String','done')
    guidata(hObject, handles);
    
catch errorObj
%     errordlg(errorObj.message,'Error');
    errordlg(getReport(errorObj,'extended','hyperlinks','off'),'Error');
    
    set(handles.status_hObject,'String','error')
    guidata(hObject, handles);
end




% --- Executes on selection change in popupmenu_kSelect.
function popupmenu_kSelect_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kSelect


contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
kSelect=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.kSelect=kSelect;
if strcmp(kSelect,'kSelectCpv')==1,
    set(handles.text_cpvThresh_hObject,'Visible','On'); 
    set(handles.text_cpvPerc_hObject,'Visible','On'); 
    set(handles.cpvThresh_hObject,'Visible','On'); 
else
    set(handles.text_cpvThresh_hObject,'Visible','Off');
    set(handles.text_cpvPerc_hObject,'Visible','Off'); 
    set(handles.cpvThresh_hObject,'Visible','Off');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_kSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','-- Method to select LV --')
handles.kSelect_hObject=hObject;
guidata(hObject, handles);


function edit_FDR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FDR as text
%        str2double(get(hObject,'String')) returns contents of edit_FDR as a double

handles.aG=str2double(get(hObject,'String'));
guidata(hObject, handles);
       

% --- Executes during object creation, after setting all properties.
function edit_FDR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FDR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.FDR_hObject=hObject;
guidata(hObject, handles);

% --- Executes on button press in checkbox_onlinePlot.
function checkbox_onlinePlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_onlinePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_onlinePlot
do_plot=get(hObject,'Value');
if do_plot==1,
    handles.monType='on-line';
else
    handles.monType='off-line';
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_statName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_statName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function RPCA_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RPCA_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% contents = get(hObject,'String'); 
set(hObject,'Visible','Off'); 
handles.RPCA_panel=hObject;
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_Method.
function popupmenu_Method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Method
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
method=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.method=method;

set(handles.PCA_panel,'Visible','Off'); 
set(handles.DPCA_panel,'Visible','Off'); 
set(handles.DPCADR_panel,'Visible','Off'); 
set(handles.RPCA_panel,'Visible','Off'); 
set(handles.MWPCA_panel,'Visible','Off'); 
set(handles.Saved_panel,'Visible','Off'); 
set(handles.MCUSUM_panel,'Visible','Off'); 
set(handles.FF_panel,'Visible','Off'); 
set(handles.text_cpvThresh_hObject,'Visible','Off'); 
set(handles.text_cpvPerc_hObject,'Visible','Off'); 
set(handles.cpvThresh_hObject,'Visible','Off');

switch method
    case {'PCA','DPCA','DPCADR'}
        
        set(handles.([method,'_panel']),'Visible','On');
        kSelectMethods={'kSelectCpv'; 'kSelectCrossVal'; 'kSelectParAnl'};
        preproc={'None'; 'Mean centering'; 'Auto-scaling'};
        typeCL={'Theoretic'; 'Empiric'};
        kSelect_visible='On';
        preproc_visible='On';
        FDR_visible='On';
        typeCL_visible='On';
        
    case {'RPCA','MWPCA'}
        
        set(handles.([method,'_panel']),'Visible','On');
        kSelectMethods={'kSelectCpv'};
        preproc=[];
        typeCL={'Theoretic'; 'Empiric'};
        kSelect_visible='On';
        preproc_visible='Off';
        FDR_visible='On';
        typeCL_visible='On';
        
    case 'Saved model'
        
        set(handles.Saved_panel,'Visible','On');
        kSelectMethods=[];
        preproc=[];
        typeCL=[];
        kSelect_visible='Off';
        preproc_visible='Off';
        FDR_visible='Off';
        typeCL_visible='Off';
        
        % update popupmenu
        vars=evalin('base','who');
        set(handles.savedModel_hObject,'String',['-- Saved model --';vars]);
        index=find(strcmp(handles.savedModel_nameModel, vars))+1;
        if isempty(index)==1,
            index=1;
        end
        set(handles.savedModel_hObject,'Value',index)
        
    case {'MCUSUM'}
        
        set(handles.([method,'_panel']),'Visible','On');
        kSelectMethods=[];
        preproc=[];
        typeCL={'Empiric'};
        kSelect_visible='Off';
        preproc_visible='Off';
        FDR_visible='On';
        typeCL_visible='On';
        
    case {'MEWMA','M1Z2','MAXD','LRT'}
        
        set(handles.FF_panel,'Visible','On');
        set(handles.FF_panel,'Title',method);
        kSelectMethods=[];
        preproc=[];
        typeCL={'Empiric'};
        kSelect_visible='Off';
        preproc_visible='Off';
        FDR_visible='On';
        typeCL_visible='On';
        
end
set(handles.kSelect_hObject,'String',['-- Method to select LV --'; kSelectMethods])
set(handles.kSelect_hObject,'Value',1)
handles.kSelect='-- Method to select LV --';
set(handles.preproc_hObject,'String',['-- Pre-processing --'; preproc])
set(handles.preproc_hObject,'Value',1)
handles.preproc='-- Pre-processing --';
set(handles.typeCL_hObject,'String',['-- Type of limits --'; typeCL])
set(handles.typeCL_hObject,'Value',1)
handles.typeCL='-- Type of limits --';
set(handles.kSelect_hObject,'Visible',kSelect_visible);
set(handles.preproc_hObject,'Visible',preproc_visible);
set(handles.FDR_hObject,'Visible',FDR_visible);
set(handles.FDR_t1_hObject,'Visible',FDR_visible);
set(handles.FDR_t2_hObject,'Visible',FDR_visible);
set(handles.typeCL_hObject,'Visible',typeCL_visible);

    
% set(handles.maxLag_DPCA_hObject,'String',handles.maxLag);
% set(handles.maxLag_DPCADR_hObject,'String',handles.maxLag);

lagMethod={'auto'; 'manual'};
set(handles.lagMethod_DPCA_hObject,'String',['-- Type of selection --'; lagMethod]);
set(handles.lagMethod_DPCADR_hObject,'String',['-- Type of selection --'; lagMethod]);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_Method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
methods={'PCA'; 'DPCA'; 'DPCADR';'RPCA';'MWPCA';'Saved model';'MCUSUM';'MEWMA';'M1Z2';'MAXD';'LRT'};
set(hObject,'String',['-- PCA method --'; methods])
handles.method_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PCA_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PCA_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.PCA_panel=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function DPCA_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DPCA_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.DPCA_panel=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function DPCADR_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DPCADR_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.DPCADR_panel=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MWPCA_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MWPCA_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.MWPCA_panel=hObject;
guidata(hObject, handles);



function edit_maxLag_DPCA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxLag_DPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxLag_DPCA as text
%        str2double(get(hObject,'String')) returns contents of edit_maxLag_DPCA as a double
handles.maxLag=str2double(get(hObject,'String'));
set(handles.maxLag_DPCADR_hObject,'String',handles.maxLag);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_maxLag_DPCA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxLag_DPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.maxLag_DPCA_hObject=hObject;
guidata(hObject, handles);



function edit_maxLag_DPCADR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxLag_DPCADR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxLag_DPCADR as text
%        str2double(get(hObject,'String')) returns contents of edit_maxLag_DPCADR as a double
handles.maxLag=str2double(get(hObject,'String'));
set(handles.maxLag_DPCA_hObject,'String',handles.maxLag);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_maxLag_DPCADR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxLag_DPCADR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.maxLag_DPCADR_hObject=hObject;
guidata(hObject, handles);



function edit_forgettingFactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_forgettingFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_forgettingFactor as text
%        str2double(get(hObject,'String')) returns contents of edit_forgettingFactor as a double
handles.eta=str2num(get(hObject,'String'));
if isnan(handles.eta)==1,
    handles.eta=[];
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_forgettingFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_forgettingFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.ff_hObject=hObject;
guidata(hObject, handles);



function edit_window_Callback(hObject, eventdata, handles)
% hObject    handle to edit_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_window as text
%        str2double(get(hObject,'String')) returns contents of edit_window as a double
handles.h=str2double(get(hObject,'String'));
if isnan(handles.h)==1,
    handles.h=[];
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.h_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanel_savedModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_savedModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.Saved_panel=hObject;
guidata(hObject, handles);



% --- Executes on selection change in popupmenu_lagModeDPCA.
function popupmenu_lagModeDPCA_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_lagModeDPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_lagModeDPCA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_lagModeDPCA
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
lagMethod=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.lagMethod=lagMethod;
% set(handles.lagMethod_DPCADR_hObject,'String',handles.lagMethod);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_lagModeDPCA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_lagModeDPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.lagMethod_DPCA_hObject=hObject;
guidata(hObject, handles);

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
lagMethod={'auto'; 'manual'};
set(hObject,'String',['-- Type of selection --'; lagMethod])
handles.lagMethod_DPCA_hObject=hObject;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_lagModeDPCADR.
function popupmenu_lagModeDPCADR_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_lagModeDPCADR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_lagModeDPCADR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_lagModeDPCADR
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
lagMethod=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.lagMethod=lagMethod;
% set(handles.lagMethod_DPCA_hObject,'String',handles.lagMethod);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_lagModeDPCADR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_lagModeDPCADR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
lagMethod={'auto'; 'manual'};
set(hObject,'String',['-- Type of selection --'; lagMethod])
handles.lagMethod_DPCADR_hObject=hObject;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_savedModel.
function popupmenu_savedModel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_savedModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_savedModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_savedModel

contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
index=get(hObject,'Value');
handles.savedModel_nameModel=contents{index};% returns selected item from popupmenu1
if index~=1,
    handles.savedModel=evalin('base',handles.savedModel_nameModel);
end

% update popupmenus
vars=evalin('base','who');

set(handles.Xtrain_hObject,'String',['-- Train data set --';vars])
index=find(strcmp(handles.nameXtrain, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtrain_hObject,'Value',index)

set(handles.Xtest_hObject,'String',['-- Test data set --';vars])
index=find(strcmp(handles.nameXtest, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtest_hObject,'Value',index)

set(handles.savedModel_hObject,'String',['-- Saved model --';vars])
index=find(strcmp(handles.savedModel_nameModel, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.savedModel_hObject,'Value',index)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_savedModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_savedModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
vars=evalin('base','who');
set(hObject,'String',['-- Saved model --';vars])
handles.savedModel_hObject=hObject;
handles.savedModel_nameModel='';
guidata(hObject, handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','stat',handles.stat);
assignin('base','ucl',handles.ucl);
assignin('base','model',handles.model);


% --- Executes during object creation, after setting all properties.
function text_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.status_hObject=hObject;
guidata(hObject, handles);


% --- Executes on button press in checkbox_finalPlot.
function checkbox_finalPlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_finalPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_finalPlot
handles.finalPlot=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure('name','Plot stat','NumberTitle','off')
subplot(2,1,1)
plotChart( [], handles.stat.sT, handles.ucl.sT );
if handles.plotMode==1,
    set(gca,'YScale','log');
end
ylabel('sT')
subplot(2,1,2)
plotChart( [], handles.stat.sE, handles.ucl.sE );
ylabel('sE')
if handles.plotMode==1,
    set(gca,'YScale','log');
end


% --- Executes on selection change in popupmenu_preproc.
function popupmenu_preproc_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_preproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_preproc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_preproc
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
preproc=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.preproc=preproc;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_preproc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_preproc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
preproc={'None'; 'Mean centering'; 'Auto-scaling'};
set(hObject,'String',['-- Pre-processing --'; preproc])
handles.preproc_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_FDR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_FDR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.FDR_t1_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_FDR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_FDR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.FDR_t2_hObject=hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton_genData.
function pushbutton_genData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_genData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiTDseries

% % update popupmenus
% vars=evalin('base','who');
% 
% set(handles.Xtrain_hObject,'String',['-- Train data set --';vars])
% index=find(strcmp(handles.nameXtrain, vars))+1;
% if isempty(index)==1,
%     index=1;
% end
% set(handles.Xtrain_hObject,'Value',index)
% 
% set(handles.Xtest_hObject,'String',['-- Test data set --';vars])
% index=find(strcmp(handles.nameXtest, vars))+1;
% if isempty(index)==1,
%     index=1;
% end
% set(handles.Xtest_hObject,'Value',index)
% 
% set(handles.savedModel_hObject,'String',['-- Saved model --';vars])
% index=find(strcmp(handles.savedModel_nameModel, vars))+1;
% if isempty(index)==1,
%     index=1;
% end
% set(handles.savedModel_hObject,'Value',index)
% 
% guidata(hObject, handles);

% --- Executes on button press in pushbutton_loadData.
function pushbutton_loadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FILENAME, PATHNAME] = uigetfile('*.mat', 'Load data');

X=load([PATHNAME,FILENAME]);
vars=fieldnames(X);
for i=1:length(vars),
    assignin('base',vars{i},X.(vars{i}))
end

% update names
vars=evalin('base','who');

set(handles.Xtrain_hObject,'String',['-- Train data set --';vars])
index=find(strcmp(handles.nameXtrain, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtrain_hObject,'Value',index)

set(handles.Xtest_hObject,'String',['-- Test data set --';vars])
index=find(strcmp(handles.nameXtest, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.Xtest_hObject,'Value',index)

set(handles.savedModel_hObject,'String',['-- Saved model --';vars])
index=find(strcmp(handles.savedModel_nameModel, vars))+1;
if isempty(index)==1,
    index=1;
end
set(handles.savedModel_hObject,'Value',index)

guidata(hObject, handles);



function ff_Callback(hObject, eventdata, handles)
% hObject    handle to ff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ff as text
%        str2double(get(hObject,'String')) returns contents of ff as a double
handles.ff=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.ffcov_hObject=hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FF_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FF_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.FF_panel=hObject;
guidata(hObject, handles);



function kRef_Callback(hObject, eventdata, handles)
% hObject    handle to kRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kRef as text
%        str2double(get(hObject,'String')) returns contents of kRef as a double
handles.kRef=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function kRef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.kRef_hObject=hObject;
guidata(hObject, handles);



function edit_trainValRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trainValRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trainValRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_trainValRatio as a double
handles.trainValRatio=str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_trainValRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trainValRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.trainValRatio_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MCUSUM_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MCUSUM_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','Off'); 
handles.MCUSUM_panel=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Data_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Data_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Method_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Method_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function checkbox_onlinePlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_onlinePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.monType='off-line';
guidata(hObject, handles);


function cpvThresh_Callback(hObject, eventdata, handles)
% hObject    handle to cpvThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpvThresh as text
%        str2double(get(hObject,'String')) returns contents of cpvThresh as a double
handles.cpvThresh=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cpvThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpvThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.cpvThresh_hObject=hObject;
set(hObject,'Visible','Off'); 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_cpv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_cpv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.text_cpvThresh_hObject=hObject;
set(hObject,'Visible','Off'); 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text_cpvPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_cpvPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.text_cpvPerc_hObject=hObject;
set(hObject,'Visible','Off'); 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function checkbox_finalPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_finalPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_logscale.
function checkbox_logscale_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_logscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_logscale
handles.plotMode=get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on selection change in popupmenuTypeCL.
function popupmenuTypeCL_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuTypeCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuTypeCL contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuTypeCL
contents = cellstr(get(hObject,'String'));% returns popupmenu1 contents as cell array
typeCL=contents{get(hObject,'Value')};% returns selected item from popupmenu1
handles.typeCL=typeCL;
if strcmp(typeCL,'Theoretic')==1,
    set(handles.textTrainVal_t1_hObject,'Visible','Off');
    set(handles.textTrainVal_t2_hObject,'Visible','Off');
    set(handles.trainValRatio_hObject,'Visible','Off');
else
    set(handles.textTrainVal_t1_hObject,'Visible','On');
    set(handles.textTrainVal_t2_hObject,'Visible','On');
    set(handles.trainValRatio_hObject,'Visible','On');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenuTypeCL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuTypeCL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
typeCL={'Theoretic'; 'Empiric'};
set(hObject,'String',['-- Type of limits --'; typeCL])
handles.typeCL_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function textTrainVal_t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textTrainVal_t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.textTrainVal_t1_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function textTrainVal_t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textTrainVal_t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.textTrainVal_t2_hObject=hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Options_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Options_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when Data_panel is resized.
function Data_panel_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to Data_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function checkbox_logscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_logscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.plotMode=get(hObject,'Value');
guidata(hObject, handles);