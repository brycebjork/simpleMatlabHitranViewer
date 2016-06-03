function varargout = simpleMatlabHitranViewer(varargin)
% SIMPLEMATLABHITRANVIEWER MATLAB code for simpleMatlabHitranViewer.fig
%      SIMPLEMATLABHITRANVIEWER, by itself, creates a new SIMPLEMATLABHITRANVIEWER or raises the existing
%      singleton*.
%
%      H = SIMPLEMATLABHITRANVIEWER returns the handle to a new SIMPLEMATLABHITRANVIEWER or the handle to
%      the existing singleton*.
%
%      SIMPLEMATLABHITRANVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMPLEMATLABHITRANVIEWER.M with the given input arguments.
%
%      SIMPLEMATLABHITRANVIEWER('Property','Value',...) creates a new SIMPLEMATLABHITRANVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simpleMatlabHitranViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simpleMatlabHitranViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simpleMatlabHitranViewer

% Last Modified by GUIDE v2.5 01-Jun-2016 10:42:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simpleMatlabHitranViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @simpleMatlabHitranViewer_OutputFcn, ...
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


% --- Executes just before simpleMatlabHitranViewer is made visible.
function simpleMatlabHitranViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simpleMatlabHitranViewer (see VARARGIN)

% Choose default command line output for simpleMatlabHitranViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simpleMatlabHitranViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simpleMatlabHitranViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function atmosphericpressure_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to atmosphericpressure_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atmosphericpressure_textbox as text
%        str2double(get(hObject,'String')) returns contents of atmosphericpressure_textbox as a double


% --- Executes during object creation, after setting all properties.
function atmosphericpressure_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atmosphericpressure_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function partialpressure_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to partialpressure_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of partialpressure_textbox as text
%        str2double(get(hObject,'String')) returns contents of partialpressure_textbox as a double


% --- Executes during object creation, after setting all properties.
function partialpressure_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partialpressure_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function temperature_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to temperature_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of temperature_textbox as text
%        str2double(get(hObject,'String')) returns contents of temperature_textbox as a double


% --- Executes during object creation, after setting all properties.
function temperature_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temperature_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pathlength_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to pathlength_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathlength_textbox as text
%        str2double(get(hObject,'String')) returns contents of pathlength_textbox as a double


% --- Executes during object creation, after setting all properties.
function pathlength_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathlength_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadfile_button.
function loadfile_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Let's load in a file
    [filename, pathname] = uigetfile( ...
        {'*.mat','MATLAB files (*.mat)';...
         '*.par','HITRAN files (*.par)'}, ...
           'Choose a file...');
       
     [~,~,fileext] = fileparts(filename);
       
     if isempty(filename)
         return
     elseif strcmp(fileext,'.par')
         % The file is a *.par file - we need to convert it
        disp 'Loading File'
        fid = fopen (fullfile(pathname,filename),'rt'); %open file
        A = fscanf(fid,'%161c',[161 inf]); % read values
        fclose(fid);
        disp 'Done'
        
        %Transpose A
        A=A';

        s.igas= str2num(A(:,1:2));
        s.iso= str2num(A(:,3:3));
        s.wnum= str2num(A(:,4:15));
        s.int= str2num(A(:,16:25));
        s.Acoeff=str2num(A(:,26:35));
        s.abroad=str2num(A(:,36:40));
        s.sbroad=str2num(A(:,41:45));
        s.els=str2num(A(:,46:55));
        s.abcoef=str2num(A(:,56:59));
        s.tsp=str2num(A(:,60:67));
        s.gn=str2num(A(:,157:160));
        
        % Save to file
        [oldfilepath,oldfilename,~] = fileparts(fullfile(pathname,filename));
        newfilename = fullfile(oldfilepath,oldfilename,'.mat');
        [matfilename,matfilepath] = uiputfile('*.mat','Save MATLAB file',newfilename);
        if ~isempty(matfilename)
            save(fullfile(matfilepath,matfilename),'-struct','s');
        end
     else
         handles = guidata(hObject);
         handles.hitrandata = load(fullfile(pathname,filename));
         handles.hitranfilename = fullfile(pathname,filename);
         set(handles.loadedfile_statictext,'String',handles.hitranfilename);
         guidata(hObject,handles);
     end

% --- Executes on button press in simulate_button.
function simulate_button_Callback(hObject, eventdata, handles)
% hObject    handle to simulate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Let's plot something in the main plot

    handles = guidata(hObject);
    T = str2double(get(handles.temperature_textbox,'String')); % set temperature of gas for simulation
    pressure_torr = str2double(get(handles.atmosphericpressure_textbox,'String')); % input pressure in Torr for simulation
    partial_pressure_torr = str2double(get(handles.partialpressure_textbox,'String'));
    
    pressure_atm = pressure_torr / 760;
    partial_pressure_atm = partial_pressure_torr / 760;
    path_length_cm = str2double(get(handles.pathlength_textbox,'String')); % absorption pathlength
    N = str2double(get(handles.simulationnumpoints_textbox,'String')); % set this number higher to increase spectral resolution
    wavenumber_max = str2double(get(handles.simulationmax_textbox,'String'));
    wavenumber_min = str2double(get(handles.simulationmin_textbox,'String'));
    guidata(hObject,handles);
    
    wavelength_max_nm = 1e7/wavenumber_min ;
    wavelength_min_nm = 1e7/wavenumber_max;

    frequency_samples_wavenumber = 1e7./([wavelength_max_nm wavelength_min_nm]);

    df_wavenumber = (frequency_samples_wavenumber(2) - frequency_samples_wavenumber(1)) / (N - 1);

    frequency_samples_wavenumber = frequency_samples_wavenumber(1): df_wavenumber:frequency_samples_wavenumber(2);

    wavelength_samples_nm = 1e7 ./ frequency_samples_wavenumber;

    wavenumber_samples_1_per_cm = 1 ./ wavelength_samples_nm * 1e7;

    f0_wavenumber = min(frequency_samples_wavenumber) + (df_wavenumber * N / 2);


    % 14N2O mass

	isotopologues_array_ = [1];
	molecular_weight_array_amu = [14 + 14 + 16];

    % % 12CH4 mass
    % 	isotopologues_array_ = [1];
    % 	molecular_weight_array_amu = [14 + 14 + 16];
    % 
    % % 13CH4 mass
    % 	isotopologues_array_ = [2];
    % 	molecular_weight_array_amu = [14 + 14 + 16];

    %%%
    % Execute
    handles = guidata(hObject);
        absorbance = load_hitran_mat(handles.hitrandata, wavenumber_samples_1_per_cm, ...
                pressure_atm, partial_pressure_atm, ...
                path_length_cm, isotopologues_array_, molecular_weight_array_amu, T);
    % Plot vs. wavelength

    set(handles.mainplot,'XData',frequency_samples_wavenumber);
    set(handles.mainplot,'YData',absorbance);
    %set(handles.mainplot,'YData',frequency_samples_wavenumber);

% --- Executes during object creation, after setting all properties.
function mainaxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainaxes
    % Let's make a dummy plot - we'll replace it with data later
    handles = guidata(hObject);
    handles.mainplot = plot(NaN,NaN);
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function loadfile_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
         % Generate necessary variables
         disp test
         handles = guidata(hObject);
         handles.hitranfilename = [];
         guidata(hObject,handles);



function isotopenumber_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to isotopenumber_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isotopenumber_textbox as text
%        str2double(get(hObject,'String')) returns contents of isotopenumber_textbox as a double


% --- Executes during object creation, after setting all properties.
function isotopenumber_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isotopenumber_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isotopemass_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to isotopemass_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isotopemass_textbox as text
%        str2double(get(hObject,'String')) returns contents of isotopemass_textbox as a double


% --- Executes during object creation, after setting all properties.
function isotopemass_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isotopemass_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function simulationmin_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to simulationmin_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simulationmin_textbox as text
%        str2double(get(hObject,'String')) returns contents of simulationmin_textbox as a double


% --- Executes during object creation, after setting all properties.
function simulationmin_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulationmin_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function simulationmax_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to simulationmax_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simulationmax_textbox as text
%        str2double(get(hObject,'String')) returns contents of simulationmax_textbox as a double


% --- Executes during object creation, after setting all properties.
function simulationmax_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulationmax_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function simulationnumpoints_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to simulationnumpoints_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simulationnumpoints_textbox as text
%        str2double(get(hObject,'String')) returns contents of simulationnumpoints_textbox as a double


% --- Executes during object creation, after setting all properties.
function simulationnumpoints_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulationnumpoints_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
