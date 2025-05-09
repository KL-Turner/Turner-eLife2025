function varargout = SelectBehavioralStateGUI_IOS_eLife2025(varargin)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: sleep scoring GUI button
%________________________________________________________________________________________________________________________
% SELECTBEHAVIORALSTATEGUI_IOS_eLife2025 MATLAB code for SelectBehavioralStateGUI_IOS_eLife2025.fig
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2025, by itself, creates a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2025 or raises the existing
%      singleton*.
%
%      H = SELECTBEHAVIORALSTATEGUI_IOS_eLife2025 returns the handle to a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2025 or the handle to
%      the existing singleton*.
%
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2025('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTBEHAVIORALSTATEGUI_IOS_eLife2025.M with the given input arguments.
%
%      SELECTBEHAVIORALSTATEGUI_IOS_eLife2025('Property','Value',...) creates a new SELECTBEHAVIORALSTATEGUI_IOS_eLife2025 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectBehavioralStateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectBehavioralStateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectBehavioralStateGUI_IOS_eLife2025

% Last Modified by GUIDE v2.5 18-Apr-2022 13:12:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SelectBehavioralStateGUI_OpeningFcn, ...
    'gui_OutputFcn',  @SelectBehavioralStateGUI_OutputFcn, ...
    'gui_LayoutFcn',  @SelectBehavioralStateGUI_IOS_eLife2025_LayoutFcn, ...
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


% --- Executes just before SelectBehavioralStateGUI is made visible.
function SelectBehavioralStateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectBehavioralStateGUI_IOS_eLife2025 (see VARARGIN)

% Choose default command line output for SelectBehavioralStateGUI_IOS_eLife2025
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectBehavioralStateGUI_IOS_eLife2025 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SelectBehavioralStateGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2025
% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2025
% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ButtonSelect_IOS_eLife2025
% Hint: get(hObject,'Value') returns toggle state of togglebutton3


% --- Creates and returns a handle to the GUI figure. 
function h1 = SelectBehavioralStateGUI_IOS_eLife2025_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end
load SelectBehavioralStateGUI_IOS_eLife2025.mat


appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', [], ...
    'text', 3, ...
    'togglebutton', 4), ...
    'override', 0, ...
    'release', [], ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'C:\Users\klt8\OneDrive - The Pennsylvania State University\Documents\Github Repositories\Data-Analysis\Pre-Processing Functions\IOS Sleep Scoring\SelectBehavioralStateGUI_IOS_eLife2025.m', ...
    'lastFilename', 'C:\Users\klt8\OneDrive - The Pennsylvania State University\Documents\Github Repositories\Data-Analysis\Pre-Processing Functions\IOS Sleep Scoring\SelectBehavioralStateGUI_IOS_eLife2025.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = mat{1};

h1 = figure(...
'Units','characters',...
'Position',[135.8 64.3076923076923 116.2 14.0769230769231],...
'Visible',get(0,'defaultfigureVisible'),...
'Color',get(0,'defaultfigureColor'),...
'IntegerHandle','off',...
'Colormap',get(0,'defaultfigureColormap'),...
'MenuBar','none',...
'Name','SelectBehavioralStateGUI_IOS_eLife2025',...
'NumberTitle','off',...
'HandleVisibility','callback',...
'CreateFcn', {@local_CreateFcn, @(fig,~)addToolbarExplorationButtons(fig), appdata} ,...
'CreateFcnMode',get(0,'defaultfigureCreateFcnMode'),...
'Tag','figure1',...
'Resize','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'ScreenPixelsPerInchMode','manual');

appdata = [];
appdata.lastValidTag = 'togglebutton1';

h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'String','Not Sleep',...
'Style','togglebutton',...
'Position',[2.6 2.46153846153846 36 4],...
'BackgroundColor',[0.0745098039215686 0.623529411764706 1],...
'Callback',@(hObject,eventdata)SelectBehavioralStateGUI_IOS_eLife2025('togglebutton1_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'ForegroundColor',get(0,'defaultuicontrolForegroundColor'),...
'Tag','togglebutton1',...
'FontSize',20,...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'togglebutton2';

h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'String','NREM Sleep',...
'Style','togglebutton',...
'Position',[39.6 2.46153846153846 36 4],...
'BackgroundColor',[0.5 1 0.25],...
'Callback',@(hObject,eventdata)SelectBehavioralStateGUI_IOS_eLife2025('togglebutton2_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'ForegroundColor',get(0,'defaultuicontrolForegroundColor'),...
'Tag','togglebutton2',...
'FontSize',20,...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'togglebutton3';

h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'String','REM Sleep',...
'Style','togglebutton',...
'Position',[76.6 2.46153846153846 36 4],...
'BackgroundColor',[1 0.75 0],...
'Callback',@(hObject,eventdata)SelectBehavioralStateGUI_IOS_eLife2025('togglebutton3_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','togglebutton3',...
'FontSize',20,...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text2';

h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'String','Select behavioral state',...
'Style','text',...
'Position',[24.6 7.38461538461539 66.2 3.92307692307692],...
'Children',[],...
'Tag','text2',...
'FontSize',24,...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error(message('MATLAB:guide:StateFieldNotFound', gui_StateFields{ i }, gui_Mfile));
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % SELECTBEHAVIORALSTATEGUI_IOS_eLife2025
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % SELECTBEHAVIORALSTATEGUI_IOS_eLife2025(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallback(gui_State, varargin{:})
    % SELECTBEHAVIORALSTATEGUI_IOS_eLife2025('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % SELECTBEHAVIORALSTATEGUI_JNEUROSCI2022(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~ishghandle(fig,'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || (isscalar(fig)&&isprop(fig,'GUIDEFigure'));
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishghandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI MATLAB code file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    % Call version of openfig that accepts 'auto' option"
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton, visible);  
%     %workaround for CreateFcn not called to create ActiveX
%         peers=findobj(findall(allchild(gui_hFigure)),'type','uicontrol','style','text');    
%         for i=1:length(peers)
%             if isappdata(peers(i),'Control')
%                 actxproxy(peers(i));
%             end            
%         end
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallback(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishghandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end
