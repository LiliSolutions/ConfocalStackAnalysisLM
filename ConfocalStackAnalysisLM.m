function varargout = ConfocalStackAnalysisLM(varargin)
%% Lilia Mesina, Polaris/CCBN, Sptember2014
% Last Modified by LM, 23Sept2014

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConfocalStackAnalysisLM_OpeningFcn, ...
                   'gui_OutputFcn',  @ConfocalStackAnalysisLM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   @ConfocalStackAnalysisLM_Callback);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before ConfocalStackAnalysisLM is made visible.
function ConfocalStackAnalysisLM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConfocalStackAnalysisLM (see VARARGIN)

% Choose default command line output for ConfocalStackAnalysisLM
global step1 step2 step3 batchModeOn
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%set(hObject, 'KeyPressFcn',@fKeyPress);
movegui(hObject,'northwest');

axis(handles.axes1, 'ij', 'on','square');
step1 = 0; step2 = 0; step3 = 0;
batchModeOn = 0;

home = pwd;
addpath(genpath(home));
fprintf(2, 'ConfocalStackAnalysisLM Path Loaded!');


%% --- Outputs from this function are returned to the command line.
function varargout = ConfocalStackAnalysisLM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function ConfocalStackAnalysisLM_Callback(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConfocalStackAnalysisLM (see VARARGIN)

% Choose default command line output for ConfocalStackAnalysisLM

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%set(hObject, 'KeyPressFcn', @fKeyPress);


% --- Executes on button press in loadImageStack.
function loadImageStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadImageStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global activeFrame step1 step2 step3
global fileName pathName fullFN batchModeOn
global elems statFr zxy_stack num_images
global newGreenS newBlueS 

sel = get(hObject,'Value');
if sel
    if ~batchModeOn
         [fileName, pathName] = uigetfile({'*.tif';'*.*'},'Stack Selector');
         if isequal(fileName,0)
           disp('User selected Cancel');
           set(hObject,'Value',0);
           return
         end
    end
    fullFN = fullfile(pathName, fileName);
    disp('==');     
    disp(['==Processing ' fullFN]);     
    caption = sprintf('ConfocalStackAnalysisLMv2.0 >>%s', fullFN);
    set(gcf, 'Name', caption );

    %% load stack
    info = imfinfo(fullFN);
    [img_size1 img_size2 img_size3]=size(imread(fullFN, 1, 'Info', info));
    img_type=['uint' num2str(info(1).BitDepth)];
    num_images = numel(info);
    elems = 1024.^2;
    statFr = zeros(8,num_images);
    zxy_stack = zeros([num_images img_size1 img_size2 img_size3],'uint8'); %, img_type);
    newGreenS = zeros([num_images img_size1 img_size2]);
    newBlueS = zeros([num_images img_size1 img_size2]);

    for k = 1:num_images
        A = imread(fullFN, k, 'Info', info);
        zxy_stack(k,:,:,:) = A;
        activeFrame = k;
        draw_frame(hObject, eventdata, handles);
        %t = isequal(A,frame);
        getStats1(A,k);
        set(handles.uitable1,'data',statFr);
        pause(0.01);
    end
    set_tresh;
    activeFrame = 1;
    step1 = 1;
    draw_frame(hObject, eventdata, handles);
    h = msgbox('Load confocal stack Done.');
    pause(1); close(h);
end


function runAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to runAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global statFr zxy_stack num_images
global step1 step2

sel = get(hObject,'Value');
if sel
    if ~step1
        msgbox('Run step1 first.');
        return
    end
    hB = waitbar(0, 'Image Analysis in progress...'); 
    for i = 1:num_images
        fr = zxy_stack(i,:,:,:);
        fr = squeeze(fr);
        getStats2(fr,i);
        set(handles.uitable1,'data',statFr);
        waitbar(i/num_images, hB, 'Image Analysis in progress...');    
    end
    stdDev = std(statFr,0,2);
    statFr(:,num_images+1:end) = [];
    statFr = horzcat(statFr,stdDev);
    rowMean = mean(statFr(:,1:num_images),2);
    statFr = horzcat(statFr,rowMean);
    set(handles.uitable1,'data',statFr);    
    step2 = 1;
    waitbar(1, hB, 'Image Analysis in progress...');   
    close(hB);
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global step2 statFr fullFN num_images

sel = get(hObject,'Value');
if sel
    if ~step2
        msgbox('Run step2 first to get image analysis done.');
        return
    end
    [dn imgName ext] = fileparts(fullFN);
    rows1 = {'_','','';'','Frames:',''}; 
    rows1{2,1} = imgName;
    rows1{2,3} = num_images;
    rows2 = {'MinGreen';'MaxGreen';'MinBlue';'MaxBlue';'MeanGreen';...
             'GreenPixels';'BluePixels';'G/B'};
    mat2 = num2cell(statFr);
    rows2 = horzcat(rows2,mat2);
    cellArray = strcat({'Frame'},int2str((1:num_images).')).';
    colName = horzcat({'_',},cellArray);
    colName = horzcat(colName,{'stdDev'});
    colName = horzcat(colName,{'stackMean'});
    rows2 = vertcat(colName,rows2);
    fn = [dn '\ConfocalStats.xls'];
    if ~exist(fn,'file')    
        xlswrite(fn,rows1);
        xlsappend(fn,rows2,1);
    else
        %rows2 = vertcat(rows1,rows2);
        xlsappend(fn,rows1,1);
        xlsappend(fn,rows2,1);   
    end
    disp(['Save results to excel Done. File ' fn]);
end


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear all;
clc;
close all;
h = msgbox('Close Analysis Done.');
pause(1);
close(h);
run ConfocalStackAnalysisLM;

function first_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global activeFrame
    
activeFrame = 1;
draw_frame(hObject, eventdata, handles);


function last_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global activeFrame num_images
    
activeFrame = num_images;
draw_frame(hObject, eventdata, handles);


function prev_frame_Callback(hObject, eventdata, handles)
% hObject    handle to prev_frame_ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global activeFrame num_images
    
activeFrame = activeFrame - 1;
if activeFrame == 0
    activeFrame = num_images;
end
draw_frame(hObject, eventdata, handles);


function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global activeFrame num_images
    
activeFrame = activeFrame + 1;
if activeFrame > num_images
    activeFrame = 1;
end
draw_frame(hObject, eventdata, handles);




%% Run GUI funcs
function draw_frame(hObject, eventdata, handles)
global zxy_stack activeFrame num_images

axis(handles.axes1, 'ij', 'on','square');
axes(handles.axes1);
disp(['activeFrame: ' num2str(activeFrame)]);
thisFrame = zxy_stack(activeFrame,:,:,:);
thisFrame = squeeze(thisFrame);
image(thisFrame);
drawnow;
caption = ['Image ' num2str(activeFrame) '/' num2str(num_images)];
title(handles.axes1,caption,'Interpreter','none','Units', 'normalized', ...
'Position', [0 1], 'HorizontalAlignment', 'left');
update_fig(hObject, eventdata, handles);
pause(0.1);


function update_fig(hObject, eventdata, handles)
hFig = findobj('Tag','newGreenFig');
if ~isempty(hFig)
    ConfocalStackAnalysisLM('ShowNewGreen_Callback',hObject,eventdata,guidata(hObject));
end
hFig = findobj('Tag','newBlueFig');
if ~isempty(hFig)
    ConfocalStackAnalysisLM('ShowNewBlue_Callback',hObject,eventdata,guidata(hObject));
end


function getStats1(fr,k)
global elems statFr

g = fr(:,:,2); 
b = fr(:,:,3);

%% get image profile
bM = reshape(b,[1,elems]);
gM = reshape(g,[1,elems]);
bM = double(bM);
gM = double(gM);
%% get stats
statFr(1,k) = min(gM);
statFr(2,k) = max(gM);
statFr(3,k) = min(bM);
statFr(4,k) = max(bM);


function getStats2(fr,k)
global statFr
global newGreenS newBlueS
global greenMinTresh
global blueMinTresh blueMaxTresh

r = fr(:,:,1); g = fr(:,:,2); b = fr(:,:,3);

% bM = reshape(b,[1,elems]);
% %m0(m0==0) = [];
% gM = reshape(g,[1,elems]);
% %m0(m0==0) = [];

%% get newGreen
 [n1, n2] = size(g(:,:));
 newGreen = g;
 for i = 1:n1
    for j = 1:n2
        if (b(i,j) < blueMinTresh) || (g(i,j) <= greenMinTresh) %(blue(i,j) < 20.4) || (green(i,j) <= 12)
            newGreen(i,j) = 0;
        end
    end
 end
 
 %image(newGreen(:,:));
 green2 = find(newGreen>0);
 nG = length(green2);
 nG2 = nG; %/elems;
 statFr(6,k) = nG2;

 %% get newBlue
 [n1, n2] = size(b(:,:));
 newBlue = b;
 for i = 1:n1
    for j = 1:n2
        if (b(i,j) < blueMinTresh) || (b(i,j) > blueMaxTresh)%|| (green(i,j) <= 12)
            newBlue(i,j) = 0;
        end
    end
 end
 %figure; imagesc(newBlue(:,:));
blue2 = find(newBlue>0);
nB = length(blue2);
nB2 = nB;%/elems;
statFr(7,k) = nB2;
statFr(8,k) = nG2/nB2;
newGreenS(k,:,:) = newGreen;
newBlueS(k,:,:) = newBlue;
%getStats1(fr,k);


function set_tresh()
global elems num_images statFr zxy_stack
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip

modeMultip = 1.2;       hObj = findobj('Tag', 'modeM');     set(hObj, 'String', num2str(modeMultip));
greenMinTresh = 12;     hObj = findobj('Tag', 'gMinTresh'); set(hObj, 'String', num2str(greenMinTresh));
greenMaxTresh = 255;    hObj = findobj('Tag', 'gMaxTresh'); set(hObj, 'String', num2str(greenMaxTresh));
blueMinTresh = 30;      hObj = findobj('Tag', 'bMinTresh'); set(hObj, 'String', num2str(blueMinTresh));
blueMaxTresh = 230;     hObj = findobj('Tag', 'bMaxTresh'); set(hObj, 'String', num2str(blueMaxTresh));

gMinTresh = zeros(num_images,1);
bMinTresh = zeros(num_images,1);
hHist = figure();

for i = 1:num_images
    fr = zxy_stack(i,:,:,:);
    fr = squeeze(fr);
    %fr3 = imagesc(fr);
    g = fr(:,:,2); %imagesc(
    b = fr(:,:,3); %imagesc(
    
%     %% clean Glia
% figure
% h2 =reshape(b,[1,elems]);
% h2 = double(h2);
% a1 = bar(h2); %,64);
% plot(a1);
% %bar(log(a1));
% figure
% bar(a1);
% shg


    h2 = reshape(g,[1,elems]);
    %? make it double
    statFr(5,i) = mean(h2);
    h3 = imhist(g);
    figure(hHist);
    bar(h3); pause(0.5);
    gMinTresh(i) = modeMultip * mode(h2);
    h2 = reshape(b,[1,elems]);
    h3 = imhist(b);
    figure(hHist);
    bar(h3); pause(0.5);
    bMinTresh(i) = modeMultip * mode(h2);
end
close(hHist);
greenMinTresh = mean(gMinTresh);
blueMinTresh = mean(bMinTresh);

hObj = findobj('Tag', 'gMinTresh'); set(hObj, 'String', num2str(greenMinTresh));
%greenMaxTresh = 255;    hObj = findobj('Tag', 'gMaxTresh'); set(hObj, 'String', num2str(greenMaxTresh));
hObj = findobj('Tag', 'bMinTresh'); set(hObj, 'String', num2str(blueMinTresh));
%hObj = findobj('Tag', 'bMaxTresh'); set(hObj, 'String', num2str(blueMaxTresh));

%% bMaxTresh = ; ??  to skip glia
% % h2 = double(h2); 
% % a1 = hist(h2,64); %,elems); %
% % figure; bar(log(a1));
% % shg
%
% % blue2 = find(newBlue>0);
% % nB = length(blue2);
% % nB2 = nB/elems;


function playB_Callback(hObject, eventdata, handles)
% hObject    handle to playB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global activeFrame num_images
    
if isempty(activeFrame)
    startFr = 1;
else
    startFr = activeFrame;
end
for i = startFr:num_images
    activeFrame = i;
    draw_frame(hObject, eventdata, handles);
    pause(0.5);
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object deletion, before destroying properties.
function uitable1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ShowNewGreen.
function ShowNewGreen_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNewGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newGreenS activeFrame statFr

sel = get(hObject,'Value');
if sel
    if ~mean(statFr(5,:))
        hM = msgbox('First run Step 2: Image Analysis.');
        set(hObject,'Value',0);
        return
    end
    hFig = findobj('Tag','newGreenFig');
    if isempty(hFig)
        hFig = figure('Tag','newGreenFig','Name','New Green','NumberTitle','off');
    end
    figure(hFig);
    caption = ['New Green: Frame ' num2str(activeFrame)];
    set(hFig,'Name',caption);
    im = newGreenS(activeFrame,:,:);
    im = squeeze(im);
    imagesc(im);
end

% --- Executes on button press in ShowNewBlue.
function ShowNewBlue_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNewBlue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newBlueS activeFrame statFr

sel = get(hObject,'Value');
if sel
    if ~mean(statFr(5,:))
        hM = msgbox('First run Step 2: Image Analysis.');
        set(hObject,'Value',0);
        return
    end    
    hFig = findobj('Tag','newBlueFig');
    if isempty(hFig)
        hFig = figure('Tag','newBlueFig','Name','New Blue','NumberTitle','off');
    end
    figure(hFig);
    caption = ['New Blue: Frame ' num2str(activeFrame)];
    set(hFig,'Name',caption);
    im = newBlueS(activeFrame,:,:);
    im = squeeze(im);
    imagesc(im);
end


function gMinTresh_Callback(hObject, eventdata, handles)
% hObject    handle to gMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gMinTresh as text
%        str2double(get(hObject,'String')) returns contents of gMinTresh as a double
global greenMinTresh greenMaxTresh
global blueMinTresh blueMaxTresh

val = get(handles.gMinTresh,'string');
greenMinTresh = str2double(val); 
disp(['greenMinTresh ' val]);


% --- Executes during object creation, after setting all properties.
function gMinTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gMaxTresh_Callback(hObject, eventdata, handles)
% hObject    handle to gMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gMaxTresh as text
%        str2double(get(hObject,'String')) returns contents of gMaxTresh as a double
global greenMinTresh greenMaxTresh
global blueMinTresh blueMaxTresh

val = get(handles.gMaxTresh,'string');
greenMaxTresh = str2double(val); 
disp(['greenMaxTresh ' val]);

% --- Executes during object creation, after setting all properties.
function gMaxTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bMinTresh_Callback(hObject, eventdata, handles)
% hObject    handle to bMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bMinTresh as text
%        str2double(get(hObject,'String')) returns contents of bMinTresh as a double
global greenMinTresh greenMaxTresh
global blueMinTresh blueMaxTresh

val = get(handles.bMinTresh,'string');
blueMinTresh = str2double(val); 
disp(['blueMinTresh ' val]);

% --- Executes during object creation, after setting all properties.
function bMinTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bMaxTresh_Callback(hObject, eventdata, handles)
% hObject    handle to bMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bMaxTresh as text
%        str2double(get(hObject,'String')) returns contents of bMaxTresh as a double
global greenMinTresh greenMaxTresh
global blueMinTresh blueMaxTresh

val = get(handles.bMaxTresh,'string');
blueMaxTresh = str2double(val); 
disp(['blueMaxTresh ' val]);

% --- Executes during object creation, after setting all properties.
function bMaxTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function modeM_Callback(hObject, eventdata, handles)
% hObject    handle to modeM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of modeM as text
%        str2double(get(hObject,'String')) returns contents of modeM as a double
global modeMultip

val = get(handles.modeM,'string');
modeMultip = str2double(val); 
disp(['modeMultip ' val]);

% --- Executes during object creation, after setting all properties.
function modeM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modeM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runAll.
function runAll_Callback(hObject, eventdata, handles)
% hObject    handle to runAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fullFN

ConfocalStackAnalysisLM('loadImageStack_Callback',hObject,eventdata,guidata(hObject));
ConfocalStackAnalysisLM('runAnalysis_Callback',hObject,eventdata,guidata(hObject));
ConfocalStackAnalysisLM('save_Callback',hObject,eventdata,guidata(hObject));
disp('Run all steps Done.');
%pause(0.2);close(h);


% --- Executes on button press in batchMode.
function batchMode_Callback(hObject, eventdata, handles)
% hObject    handle to batchMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fileName pathName fullFN batchModeOn batchRoot home

% select batch root global batchRoot
%
batchModeOn = 1;
batchRoot = uigetdir(home,'Select Batch Folder');
if isequal(batchRoot,0)
    disp('User selected Cancel');
    set(hObject,'Value',0);
    return
else   
    disp(['User selected ', batchRoot])
end

pathName = batchRoot;
tifFiles = dir(fullfile(batchRoot,'*.tif'));
n = length(tifFiles);
for i = 1:n
    fileName = tifFiles(i).name;
    fullFN = fullfile(pathName, fileName);
    ConfocalStackAnalysisLM('runAll_Callback',hObject,eventdata,guidata(hObject));
end
disp('==Batch Processing Done.');


% --- Executes on button press in globalStats.
function globalStats_Callback(hObject, eventdata, handles)
% hObject    handle to globalStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fileName pathName fullFN batchModeOn batchRoot
global home globalRoot

%select globalRoot
globalRoot = uigetdir(home,'Select Root Folder');
if isequal(globalRoot,0)
    disp('User selected Cancel');
    set(hObject,'Value',0);
    return
else   
    disp(['User selected ', globalRoot])
end

pathName = globalRoot;
rootDirs = dir(fullfile(globalRoot));
n = length(rootDirs);
if ~n 
    disp('==No files found.');
    return
end

fn = [globalRoot '\GlobalStats.xls'];
if exist(fn,'file')   
    msgbox('GlobalStats.xls file exists.');
    return
end

hB = waitbar(0, 'Global Analysis in progress...');
rows0 = {{'StackName'},{'B/G Mean'},{'B/G Max'},{'B/G Sum'}};
rows1 = {{'_'},[],[],[]};
for i = 1:n
    fileName = rootDirs(i).name;
    if strcmp(fileName,'.') || strcmp(fileName,'..')
        continue
    end
    fullFN = fullfile(pathName, fileName);
    localFN = [fullFN '\' fileName '_GB_stats.xls'];
    if ~exist(localFN,'file')   
        disp(['==No Excel file found in ' fullFN '. Skip folder.']);
        continue
    end
    [newRows] = getStats(localFN);
    rows1 = vertcat(rows1,newRows);
    waitbar(i/n, hB, 'Global Analysis in progress...');    
end

% % write stats to excel file
% if ~exist(fn,'file')
%     xlswrite(fn,rows0); 
% end
% [success,message] = xlsappend(fn,rows1);
%%
xlswrite(fn,rows0); 
xlsappend(fn,rows1);

waitbar(1, hB, 'Global Analysis in progress...');   
close(hB);
msgbox('Global statistics saved to Excel.');


function [rowsData] = getStats(localFN)
global fullFN globalRoot

%% open each folder, read statistics and append to the global
    rowsData = [];
    [ndata, text, alldata] = xlsread(localFN);
    n = length(alldata(:,1));
    i = 1; k = 0;
    while i <= n
        var1 = alldata{i,1};
        var2 = alldata{i,2};
        if strcmp(var1,'_') && isnan(var2)
            k = k+1;
            stackName(k) = alldata(i+1,1);
            nrFrames = alldata{i+1,3};
            % Mean of G/B row
            stackMean(k) = alldata{i+10,nrFrames+3}; 
            % Max  of G/B row
            v = alldata(i+10,[2:nrFrames+1]);
            stackMax(k) = max( cell2mat(v)); 
            % Sum(G)/Sum(B) row
            vG = alldata(i+8,[2:nrFrames+1]);
            vB = alldata(i+9,[2:nrFrames+1]);
            sumG = sum( cell2mat(vG), 2);
            sumB = sum( cell2mat(vB), 2);
            stackSum(k) = sumG/sumB;
            %writeStat(globalFN,stackName,stackInt);
            %rows1 = '*';
            i = i + 10;
        end
        i = i + 1;
    end
rowsData = stackName';
c2 = num2cell(stackMean');
rowsData = horzcat(rowsData,c2);
c2 = num2cell(stackMax');
rowsData = horzcat(rowsData,c2);
c2 = num2cell(stackSum');
rowsData = horzcat(rowsData,c2);

    
