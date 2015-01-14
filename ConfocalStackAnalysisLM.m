function varargout = ConfocalStackAnalysisLM(varargin)
%% Lilia Mesina, Polaris/CCBN, Sptember2014
% Last Modified by LM, 1Dec2014

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
global ver
global step1 step2 step3 
global runAllStepsOn batchModeOn
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip



handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%set(hObject, 'KeyPressFcn',@fKeyPress);
movegui(hObject,'northwest');

%inint vars
ver = '2.3';
modeMultip = [];
redMinTresh   = [];     redMaxTresh   = [];
greenMinTresh = [];     greenMaxTresh = [];
blueMinTresh  = [];     blueMaxTresh  = [];

runAllStepsOn = 0;
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
global fileName pathName fullFN batchModeOn runAllStepsOn
global elems statFr zxy_stack num_images
global newRedS newGreenS newBlueS 

global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip
global ver


sel = get(hObject,'Value');
if sel
    if ~batchModeOn || ~runAllStepsOn
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
    caption = sprintf('ConfocalStackAnalysisLMv%s >>%s', ver, fullFN);
    set(gcf, 'Name', caption );

    %% load stack
    info = imfinfo(fullFN);
    [img_size1 img_size2 img_size3]=size(imread(fullFN, 1, 'Info', info));
    img_type=['uint' num2str(info(1).BitDepth)];
    num_images = numel(info);
    elems  = img_size1 * img_size2; %1024.^2;
    statFr    = zeros(13,num_images); 
    zxy_stack = zeros([num_images img_size1 img_size2 img_size3],'uint8'); %, img_type);
    newRedS   = zeros([num_images img_size1 img_size2]);
    newGreenS = zeros([num_images img_size1 img_size2]);
    newBlueS  = zeros([num_images img_size1 img_size2]);

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
    if isempty(modeMultip)
        set_tresh; 
    end
    activeFrame = 1;
    step1 = 1;
    draw_frame(hObject, eventdata, handles);
    set(hObject,'Value',1);
    h = msgbox('Load confocal stack Done.');
    pause(1); close(h);
end


function runAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to runAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global statFr zxy_stack num_images
global step1 step2
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip
global newRedS newGreenS newBlueS


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
    stdDev = nanstd(statFr,0,2);
    statFr(:,num_images+1:end) = [];
    statFr = horzcat(statFr,stdDev);
    rowMean = nanmean(statFr(:,1:num_images),2);
    statFr = horzcat(statFr,rowMean);
    set(handles.uitable1,'data',statFr);    
    set(hObject,'Value',1);
    step2 = 1;
    waitbar(1, hB, 'Image Analysis in progress...');   
    close(hB);
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ver
global step2 statFr fullFN num_images
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip


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
    rows1{2,4}  = 'redMinTresh';     rows1{2,5}  = redMinTresh;
    rows1{2,6}  = 'redMaxTresh';     rows1{2,7}  = redMaxTresh;
    rows1{2,8}  = 'greenMinTresh';   rows1{2,9}  = greenMinTresh;
    rows1{2,10} = 'greenMaxTresh';   rows1{2,11} = greenMaxTresh;
    rows1{2,12} = 'blueMinTresh';    rows1{2,13} = blueMinTresh;
    rows1{2,14} = 'blueMaxTresh';    rows1{2,15} = blueMaxTresh;
    rows1{2,16} = 'modeMultip';      rows1{2,17} = modeMultip;
    rows1{2,18} = 'AnalysisDate';    rows1{2,19} = datestr(now);
    rows1{2,20} = 'GUI Version';     rows1{2,21} = ver;
    rows2 = {'MinRed';'MaxRed';'MinGreen';'MaxGreen';...
             'MinBlue';'MaxBlue';'MeanRed';'MeanGreen';...
             'RedPixels';'GreenPixels';'BluePixels';...
             'R/B';'G/B'};
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
hFig = findobj('Tag','newRedFig');
if ~isempty(hFig)
    ConfocalStackAnalysisLM('ShowNewRed_Callback',hObject,eventdata,guidata(hObject));
end
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

r = fr(:,:,1); g = fr(:,:,2); b = fr(:,:,3);

%% get image profile
rM = reshape(r,[1,elems]); rM = double(rM);
gM = reshape(g,[1,elems]); gM = double(gM);
bM = reshape(b,[1,elems]); bM = double(bM);

%% get stats
statFr(1,k) = min(rM);
statFr(2,k) = max(rM);
statFr(3,k) = min(gM);
statFr(4,k) = max(gM);
statFr(5,k) = min(bM);
statFr(6,k) = max(bM);


function getStats2(fr,k)
global statFr
global newRedS newGreenS newBlueS
global redMinTresh greenMinTresh
global blueMinTresh blueMaxTresh

r = fr(:,:,1); g = fr(:,:,2); b = fr(:,:,3);

% bM = reshape(b,[1,elems]);
% %m0(m0==0) = [];
% gM = reshape(g,[1,elems]);
% %m0(m0==0) = [];

%% get newRed
 [n1, n2] = size(r(:,:));
 newRed = r;
 for i = 1:n1
    for j = 1:n2
        if (b(i,j) < blueMinTresh) || (r(i,j) <= redMinTresh) %(blue(i,j) < 20.4) || (green(i,j) <= 12)
            newRed(i,j) = 0;
        end
    end
 end
 
 %image(newGreen(:,:));
 red2 = find(newRed>0);
 nR = length(red2);
 nR2 = nR; %/elems;
 statFr(9,k) = nR2;

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
 statFr(10,k) = nG2;

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
statFr(11,k) = nB2;
statFr(12,k) = nR2/nB2;
statFr(13,k) = nG2/nB2;
newRedS(k,:,:)   = newRed;
newGreenS(k,:,:) = newGreen;
newBlueS(k,:,:)  = newBlue;
%getStats1(fr,k);


function set_tresh()
global elems num_images statFr zxy_stack
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip

modeMultip = 1.2;       
hObj = findobj('Tag', 'modeM');     set(hObj, 'String', num2str(modeMultip));

rMinTresh = zeros(num_images,1);
gMinTresh = zeros(num_images,1);
bMinTresh = zeros(num_images,1);
hHist = figure();
for i = 1:num_images
    fr = zxy_stack(i,:,:,:);
    fr = squeeze(fr);
    %fr3 = imagesc(fr);
    r = fr(:,:,1);    g = fr(:,:,2);     b = fr(:,:,3); %imagesc(
%%     %% clean Glia
% figure
% h2 =reshape(b,[1,elems]);
% h2 = double(h2);
% a1 = bar(h2); %,64);
% plot(a1);
% %bar(log(a1));
% figure
% bar(a1);
% shg
%%
    h2r = reshape(r,[1,elems]);  statFr(7,i) = mean(h2r);
    h2g = reshape(g,[1,elems]);  statFr(8,i) = mean(h2g); %? make it double
    h2b = reshape(b,[1,elems]);
    
    h3 = imhist(r);
    figure(hHist);
    bar(h3); pause(0.5);
    rMinTresh(i) = modeMultip * mode(h2r);
    
    h3 = imhist(g);
    figure(hHist);
    bar(h3); pause(0.5);
    gMinTresh(i) = modeMultip * mode(h2g);    

    h3 = imhist(b);
    figure(hHist);
    bar(h3); pause(0.5);
    bMinTresh(i) = modeMultip * mode(h2b);
end
close(hHist);
if isempty(redMinTresh), redMinTresh = mean(rMinTresh); end %min 12      
hObj = findobj('Tag', 'rMinTresh'); set(hObj, 'String', num2str(redMinTresh));
redMaxTresh = 255;      
hObj = findobj('Tag', 'rMaxTresh'); set(hObj, 'String', num2str(redMaxTresh));
if isempty(greenMinTresh), greenMinTresh = mean(gMinTresh); end %min 12
hObj = findobj('Tag', 'gMinTresh'); set(hObj, 'String', num2str(greenMinTresh));
greenMaxTresh = 255;    
hObj = findobj('Tag', 'gMaxTresh'); set(hObj, 'String', num2str(greenMaxTresh));
if isempty(blueMinTresh), blueMinTresh = mean(bMinTresh); end %min 30 to skip glia
hObj = findobj('Tag', 'bMinTresh'); set(hObj, 'String', num2str(blueMinTresh));
blueMaxTresh = 230;     
hObj = findobj('Tag', 'bMaxTresh'); set(hObj, 'String', num2str(blueMaxTresh));


hObj = findobj('Tag', 'rMinTresh'); set(hObj, 'String', num2str(redMinTresh));
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


% --- Executes on button press in ShowNewRed.
function ShowNewRed_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNewRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowNewRed
global newRedS activeFrame statFr

sel = get(hObject,'Value');
if sel
    if ~mean(statFr(7,:))
        hM = msgbox('First run Step 2: Image Analysis.');
        set(hObject,'Value',0);
        return
    end
    hFig = findobj('Tag','newRedFig');
    if isempty(hFig)
        hFig = figure('Tag','newRedFig','Name','NewRed','NumberTitle','off');
    end
    figure(hFig);
    caption = ['New Red: Frame ' num2str(activeFrame)];
    set(hFig,'Name',caption);
    im = newRedS(activeFrame,:,:);
    im = squeeze(im);
    imagesc(im);
end


% --- Executes on button press in ShowNewGreen.
function ShowNewGreen_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNewGreen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newGreenS activeFrame statFr

sel = get(hObject,'Value');
if sel
    if ~mean(statFr(8,:))
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
    if ~mean(statFr(8,:))
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
global greenMaxTresh

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
global blueMinTresh 

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
global blueMaxTresh

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
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip
global runAllStepsOn

runAllStepsOn = 1;
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
global redMinTresh redMaxTresh 
global greenMinTresh greenMaxTresh 
global blueMinTresh blueMaxTresh modeMultip


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
global home globalRoot pathName

%if use fullFN, might need to update it first
%select globalRoot
globalRoot = uigetdir(home,'Select Root Folder');
if isequal(globalRoot,0)
    disp('User selected Cancel');
    set(hObject,'Value',0);
    return
else   
    disp(['User selected ', globalRoot])
end

fn = [globalRoot '\GlobalStats.xls'];
if exist(fn,'file')   
    msgbox('GlobalStats.xls file already exists.');
    return
end


pathName = globalRoot;
rootDirs = subdir(fullfile(globalRoot, '*ConfocalStats.xls'));

n = length(rootDirs);
if ~n 
    disp('==No Excel files found.');
    return
end

hB = waitbar(0, 'Global Analysis in progress...');

rows0 = {'StackName','R/B Mean','R/B Max','R/B Sum', ...
                     'G/B Mean','G/B Max','G/B Sum'};
rows1 = [];
for i = 1:n
    localFN = rootDirs(i).name;
%     fullFN = fullfile(pathName, fileName);
%     localFN = [fullFN '\' fileName '_ConfocalStats.xls'];
    [newRows] = getStats(localFN);
    rows1 = vertcat(rows1,newRows);
    waitbar(i/n, hB, 'Global Analysis in progress...');    
end

%%
xlswrite(fn,rows0); 
xlsappend(fn,rows1);

waitbar(1, hB, 'Global Analysis in progress...');   
close(hB);
msgbox('Global statistics saved to Excel.');


function [rowData] = getStats(localFN)

%%  open each file, read statistics and append to the global
    rowData = [];
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
            
            %% Mean of R/B row
            stackMeanR(k) = alldata{i+14,nrFrames+3}; 
            % Max  of R/B row
            v = alldata(i+14,[2:nrFrames+1]);
            stackMaxR(k) = max( cell2mat(v)); 
            % Sum(R)/Sum(B) row
            vR = alldata(i+11,[2:nrFrames+1]);
            vB = alldata(i+13,[2:nrFrames+1]);
            sumR = sum( cell2mat(vR), 2);
            sumB = sum( cell2mat(vB), 2);
            stackSumR(k) = sumR/sumB;
            
            %% Mean of G/B row
            stackMeanG(k) = alldata{i+15,nrFrames+3}; 
            % Max  of G/B row
            v = alldata(i+15,[2:nrFrames+1]);
            stackMaxG(k) = max( cell2mat(v)); 
            % Sum(G)/Sum(B) row
            vG = alldata(i+12,[2:nrFrames+1]);
            sumG = sum( cell2mat(vG), 2);
            stackSumG(k) = sumG/sumB;
            %writeStat(globalFN,stackName,stackInt);
            %rows1 = '*';
            i = i + 15;
        end
        i = i + 1;
    end
rowData = stackName';
%% R/B
c2 = num2cell(stackMeanR');
rowData = horzcat(rowData,c2);
c2 = num2cell(stackMaxR');
rowData = horzcat(rowData,c2);
c2 = num2cell(stackSumR');
rowData = horzcat(rowData,c2);
%% G/B
c2 = num2cell(stackMeanG');
rowData = horzcat(rowData,c2);
c2 = num2cell(stackMaxG');
rowData = horzcat(rowData,c2);
c2 = num2cell(stackSumG');
rowData = horzcat(rowData,c2);



function rMinTresh_Callback(hObject, eventdata, handles)
% hObject    handle to rMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rMinTresh as text
%        str2double(get(hObject,'String')) returns contents of rMinTresh as a double
global redMinTresh redMaxTresh

val = get(handles.rMinTresh,'string');
redMinTresh = str2double(val); 
disp(['redMinTresh ' val]);

% --- Executes during object creation, after setting all properties.
function rMinTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rMinTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rMaxTresh_Callback(hObject, eventdata, handles)
% hObject    handle to rMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rMaxTresh as text
%        str2double(get(hObject,'String')) returns contents of rMaxTresh as a double
global redMaxTresh

val = get(handles.rMaxTresh,'string');
redMaxTresh = str2double(val); 
disp(['redMaxTresh ' val]);

% --- Executes during object creation, after setting all properties.
function rMaxTresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rMaxTresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%{
% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%}
