function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 25-Nov-2014 16:22:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%set(jCBList, 'ValueChangedCallback', @myMatlabCallbackFcn);
% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.main_window);



% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function cell_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%columnname={'Include','ID','type','Reads','Genes tagged','COV','Simpson diversity','Lorenz dynamic range','Turing coverage'};
%columnformat={'logical','char','char','numeric','numeric','numeric','numeric','numeric','numeric'};
%set(hObject,'ColumnName',columnname,'ColumnFormat',columnformat,'ColumnEditable',[true false false false false false false false false]);
set(hObject,'Data',{});


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uigetfile('*.*','Select FeatureCounts output...');
else
    [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select FeatureCounts output...');
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.main_window,'UserData',main_data);
    h=waitbar(.5,['Loading ' strrep(fname,'_','\_')]);
    d=load_fc_out(fullfile(pname,fname),'hum');
    d.qc=false;
    main_data.d=d;
end
set(handles.main_window,'UserData',main_data);
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    ct_dat{i,1}=true;
    ct_dat{i,2}=d.slbls{i};
    ct_dat{i,3}=sum(d.counts(:,i));
    ct_dat{i,4}=nnz(d.counts(:,i));
end
set(handles.cell_table,'Data',ct_dat);
delete(h)



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in restore_button.
function restore_button_Callback(hObject, eventdata, handles)
% hObject    handle to restore_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in qc_button.
function qc_button_Callback(hObject, eventdata, handles)
% hObject    handle to qc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
M=zeros(size(d.counts,2),4);%matrix to hold qc metrics, used for pareto rank
preseq_dir='./';
%for each cell call preseq
j=1;
h=waitbar(0,['processing cell ' num2str(j) ' of ' num2str(size(d.counts,2))]);
while ~d.qc&&j<=size(d.counts,2)
    fname=tempname;
    f=fopen(fname,'w');
    for i=1:size(d.counts,1)
        if d.counts(i,j)>0,fprintf(f,'%d\n',d.counts(i,j));end
    end
    [status,result]=system([preseq_dir 'preseq lc_extrap -V ' fname]);
    if status~=0, d.preseq(j)=0;
    else
        D=textscan(result,'%n%n%n%n','Headerlines',1);
        D=D{2};
        d.preseq(j)=nnz(d.counts(:,j))/median(D(floor(length(D)*.75):end));
    end
    M(j,1)=d.preseq(j);
    fclose(f);
    %turing
    d.turing(j)=1-sum(d.counts(:,j)==1)/sum(d.counts(d.counts(:,j)>0,j));
    t=find(d.counts(:,j)>0);
    pt=d.counts(t,j)/sum(d.counts(t,j));
    M(j,2)=d.turing(j);
    %simpson
    d.simpson(j)=1/(pt'*pt);
    M(j,3)=d.simpson(j);
    waitbar(j/size(d.counts,2),h,['processing cell ' num2str(j) ' of ' num2str(size(d.counts,2))]);
    j=j+1;
end
%compute lorenz outlier detection
[~,sf,~,lorenzh,pval,sidx,cxi]=normalize_samples(d.counts,[],1);
delete(h);
figure
set(gcf,'color','w');
subplot(2,1,1);
set(gca,'FontSize',18);
ylabel('Simpson diversity','FontSize',18)
notBoxPlot(d.simpson(pval<0.05),1);
hold
notBoxPlot(d.simpson(pval>=0.05),2);
set(gca,'XTick',1:2,'XTickLabel',{'QC pass','QC fail'});
subplot(2,1,2);
set(gca,'FontSize',18);
ylabel('Preseq coverage','FontSize',18)
notBoxPlot(d.preseq(pval<0.05),1);
hold
notBoxPlot(d.preseq(pval>=0.05),2);
set(gca,'XTick',1:2,'XTickLabel',{'QC pass','QC fail'});
d.lorenz=pval;d.lorenzh=lorenzh;d.sf=sf;
M(:,4)=d.lorenz;
%compute pareto ranking
[~,f]=paretofronts(M,[1 1 1 1]);
d.pareto=f;
d.sidx=sidx;%index vector which sorts order stats of geo-mean
d.cxi=cxi;%index scalar of point of maximal separation of geomean to poisson noise
d.qc=true;
main_data.d=d;
set(handles.main_window,'UserData',main_data);
%update sample list
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    ct_dat{i,5}=d.simpson(i);
    ct_dat{i,6}=d.preseq(i);
    ct_dat{i,7}=d.turing(i);
    ct_dat{i,8}=d.lorenz(i);
    ct_dat{i,9}=d.pareto(i);
end
set(handles.cell_table,'Data',ct_dat);


% --- Executes on button press in norm_button.
function norm_button_Callback(hObject, eventdata, handles)
% hObject    handle to norm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in analyze_button.
function analyze_button_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
pca_tool(main_data.d);

% --- Executes during object creation, after setting all properties.
function cell_info_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_info_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected cell(s) is changed in cell_table.
function cell_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to cell_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
mdat=get(handles.main_window,'UserData');
d=mdat.d;
tdat=get(hObject,'Data');
t=find(strcmp(d.slbls,tdat{eventdata.Indices(1),2}));%second column is sample name
m=get(handles.disp_table,'Data');
m{1}=d.slbls{t};
m{2}=sum(d.counts(:,t));
m{3}=nnz(d.counts(:,t));
if isfield(d,'turing')&&~isempty(d.turing)%then we've done QC
    m{4}=d.simpson(t);
    m{5}=d.preseq(t);
    m{6}=d.turing(t);
    m{7}=d.lorenz(t);
    m{8}=d.pareto(t);
else
    for i=4:8, m{i}='NA'; end
end
set(handles.disp_table,'Data',m);



% --- Executes during object creation, after setting all properties.
function disp_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to disp_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Data',cell(8,1));
