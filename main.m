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

% Last Modified by GUIDE v2.5 01-Sep-2015 14:05:31

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
if isempty(gcp('nocreate')), parpool; end
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
if isempty(main_data), main_data=struct('d',[]); end
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
    d_new=load_counts(fullfile(pname,fname),'hum','ct');
%     load name_chg.mat
%     for i=1:length(d.slbls)
%         t=find(strcmp(d.slbls{i},name_chg.old));
%         if isempty(t), keyboard(), end
%         d.slbls{i}=name_chg.new{t};
%     end
    d_old=main_data.d;
    if isempty(d_old)
        d=d_new;
        d.qc=false;
        d.cidx=ones(size(d.slbls));%index of cells to include in analysis
        d.gidx=ones(size(d.gsymb));%index of genes to include in analysis
        d.ctype=cell(size(d.slbls));
        for i=1:length(d.ctype), d.ctype{i}='NA'; end
        d.mapped=-1*ones(size(d.slbls));%mapped reads
        d.unmapped=-1*ones(size(d.slbls));%unmapped reads
        d.ld_call=cell(size(d.slbls));%live-dead call string
        for i=1:length(d.slbls), d.ngns(i)=nnz(d.counts(:,i)); end
    else
        d=d_old;
        keyboard()
        for i=1:length(d_new.slbls), d.slbls{end+1}=d_new.slbls{i}; end
        d.counts=[d.counts,d_new.counts];
        d.cpm=[d.cpm,d_new.cpm];
        d.qc=false;
        d.cidx=ones(size(d.slbls));%index of cells to include in analysis
        d.gidx=ones(size(d.gsymb));%index of genes to include in analysis
        for i=1:length(d_new.slbls), d.ctype{end+1}='NA'; end
        d.mapped=[d.mapped,ones(size(d_new.slbls))];
        d.unmapped=[d.unmapped,ones(size(d_new.slbls))];
        for i=1:length(d_new.slbls), d.ld_call{end+1}=''; end
        for i=1:length(d_new.slbls), d.ngns(end+1)=nnz(d_new.counts(:,i)); end
    end
    main_data.d=d;
end
set(handles.main_window,'UserData',main_data);
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    ct_dat{i,1}=true;
    ct_dat{i,2}=d.slbls{i};
    ct_dat{i,3}=nansum(d.counts(:,i));
    ct_dat{i,7}=nnz(d.counts(:,i));
end
set(handles.cell_table,'Data',ct_dat);
delete(h)



% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
[fname pname]=uiputfile('*.mat','Save data as...','my_session');
save(fullfile(pname,fname),'d','-v7.3');

% --- Executes on button press in restore_button.
function restore_button_Callback(hObject, eventdata, handles)
% hObject    handle to restore_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.main_window,'UserData');
if isempty(main_data), main_data=struct('d',[]); end
[fname pname]=uigetfile('*.mat','Restore session from matfile...');
load(fullfile(pname,fname));
main_data.d=d;
set(handles.main_window,'UserData',main_data);
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    if d.cidx(i), ct_dat{i,1}=true; else, ct_dat{i,1}=false; end
    ct_dat{i,2}=d.slbls{i};
    ct_dat{i,3}=nansum(d.counts(:,i));
    if isfield(d,'mapped')&&~isempty(d.mapped),ct_dat{i,4}=d.mapped(i);end
    if isfield(d,'unmapped')&&~isempty(d.unmapped),ct_dat{i,5}=d.mapped(i);end
    if isfield(d,'ld_call')&&~isempty(d.ld_call),ct_dat{i,6}=d.ld_call{i};end
    ct_dat{i,7}=nnz(d.counts(:,i));
    if isfield(d,'simpson')&&~isempty(d.simpson),ct_dat{i,8}=d.simpson(i);end
    if isfield(d,'preseq')&&~isempty(d.preseq),ct_dat{i,9}=d.preseq(i);end
    if isfield(d,'turing')&&~isempty(d.turing),ct_dat{i,10}=d.turing(i);end
    if isfield(d,'lorenz')&&~isempty(d.lorenz),ct_dat{i,11}=d.lorenz(i);end
    if isfield(d,'preseq_mar')&&~isempty(d.preseq_mar),ct_dat{i,12}=d.preseq_mar(i);end
end
set(handles.cell_table,'Data',ct_dat);
m=get(handles.disp_table,'Data');
m{1}=d.slbls{1};
m{2}=nansum(d.counts(:,1));
m{3}=nnz(d.counts(:,1));
if isfield(d,'simpson')&&~isempty(d.simpson),m{4}=d.simpson(1);end
if isfield(d,'preseq')&&~isempty(d.preseq),m{5}=d.preseq(1);end
if isfield(d,'turing')&&~isempty(d.turing),m{6}=d.turing(1);end
if isfield(d,'lorenz')&&~isempty(d.lorenz),m{7}=d.lorenz(1);end
if isfield(d,'preseq_mar')&&~isempty(d.preseq_mar),m{8}=d.preseq_mar(1);end
set(handles.disp_table,'Data',m);


% --- Executes on button press in qc_button.
function qc_button_Callback(hObject, eventdata, handles)
% hObject    handle to qc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
j=1;
h=waitbar(0,['processing cell ' num2str(j) ' of ' num2str(size(d.counts,2))]);
while ~d.qc&&j<=size(d.counts,2)
    %don't recompute metrics for samples already done
    comp_turing=~isfield(d,'turing')||isempty(d.turing)||j>length(d.turing)||d.turing(j)==0;
    comp_simp=~isfield(d,'simpson')||isempty(d.simpson)||j>length(d.simpson)||d.simpson(j)==0;
    %turing
    if comp_turing
        d.turing(j)=1-nansum(d.counts(:,j)==1)/nansum(d.counts(d.counts(:,j)>0,j));
        t=find(d.counts(:,j)>0);
        pt=d.counts(t,j)/nansum(d.counts(t,j));
    end
    %simpson
    if comp_simp
        pt=d.counts(t,j)/nansum(d.counts(t,j));
        d.simpson(j)=1/(pt'*pt);
    end
    waitbar(j/size(d.counts,2),h,['processing cell ' num2str(j) ' of ' num2str(size(d.counts,2))]);
    j=j+1;
end
delete(h);
comp_preseq=~isfield(d,'preseq')||isempty(d.preseq);
if comp_preseq
    if isdeployed
        if ~isfield(main_data,'last_dir')
            [fname, pname]=uigetfile('*.*','Select PRESEQ executable...');
        else
            [fname, pname]=uigetfile([main_data.last_dir,'*.*'],'Select PRESEQ executable...');
        end
    else
        pname='./';
        fname='preseq';
    end
    if ~isempty(fname)
        h=waitbar(0.5,'estimating library complexity');
        preseq=zeros(size(d.counts,2),1);
        preseq_mar=zeros(size(d.counts,2),1);
        counts=d.counts;
        parfor i=1:length(preseq)
            tn=tempname;
            f=fopen(tn,'w');
            for k=1:size(d.counts,1)
                if counts(k,i)>0,fprintf(f,'%d\n',counts(k,i));end
            end
            [status,result]=system([fullfile(pname,fname) ' lc_extrap -V ' tn]);
            if status==0
                D=textscan(result,'%n%n%n%n','Headerlines',1);
                if ~isempty(D)&&~isempty(D{2})
                    D=D{2};
                    if numel(D)>2
                        preseq_mar(i)=D(3)-D(2);
                        preseq(i)=nnz(counts(:,i))/median(D(floor(length(D)*.75):end));
                    end
                end
            end
            fclose(f);
        end
        waitbar(1,h,'done');
        d.preseq=preseq;
        d.preseq_mar=preseq_mar;
        delete(h);
    end
end
h=waitbar(0.5,'identifying outliers');
%compute lorenz outlier detection
[~,sf,~,lorenzh,pval,sidx,cxi,bak_idx]=normalize_samples(d.counts,[],1,d.slbls);
waitbar(1,h,'done');
[~,q]=mafdr(pval);
d.lorenz=q;d.lorenzh=lorenzh;d.sf=sf;d.bak_idx=bak_idx;
delete(h);
figure
set(gcf,'color','w');
subplot(2,1,1);
set(gca,'FontSize',18);
ylabel('Simpson diversity','FontSize',18)
H=notBoxPlot(d.simpson(d.lorenz>=0.05),1);
set([H.data],'markersize',2);
hold
H=notBoxPlot(d.simpson(d.lorenz<0.05),2);
set([H.data],'markersize',2);
set(gca,'XTick',1:2,'XTickLabel',{'QC pass','QC fail'});
subplot(2,1,2);
set(gca,'FontSize',18);
ylabel('Preseq coverage','FontSize',18)
H=notBoxPlot(d.preseq(d.lorenz>=0.05),1);
set([H.data],'markersize',2);
hold
H=notBoxPlot(d.preseq(d.lorenz<0.05),2);
set([H.data],'markersize',2);
set(gca,'XTick',1:2,'XTickLabel',{'QC pass','QC fail'});
d.sidx=sidx;%index vector which sorts order stats of geo-mean
d.cxi=cxi;%index scalar of point of maximal separation of geomean to poisson noise
d.qc=true;
main_data.d=d;
set(handles.main_window,'UserData',main_data);
%update sample list
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    ct_dat{i,8}=d.simpson(i);
    ct_dat{i,9}=d.preseq(i);
    ct_dat{i,10}=d.turing(i);
    ct_dat{i,11}=d.lorenz(i);
    ct_dat{i,12}=d.preseq_mar(i);
end
set(handles.cell_table,'Data',ct_dat);
m=get(handles.disp_table,'Data');
m{1}=d.slbls{1};
m{2}=nansum(d.counts(:,1));
m{3}=nnz(d.counts(:,1));
m{4}=d.simpson(1);
m{5}=d.preseq(1);
m{6}=d.turing(1);
m{7}=d.lorenz(1);
m{8}=d.preseq_mar(1);
set(handles.disp_table,'Data',m);

% --- Executes on button press in norm_button.
function norm_button_Callback(hObject, eventdata, handles)
% hObject    handle to norm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
if isempty(main_data)
    alert('String','Load data first');
    return;
end
d=main_data.d;
kidx=find(d.cidx);
d.slbls=d.slbls(kidx);
d.counts=d.counts(:,kidx);
if ~isfield(d,'sf')||isempty(d.sf)
    alert('String','QC libraries first')
    return;
end
d.sf=d.sf(kidx);
norm_tool(d,handles.main_window);

% --- Executes on button press in analyze_button.
function analyze_button_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
if isempty(main_data)
    alert('String','Load data first');
    return;
end
d=main_data.d;
% sc=sum(d.counts+1)';
% for i=1:size(d.counts,2)
%     d.nrmC(:,i)=log2(1e6*(d.counts(:,i)+1)/sc(i));
% end
%d.nrmC=d.nrmC(d.gidx,find(d.cidx));
cidx=find(d.cidx);
d.gsymb_full=d.gsymb;
d.cpm_full=d.cpm;
d.cidx_full=d.cidx;
d.slbls_full=d.slbls;
d.gsymb=d.gsymb(d.gidx);
d.counts=d.counts(d.gidx,cidx);
d.cpm=d.cpm(d.gidx,cidx);
d.ctype=d.ctype(cidx);
d.ngns=d.ngns(cidx);
d.ent=d.ent(d.gidx);
d.ld_call=d.ld_call(cidx);
d.lorenz=d.lorenz(cidx);
d.lorenzh=d.lorenzh(cidx);
d.mapped=d.mapped(cidx);
d.unmapped=d.unmapped(cidx);
d.preseq=d.preseq(cidx);
d.preseq_mar=d.preseq_mar(cidx);
d.sf=d.sf(cidx);
d.simpson=d.simpson(cidx);
d.slbls=d.slbls(cidx);
d.turing=d.turing(cidx);
d.iod=d.iod(d.gidx);
d.iod_fdr=d.iod_fdr(d.gidx);
d.zinf_fdr=d.zinf_fdr(d.gidx);
d.sidx=d.sidx(cidx);
d.pnz=d.pnz(d.gidx);
d.bak_idx=intersect(d.gidx,d.bak_idx);
d.cidx=ones(length(cidx),1);
d.gidx=ones(length(d.gidx),1);
h=waitbar(0.5,'Computing PCA...');
if isempty(d.nrmC)
    alert('String','Select genes/normalize samples first!');
    delete(h);
    return;
end
[coeff,score]=pca(d.nrmC');
delete(h);
c=PcaComputeMock2(coeff, score);
p= PcaTool(c);
p.show;
c.changeD(d);

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
m{2}=nansum(d.counts(:,t));
m{3}=nnz(d.counts(:,t));
if isfield(d,'turing')&&~isempty(d.turing)%then we've done QC
    m{4}=d.simpson(t);
    m{5}=d.preseq(t);
    m{6}=d.turing(t);
    m{7}=d.lorenz(t);
    m{8}=d.preseq_mar(t);
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


% --- Executes when entered data in editable cell(s) in cell_table.
function cell_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to cell_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if isempty(eventdata.NewData),return;end
main_data=get(handles.main_window,'UserData');
d=main_data.d;
hdrs=get(handles.cell_table,'ColumnName');
if strcmp(hdrs{eventdata.Indices(2)},'Include') %update whether to include cell 
    d.cidx(eventdata.Indices(1))=eventdata.NewData;
end
if strcmp(hdrs{eventdata.Indices(2)},'ID') %change cell ID
    d.slbls{eventdata.Indices(1)}=eventdata.NewData;
end
if strcmp(hdrs{eventdata.Indices(2)},'type') %change cell type
    d.ctype{eventdata.Indices(1)}=eventdata.NewData;
end
main_data.d=d;
set(handles.main_window,'UserData',main_data);


% --- Executes on button press in featureCounts_toggle.
function featureCounts_toggle_Callback(hObject, eventdata, handles)
% hObject    handle to featureCounts_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of featureCounts_toggle


% --- Executes on button press in write_qc_button.
function write_qc_button_Callback(hObject, eventdata, handles)
% hObject    handle to write_qc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
if isempty(d.preseq)
    alert('String','Run library QC first');
    return;
end
[fname pname]=uiputfile('*.txt','Save data as...','QC_metrics.txt');
try
    f=fopen(fullfile(pname,fname),'w');
catch me
    alert('String',['Error opening ' fullfile(pname,fname)]);
    return;
end
if ~isfield(d,'mapped')
    fprintf(f,'ID\tTag_count\tGenes_tagged\tSimpson_diversity\tPRESEQ_coverage\tTuring_coverage\tLorenz_outlier_test\tOutlier_PASS-FAIL\tMarginal_return\n');
    for i=1:length(d.slbls)
        fprintf(f,'%s\t',d.slbls{i});
        fprintf(f,'%i\t',nansum(d.counts(:,i)));
        fprintf(f,'%i\t',nnz(d.counts(:,i)));
        fprintf(f,'%g\t',d.simpson(i));
        fprintf(f,'%g\t',d.preseq(i));
        fprintf(f,'%g\t',d.turing(i));
        fprintf(f,'%g\t',d.lorenz(i));
        if d.lorenz(i)<0.05,fprintf(f,'FAIL\t');else,fprintf(f,'PASS\t');end
        fprintf(f,'%g\n',d.preseq_mar(i));
    end
    fclose(f);
else
    fprintf(f,'ID\tReads_mapped\tReads_unmapped\tMapping_rate\tLive-dead_call\tTag_count\tGenes_tagged\tSimpson_diversity\tPRESEQ_coverage\tTuring_coverage\tLorenz_outlier_test\tOutlier_PASS-FAIL\tMarginal_return\n');
    for i=1:length(d.slbls)
        fprintf(f,'%s\t',d.slbls{i});
        fprintf(f,'%i\t',d.mapped(i));
        fprintf(f,'%i\t',d.unmapped(i));
        fprintf(f,'%g\t',d.mapped(i)/(d.mapped(i)+d.unmapped(i)));
        fprintf(f,'%s\t',d.ld_call{i});
        fprintf(f,'%i\t',nansum(d.counts(:,i)));
        fprintf(f,'%i\t',nnz(d.counts(:,i)));
        fprintf(f,'%g\t',d.simpson(i));
        fprintf(f,'%g\t',d.preseq(i));
        fprintf(f,'%g\t',d.turing(i));
        fprintf(f,'%g\t',d.lorenz(i));
        if d.lorenz(i)<0.05,fprintf(f,'FAIL\t');else,fprintf(f,'PASS\t');end
        fprintf(f,'%g\n',d.preseq_mar(i));
    end
    fclose(f);
end


% --- Executes on selection change in filter_menu.
function filter_menu_Callback(hObject, eventdata, handles)
% hObject    handle to filter_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filter_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filter_menu


% --- Executes during object creation, after setting all properties.
function filter_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String',{'Lorenz','PRESEQ'});
    set(hObject,'Value',1);
end



function filter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_edit as text
%        str2double(get(hObject,'String')) returns contents of filter_edit as a double


% --- Executes during object creation, after setting all properties.
function filter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filter_button.
function filter_button_Callback(hObject, eventdata, handles)
% hObject    handle to filter_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.main_window,'UserData');
d=main_data.d;
thr=str2num(get(handles.filter_edit,'String'));
if isempty(thr), alert('title','Error...','String','Enter a numeric threshold'); end
ct_dat=get(handles.cell_table,'Data');
fltr_str=get(handles.filter_menu,'String');
fltr_str=fltr_str{get(handles.filter_menu,'Value')};
switch fltr_str
    case 'Lorenz'
        x=d.lorenz; %threshold on Lorenz p-value
        f=@gt; %threshold less than thr
    case 'PRESEQ'
end
for i=1:size(ct_dat,1)
    t=find(strcmp(ct_dat{i,2},d.slbls));
    if f(x(t),thr)
        ct_dat{i,1}=true;
        d.cidx(t)=1;
    else
        ct_dat{i,1}=false;
        d.cidx(t)=0;
    end
end
main_data.d=d;
set(handles.main_window,'UserData',main_data);
set(handles.cell_table,'Data',ct_dat);
        


% --- Executes on button press in load_meta_button.
function load_meta_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_meta_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
if ~isfield(main_data,'last_dir')
    [fname, pname]=uigetfile('*.*','Select meta-data file...');
else
    [fname, pname]=uigetfile([main_data.last_dir,'*.*'],'Select meta-data file...');
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.main_window,'UserData',main_data);
    h=waitbar(.5,['Loading ' strrep(fname,'_','\_')]);
    f=fopen(fullfile(pname,fname));
    try
        D=textscan(f,'%s%n%n%s%s','HeaderLines',1,'Delimiter',char(9));
    catch me
        alert('String','File I/O error, check format.')
    end
    not_found={};
    for i=1:length(D{1})
        t=find(strcmp(D{1}{i},d.slbls));
        if isempty(t)
            not_found{end+1}=D{1}{i};
            continue;
        end
        if length(t)>1, keyboard(),end
        d.mapped(t)=D{2}(i);
        d.unmapped(t)=D{3}(i);
        d.ld_call{t}=D{4}{i};
        d.ctype{t}=D{5}{i};
    end 
    main_data.d=d;
end
set(handles.main_window,'UserData',main_data);
ct_dat=get(handles.cell_table,'Data');
for i=1:length(d.slbls)
    ct_dat{i,4}=d.mapped(i);
    ct_dat{i,5}=d.unmapped(i);
    ct_dat{i,6}=d.ld_call{i};
end
set(handles.cell_table,'Data',ct_dat);
delete(h)


% --- Executes on button press in select_all.
function select_all_Callback(hObject, eventdata, handles)
% hObject    handle to select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
ct_dat=get(handles.cell_table,'Data');
for i=1:size(ct_dat,1)
    ct_dat{i,1}=true;
    d.cidx(i)=1;
end
main_data.d=d;
set(handles.main_window,'UserData',main_data);
set(handles.cell_table,'Data',ct_dat);



% --- Executes on button press in unselect_all.
function unselect_all_Callback(hObject, eventdata, handles)
% hObject    handle to unselect_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.main_window,'UserData');
d=main_data.d;
ct_dat=get(handles.cell_table,'Data');
for i=1:size(ct_dat,1)
    ct_dat{i,1}=false;
    d.cidx(i)=0;
end
main_data.d=d;
set(handles.main_window,'UserData',main_data);
set(handles.cell_table,'Data',ct_dat);
