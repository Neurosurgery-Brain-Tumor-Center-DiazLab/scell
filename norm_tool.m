function varargout = norm_tool(varargin)
% NORM_TOOL MATLAB code for norm_tool.fig
%      NORM_TOOL, by itself, creates a new NORM_TOOL or raises the existing
%      singleton*.
%
%      H = NORM_TOOL returns the handle to a new NORM_TOOL or the handle to
%      the existing singleton*.
%
%      NORM_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NORM_TOOL.M with the given input arguments.
%
%      NORM_TOOL('Property','Value',...) creates a new NORM_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before norm_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to norm_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help norm_tool

% Last Modified by GUIDE v2.5 20-Feb-2015 14:07:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @norm_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @norm_tool_OutputFcn, ...
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


% --- Executes just before norm_tool is made visible.
function norm_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to norm_tool (see VARARGIN)

% Choose default command line output for norm_tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes norm_tool wait for user response (see UIRESUME)
% uiwait(handles.norm_tool_root);
if strcmp(get(hObject,'Visible'),'off')
    main_data=get(handles.norm_tool_root,'UserData');
    if length(varargin)>0,d=varargin{1};else,d=[];end
    if isfield(d,'factor_ids')%redo factor analysis every time tool is loaded
        d.factor_ids={};%string IDs for each factor in the regression model
        d.fac_varexp={};%vectors of variance explained per gene, per factor
        d.fac_counts={};%matrices of counts, used to generate each factor, stored samples-by-genes
        d.factor={};%matrices of the factors themselves, derived from d.fac_counts
    end
    main_data.d=d;
    set(handles.norm_tool_root,'UserData',main_data);
end


% --- Outputs from this function are returned to the command line.
function varargout = norm_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
main_data=get(handles.norm_tool_root,'UserData');
varargout{1} = main_data.d;


% --- Executes on button press in ercc_checkbox.
function ercc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ercc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ercc_checkbox


% --- Executes on button press in bak_checkbox.
function bak_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to bak_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bak_checkbox


% --- Executes on button press in cyclin_checkbox.
function cyclin_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cyclin_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cyclin_checkbox


% --- Executes on button press in user_checkbox.
function user_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to user_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of user_checkbox


% --- Executes on button press in ercc_pushbutton.
function ercc_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ercc_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname pname]=uigetfile([],'Please load the ERCC gene identifiers.');
try
    D=importdata(fullfile(pname,fname));
catch me
    alert('String','File input error');
end
%find the ids in d and generate an index to them
main_data=get(handles.norm_tool_root,'UserData');
main_data.ercc_idx=[];
not_found={};
d=main_data.d;
for i=1:length(D)
    t=min(find(strcmp(D{i},d.gsymb)));
    if ~isempty(t), main_data.ercc_idx=[main_data.ercc_idx;t];
    else, not_found{end+1}=D{i};end
end
set(handles.norm_tool_root,'UserData',main_data);
if ~isempty(not_found)
    out_str=[num2str(length(not_found)),...
            sprintf(' genes could not be identified\n'),...
            sprintf('Genes not found:\n')];
    for u=1:length(not_found)-1,out_str=[out_str,sprintf('%s,',not_found{u})];end
    out_str=[out_str,sprintf('%s',not_found{end})];
    alert('String',out_str);
end 

% --- Executes on button press in bak_pushbutton.
function bak_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to bak_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cyclin_pushbutton.
function cyclin_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cyclin_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in user_pushbutton.
function user_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to user_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in select_genes_pushbutton.
function select_genes_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_genes_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.norm_tool_root,'UserData');
d=main_data.d;
dnew=gene_select_tool(d);
d.gidx=dnew.gidx;
d.iod=dnew.iod;
d.pnz=dnew.pnz;
d.iod_fdr=dnew.iod_fdr;
d.zinf_fdr=dnew.zinf_fdr;
main_data.d=d;
set(handles.norm_tool_root,'UserData',main_data);


% --- Executes on button press in norm_button.
function norm_button_Callback(hObject, eventdata, handles)
% hObject    handle to norm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%for the ercc's and mutual background, represent regressor with reduced
%singular-vector space of raw counts
%for cyclin/cdks and user list, use cannonical factors of raw counts
%obtained through cannonical correlation analysis with the raw counts

main_data=get(handles.norm_tool_root,'UserData');
d=main_data.d;
if ~isfield(d,'gidx')
    alert('String','select genes for analysis first');
    return;
end
if ~isfield(d,'factor_ids')
    d.factor_ids={};%string IDs for each factor in the regression model
    d.fac_varexp={};%vectors of variance explained per gene, per factor
    d.fac_counts={};%matrices of counts, used to generate each factor, stored samples-by-genes
    d.factor={};%matrices of the factors themselves, derived from d.fac_counts
end
%identify which factors to use
%ERCCs
if get(handles.ercc_checkbox,'Value')&&~any(strcmp(d.factor_ids,'ERCCs')) 
    if ~isfield(main_data,'last_dir')
        [fname pname]=uigetfile('*.*','Select ERCCs readcounts...');
    else
        [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select ERCCs readcounts...');
    end
    try
        ercc_cts=load_fc_out(fullfile(pname,fname),'hum');
        main_data.last_dir=pname;
    catch me
        alert('String','Error loading ERCCs');
        return;
    end
    d.factor_ids{end+1}='ERCCs';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
    d.fac_counts{end+1}=ercc_cts.counts';
    d.factor{end+1}=[];
    [U,S,~]=svd(ercc_cts.counts');
    W=U*S;
    ds=diag(S);
    p=length(ds);
    gk=cumsum(1./[p:-1:1])/p;%broken stick criterion to select singular values
    gk=gk(end:-1:1)';
    chk=diag(S)/sum(ds);
    kpt=find(gk<chk);
    d.factor{end+1}=W(:,kpt);
    [~,rsq]=glm_reg(d.factor{end},d.counts(d.gidx,:),log(d.sf));
    t=d.fac_varexp{end};
    t(d.gidx)=rsq;
    d.fac_varexp{end}=t;
end
%mutual background
if get(handles.bak_checkbox,'Value')&&~any(strcmp(d.factor_ids,'Background'))
    if ~isfield(d,'bak_idx')
        alert('String',sprintf('No estimate of mutual background found\nUncheck mutual background, or close the normalization tool and run QC first'));
        return;
    end
    d.factor_ids{end+1}='Background';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1); 
    d.fac_counts{end+1}=d.counts(d.bak_idx,:)';
    [U,S,~]=svd(d.counts(d.bak_idx,:)');
    W=U*S;
    ds=diag(S);
    p=length(ds);
    gk=cumsum(1./[p:-1:1])/p;%broken stick criterion
    gk=gk(end:-1:1)';
    chk=diag(S)/sum(ds);
    kpt=find(gk<chk);
    d.factor{end+1}=W(:,kpt);
    [~,rsq]=glm_reg(d.factor{end},d.counts(d.gidx,:),log(d.sf));
    t=d.fac_varexp{end};
    t(d.gidx)=rsq;
    d.fac_varexp{end}=t;
end
%cyclins/CDKs
if get(handles.cyclin_checkbox,'Value')&&~any(strcmp(d.factor_ids,'Cyclins'))
    if ~isfield(main_data,'last_dir')
        [fname pname]=uigetfile('*.*','Select Cyclin/CDKs gene symbol list...');
    else
        [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select Cyclin/CDKs gene symbol list...');
    end
    try
        D=importdata(fullfile(pname,fname));
        main_data.last_dir=pname;
    catch me
        alert('String','Error loading Cyclins/CDKs gene symbol list');
        return;
    end
    tgidx=[];%find the cyclin/CDKs in the list, and in chosen gene panel
    for i=1:length(D)
        t=min(find(strcmp(D{i},d.gsymb)));
        if ~isempty(t),tgidx=[tgidx;t];end
    end
    if isempty(tgidx),alert('String','error no genes found');return;end
    tgidx=intersect(tgidx,d.gidx);%only use cyclins/CDKs in chosen gene panel
    d.factor_ids{end+1}='Cyclins';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
    d.fac_counts{end+1}=d.counts(tgidx,:)';
    X=d.counts(setdiff(d.gidx,tgidx),:)';
    Y=d.counts(tgidx,:)';
    [U,V]=comp_cca(X,Y,0.05);%keep cannonical factors of cyclins at p=0.05 cutoff  
    d.factor{end+1}=V;
    [~,rsq]=glm_reg(d.factor{end},d.counts(d.gidx,:),log(d.sf));
    t=d.fac_varexp{end};
    t(d.gidx)=rsq;
    d.fac_varexp{end}=t;   
    %write a ranking of genes in the factor, by mean canonical cross
    %correlation
    crs=mean(corr(U,Y));
    [~,cidx]=sort(abs(crs),'descend');
    [fname pname]=uiputfile('cyclin-CDK_correlation_rank.tsv','Where should I save a file of the most correlated Cyclin-CDKs?');
    flag=true;
    try
        f=fopen(fullfile(pname,fname),'w');
    catch me
        alert('String',['I couln''t open ' fname]);
        flag=false;
    end
    if flag
        fprintf(f,[fstr '_gene\tMean_correlation\n']);
        for i=1:length(cidx)
            fprintf(f,'%s\t',d.gsymb{d.gidx(cidx(i))});
            fprintf(f,'%g\n',crs(cidx(i)));
        end
        fclose(f);
    end
    %write a list of top correlated genes with the factor
    [fname pname]=uiputfile('Cyclin-CDK_correlated_genes.tsv','Where should I save genes most correlated with Cyclins/CDKs?');
    flag=true;
    try
        f=fopen(fullfile(pname,fname),'w');
    catch me
        alert('String',['I couldn''t open ' fname]);
        flag=false;
    end
    if flag
        fprintf(f,'Gene\tMean_correlation_with_Cyclin-CDKs\tCOV\n');
        crs=max(corr(X,V)');
        cut=quantile(crs,.5);
        cvs=var(X)./(mean(X).^2);%coefficient of variation
        [scvs,cvidx]=sort(cvs,'descend');
        cgidx=setdiff(d.gidx,tgidx);
        for i=1:length(cvidx)
            fprintf(f,'%s\t',d.gsymb{cgidx(cvidx(i))});
            fprintf(f,'%g\t',crs(cvidx(i)));
            fprintf(f,'%g\n',scvs(i));
        end
        fclose(f);
    end
end
%user list
if get(handles.user_checkbox,'Value')&&~any(strcmp(d.factor_ids,'User_list'))
    if ~isfield(main_data,'last_dir')
        [fname pname]=uigetfile('*.*','Select gene symbol list...');
    else
        [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select gene symbol list...');
    end
    try
        D=importdata(fullfile(pname,fname));
        main_data.last_dir=pname;
    catch me
        alert('String','Error loading gene symbol list');
        return;
    end
    tgidx=[];%find user genes in the list, and in chosen gene panel
    for i=1:length(D)
        t=min(find(strcmp(D{i},d.gsymb)));
        if ~isempty(t),tgidx=[tgidx;t];end
    end
    if isempty(tgidx),alert('String','error no genes found');return;end
    tgidx=intersect(tgidx,d.gidx);%only use user genes in chosen gene panel
    d.factor_ids{end+1}='User_list';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
    d.fac_counts{end+1}=d.counts(tgidx,:)';
    X=d.counts(setdiff(d.gidx,tgidx),:)';
    Y=d.counts(tgidx,:)';
    [U,V]=comp_cca(X,Y,0.05);%keep cannonical factors of user genes at p=0.05 cutoff  
    d.factor{end+1}=V;
    [~,rsq]=glm_reg(d.factor{end},d.counts(d.gidx,:),log(d.sf));
    t=d.fac_varexp{end};
    t(d.gidx)=rsq;
    d.fac_varexp{end}=t;   
    %write a ranking of genes in the factor, by mean canonical cross
    %correlation
    crs=mean(corr(U,Y));
    [~,cidx]=sort(abs(crs),'descend');
    [fname pname]=uiputfile('user_genes_correlation_rank.tsv','Where should I save a file of the genes on your list that most correlate with other genes?');
    flag=true;
    try
        f=fopen(fullfile(pname,fname),'w');
    catch me
        alert('String',['I couln''t open ' fname]);
        flag=false;
    end
    if flag
        fprintf(f,[fstr '_gene\tMean_correlation\n']);
        for i=1:length(cidx)
            fprintf(f,'%s\t',d.gsymb{d.gidx(cidx(i))});
            fprintf(f,'%g\n',crs(cidx(i)));
        end
        fclose(f);
    end
    %write a list of top correlated genes with the factor
    [fname pname]=uiputfile('Genes_correlated_with_user_list.tsv','Where should I save genes most correlated with your gene list?');
    flag=true;
    try
        f=fopen(fullfile(pname,fname),'w');
    catch me
        alert('String',['I couldn''t open ' fname]);
        flag=false;
    end
    if flag
        fprintf(f,'Gene\tMean_correlation_with_user_list\tCOV\n');
        crs=max(corr(X,V)');
        cut=quantile(crs,.5);
        cvs=var(X)./(mean(X).^2);%coefficient of variation
        [scvs,cvidx]=sort(cvs,'descend');
        cgidx=setdiff(d.gidx,tgidx);
        for i=1:length(cvidx)
            fprintf(f,'%s\t',d.gsymb{cgidx(cvidx(i))});
            fprintf(f,'%g\t',crs(cvidx(i)));
            fprintf(f,'%g\n',scvs(i));
        end
        fclose(f);
    end
end
keyboard()
    
