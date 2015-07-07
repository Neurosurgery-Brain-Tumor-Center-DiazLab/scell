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

% Last Modified by GUIDE v2.5 22-Feb-2015 20:40:41

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

if strcmp(get(hObject,'Visible'),'off')
    main_data=get(handles.norm_tool_root,'UserData');
    if isempty(main_data), main_data=struct('d',[]); end
    if length(varargin)>0
        d=varargin{1};
        main_window_handle=varargin{2};
    else, return; end
    if isfield(d,'factor_ids')%redo factor analysis every time tool is loaded
        d.factor_ids={};%string IDs for each factor in the regression model
        d.fac_varexp={};%vectors of variance explained per gene, per factor
        d.fac_counts={};%matrices of counts, used to generate each factor, stored samples-by-genes
        d.factor={};%matrices of the factors themselves, derived from d.fac_counts
    end
    main_data.d=d;
    main_data.main_window_handle=main_window_handle;
    set(handles.norm_tool_root,'UserData',main_data);
end

% --- Outputs from this function are returned to the command line.
function varargout = norm_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
gene_select_tool(d,handles.norm_tool_root);



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
    d.factor={};%matrices of the factors themselves, derived from d.fac_counts, e.g. via SVD
end
h=waitbar(0,'Processing factors');
%identify which factors to use
%ERCCs
if get(handles.ercc_checkbox,'Value')&&~any(strcmp(d.factor_ids,'ERCCs')) 
    if ~isfield(main_data,'last_dir')
        [fname pname]=uigetfile('*.*','Select ERCCs readcounts...');
    else
        [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select ERCCs readcounts...');
    end
    try
        ercc_cts=load_counts(fullfile(pname,fname),'hum','ct');
        main_data.last_dir=pname;
    catch me
        alert('String','Error loading ERCCs');
        return;
    end
    d.factor_ids{end+1}='ERCCs';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
    d.fac_counts{end+1}=ercc_cts.counts';
    d.factor{end+1}=[];
    [U,S,~]=svd(log(ercc_cts.counts'+1));
    W=U*S;
    ds=diag(S);
    p=length(ds);
    gk=cumsum(1./[p:-1:1])/p;%broken stick criterion to select singular values
    gk=gk(end:-1:1)';
    chk=diag(S)/sum(ds);
    d.factor{end+1}=W(:,gk<chk);
end
waitbar(0.25,h,'Processing factors');
%mutual background
if get(handles.bak_checkbox,'Value')&&~any(strcmp(d.factor_ids,'Background'))
    if ~isfield(d,'bak_idx')
        alert('String',sprintf('No estimate of mutual background found\nUncheck mutual background, or close the normalization tool and run QC first'));
        return;
    end
    d.factor_ids{end+1}='Background';
    d.fac_varexp{end+1}=zeros(length(d.gsymb),1); 
    d.fac_counts{end+1}=d.counts(d.bak_idx,:)';
    [U,S,~]=svd(log(d.counts(d.bak_idx,:)'+1));
    W=U*S;
    ds=diag(S);
    p=length(ds);
    gk=cumsum(1./[p:-1:1])/p;%broken stick criterion
    gk=gk(end:-1:1)';
    chk=diag(S)/sum(ds);
    d.factor{end+1}=W(:,gk<chk);
end
waitbar(0.5,h,'Processing factors');
%cyclins/CDKs
if get(handles.cyclin_checkbox,'Value')
    if ~any(strcmp(d.factor_ids,'Cyclins')) %haven't already done CCA analysis
        load cyclins.mat;
        tgidx1=[];%find the cyclin/CDKs in the list, and in chosen gene panel
        for i=1:length(cln_gns)
            t=min(find(strcmpi(cln_gns{i},d.gsymb)));
            if ~isempty(t),tgidx1=[tgidx1;t];end
        end
        if isempty(tgidx1),alert('String','error no genes found');return;end
        tgidx=intersect(tgidx1,d.gidx);%only use cyclins/CDKs in chosen gene panel
        d.cyclin_idx=tgidx;
        d.factor_ids{end+1}='Cyclins';
        d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
        d.fac_counts{end+1}=d.counts(tgidx,:)';
        d.non_cln_idx=setdiff(d.gidx,tgidx1);
        X=d.counts(d.non_cln_idx,:)';
        Y=d.counts(tgidx,:)';
        [U,V]=comp_cca(X,Y,0.05);%keep cannonical factors of cyclins at p=0.05 cutoff  
        if ~isempty(U)&&~isempty(V)
            d.cyclinU=U;d.cyclinV=V;
            %plot the top 20 most correlated Cyclins
            cln_crs=nanmean(corr(U,Y)');%correlations between cyclins and gene factors
            [cln_scrs,cln_cidx]=sort(abs(cln_crs),'descend');
            rdn=sum(nanmean(corr(X,V).^2))
            f1=figure;
            set(f1,'color','w');
            ax=gca;
            set(ax,'FontSize',18);
            bar(cln_scrs(1:min(length(cln_scrs),20)));
            set(ax,'XTick',1:min(length(cln_scrs),20));
            set(ax,'XTickLabel',d.gsymb(d.cyclin_idx(cln_cidx(1:min(length(cln_cidx),20)))));
            title(['Cyclin/CDK Tenenhaus redundancy = ' num2str(rdn*100) '%'],'FontSize',18);
            ylabel(sprintf('Mean correlation \nwith non-Cyclin/CDK gene factors'));
            rotateXLabels(ax,90);
            xlim([0 min(length(cln_scrs),20)+1]);
            %plot the top 20 most correlated genes

%             rdn=sum(mean(corr(X,V).^2));
%             f2=figure;
%             set(f2,'color','w');
%             ax=gca;
%             set(ax,'FontSize',18);
%             bar(gn_scrs(1:min(length(gn_scrs),20)));
%             set(ax,'XTick',1:min(length(gn_scrs),20));
%             set(ax,'XTickLabel',d.gsymb(d.non_cln_idx(gn_cidx(1:min(length(gn_cidx),20)))));
%             title(['Cyclin/CDK Tenenhaus redundancy = ' num2str(rdn*100) '%'],'FontSize',18);
%             ylabel(sprintf('Mean correlation \nwith Cyclins-CDK factors'));
%             rotateXLabels(ax,90);
%             xlim([0 min(length(gn_scrs),20)+1]);
            %write a ranking of genes in the factor, by mean canonical cross
            %correlation
            [fname pname]=uiputfile('cyclin-CDK_correlation_rank.tsv','Where should I save a file of the most correlated Cyclin-CDKs?');
            flag=true;
            try
                f=fopen(fullfile(pname,fname),'w');
            catch me
                alert('String',['I couln''t open ' fname]);
                flag=false;
            end
            if flag
                fprintf(f,'Cyclin-CDK\tMean_correlation\n');
                for i=1:length(cln_cidx)
                    if isnan(cln_crs(cln_cidx(i))), continue;end
                    fprintf(f,'%s\t',d.gsymb{d.cyclin_idx(cln_cidx(i))});
                    fprintf(f,'%g\n',cln_crs(cln_cidx(i)));
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
                gn_crs=nanmean(corr(X,V)');%correlations between cyclins and gene factors
                fprintf(f,'Gene\tMean_correlation_with_Cyclin-CDKs\n');
                for i=1:length(d.non_cln_idx)
                    if isnan(gn_crs(i)), continue;end
                    fprintf(f,'%s\t',d.gsymb{d.non_cln_idx(i)});
                    fprintf(f,'%g\n',gn_crs(i));
                end
                fclose(f);
            end
        end
        [U1,S,~]=svd(log(Y+1));
        W=U1*S;
        ds=diag(S);
        p=length(ds);
        gk=cumsum(1./[p:-1:1])/p;%broken stick criterion
        gk=gk(end:-1:1)';
        chk=diag(S)/sum(ds);
        d.factor{end+1}=W(:,gk<chk);
    elseif isfield(d,'cyclinU')
        U=d.cyclinU;V=d.cyclinV;
        tgidx=intersect(d.cyclin_idx,d.gidx);
        d.non_cln_idx=setdiff(d.gidx,tgidx)
        X=d.counts(d.non_cln_idx,:)';
        Y=d.counts(tgidx,:)';
        %plot the top 20 most correlated Cyclins
        cln_crs=max(corr(U,Y));%correlations between cyclins and gene factors
        [cln_scrs,cln_cidx]=sort(abs(cln_crs),'descend');
        f1=figure;
        set(f1,'color','w');
        ax=gca;
        set(ax,'FontSize',18);
        bar(cln_scrs(1:min(length(cln_scrs),20)));
        set(ax,'XTickLabel',d.gsymb(d.non_cln_idx(cln_cidx(1:min(length(cln_cidx),20)))));
        ylabel('Maximum correlation with non-Cyclin/CDK genes');
        rotateXLabels(ax,90);
        %plot the top 20 most correlated genes
        gn_crs=max(corr(X,V));%correlations between cyclins and gene factors
        [gn_scrs,gn_cidx]=sort(abs(gn_crs),'descend');
        rdn=sum(max(corr(X,V).^2));
        f2=figure;
        set(f1,'color','w');
        ax=gca;
        set(ax,'FontSize',18);
        bar(gn_scrs(1:min(length(gn_scrs),20)));
        set(ax,'XTickLabel',d.gsymb(d.cyclin_idx(gn_cidx(1:min(length(gn_cidx),20)))));
        title(['Cyclin/CDK Tenenhaus redundancy = ' num2str(rdn*100) '%'],'FontSize',18);
        ylabel('Maximum correlation with Cyclins/CDKs');
        rotateXLabels(ax,90);
    end
end
%user list
% if get(handles.user_checkbox,'Value')&&~any(strcmp(d.factor_ids,'User_list'))
%     if ~isfield(main_data,'last_dir')
%         [fname pname]=uigetfile('*.*','Select gene symbol list...');
%     else
%         [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select gene symbol list...');
%     end
%     try
%         D=importdata(fullfile(pname,fname));
%         main_data.last_dir=pname;
%     catch me
%         alert('String','Error loading gene symbol list');
%         return;
%     end
%     tgidx=[];%find user genes in the list, and in chosen gene panel
%     for i=1:length(D)
%         t=min(find(strcmp(D{i},d.gsymb)));
%         if ~isempty(t),tgidx=[tgidx;t];end
%     end
%     if isempty(tgidx),alert('String','error no genes found');return;end
%     tgidx=intersect(tgidx,d.gidx);%only use user genes in chosen gene panel
%     d.factor_ids{end+1}='User_list';
%     d.fac_varexp{end+1}=zeros(length(d.gsymb),1);
%     d.fac_counts{end+1}=d.counts(tgidx,:)';
%     X=d.counts(setdiff(d.gidx,tgidx),:)';
%     Y=d.counts(tgidx,:)';
%     [U,V]=comp_cca(X,Y,0.05);%keep cannonical factors of user genes at p=0.05 cutoff  
%     d.factor{end+1}=V;
%     %[~,rsq]=glm_reg(d.factor{end},d.counts(d.gidx,:)',log(d.sf));
%     %t=d.fac_varexp{end};
%     %t(d.gidx)=rsq;
%     %d.fac_varexp{end}=t;   
%     %write a ranking of genes in the factor, by mean canonical cross
%     %correlation
%     crs=max(corr(U,Y));
%     [~,cidx]=sort(abs(crs),'descend');
%     [fname pname]=uiputfile('user_genes_correlation_rank.tsv','Where should I save a file of the genes on your list that most correlate with other genes?');
%     flag=true;
%     try
%         f=fopen(fullfile(pname,fname),'w');
%     catch me
%         alert('String',['I couln''t open ' fname]);
%         flag=false;
%     end
%     if flag
%         fprintf(f,'User_gene\tMean_correlation\n');
%         for i=1:length(cidx)
%             fprintf(f,'%s\t',d.gsymb{d.gidx(cidx(i))});
%             fprintf(f,'%g\n',crs(cidx(i)));
%         end
%         fclose(f);
%     end
%     %write a list of top correlated genes with the factor
%     [fname pname]=uiputfile('Genes_correlated_with_user_list.tsv','Where should I save genes most correlated with your gene list?');
%     flag=true;
%     try
%         f=fopen(fullfile(pname,fname),'w');
%     catch me
%         alert('String',['I couldn''t open ' fname]);
%         flag=false;
%     end
%     if flag
%         fprintf(f,'Gene\tMean_correlation_with_user_list\tCOV\n');
%         crs=max(corr(X,V)');
%         cut=quantile(crs,.5);
%         cvs=var(X)./(mean(X).^2);%coefficient of variation
%         [scvs,cvidx]=sort(cvs,'descend');
%         cgidx=setdiff(d.gidx,tgidx);
%         for i=1:length(cvidx)
%             fprintf(f,'%s\t',d.gsymb{cgidx(cvidx(i))});
%             fprintf(f,'%g\t',crs(cvidx(i)));
%             fprintf(f,'%g\n',scvs(i));
%         end
%         fclose(f);
%     end
% end
delete(h);
[var_exp,R]=lm_reg_step(d);
for i=1:length(d.fac_varexp)
    t=d.fac_varexp{i};
    t(d.gidx)=var_exp{i};
    d.fac_varexp{i}=t;
end
d.nrmC=R';

main_data.d=d;
set(handles.norm_tool_root,'UserData',main_data);

    
% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.norm_tool_root,'UserData');
dnew=main_data.d;
main_window_data=get(main_data.main_window_handle,'UserData');
d=main_window_data.d;
d.gidx=dnew.gidx;
d.iod=dnew.iod;
d.pnz=dnew.pnz;
d.iod_fdr=dnew.iod_fdr;
d.zinf_fdr=dnew.zinf_fdr;
if isfield(dnew,'factor_ids')&&~isempty(dnew.factor_ids)
    d.factor_ids=dnew.factor_ids;%string IDs for each factor in the regression model
    d.fac_varexp=dnew.fac_varexp;%vectors of variance explained per gene, per factor
    d.fac_counts=dnew.fac_counts;%matrices of counts, used to generate each factor, stored samples-by-genes
    d.factor=dnew.factor;
    d.nrmC=dnew.nrmC;
end
if ~isfield(d,'nrmC')||isempty(d.nrmC)
    sf=sum(d.counts'+1)';
    for i=1:size(d.counts,2)
        d.nrmC(:,i)=log2((d.counts(:,i)+1)*1e6/sf(i));
    end
    d.nrmC=d.nrmC(d.gidx,find(d.cidx));
end
main_window_data.d=d;
set(main_data.main_window_handle,'UserData',main_window_data);
close(handles.norm_tool_root);
