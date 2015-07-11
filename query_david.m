function c=query_david(ids,eml)
%function c=query_david(ids,eml)
%
%IN:ids is a vector of ENTREZ gene ids 
%   eml is a string containing the email address to authenticate to use the
%   DAVID server
%
%OUT x is an object containing a david gene ontology term
%clustering report


%warning('off','all')
h=waitbar(0.1,'Initializing query...');
stub=sample.session.client.stub.DAVIDWebServiceStub('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/');
tmp=javaMethod('_getServiceClient',stub);
tmp.getOptions.setManageSession(true);
waitbar(0.3,h,'authenticating...')
usr_ver=stub.authenticate(eml);
if ~strcmpi(char(usr_ver),'true')
    delete(h);
    alert('Title','Authentication error','String',{['failed to authenticate ' eml],'Make sure you are registered at http://david.abcc.ncifcrf.gov/webservice/register.htm'});
    return;
else
    waitbar(0.6,h,'submitting gene list...')
    %if strcmpi(gtype,'ENTREZ_GENE_ID')
        lst=[];
        for i=1:length(ids)-1,lst=[lst num2str(ids(i)) ','];end
        lst=[lst num2str(ids(end))];
    %else
    %    lst=[];
    %    for i=1:length(ids)-1,lst=[lst ids{i} ','];end
    %    lst=[lst ids{end}];
    %end
    lst_name='user_list';
    lst_type=0;
    vld=stub.setCategories('');
    addlst_out=stub.addList(lst,'ENTREZ_GENE_ID',lst_name,lst_type);
    %stub.getSpecies
    %all_lst_name=stub.getAllListNames();
    overlap = 3;
    initialSeed = 2;
    finalSeed = 2;
    linkage = 0.5;
    kappa = 20;
    waitbar(0.9,h,'requesting the GO term cluster report...')
    c=stub.getTermClusterReport(overlap,initialSeed,finalSeed,linkage,kappa);
    delete(h);
end
%for i=1:length(d),javarmpath(fullfile(['david_java_client/lib/' d(i).name]));end