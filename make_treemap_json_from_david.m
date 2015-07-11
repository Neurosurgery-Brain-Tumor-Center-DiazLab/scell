function js=make_treemap_json_from_david(x)
%function js=make_treemap_json_from_david(x)
%
%IN:
%   x is a struct contining information for a JAVA string encoding a JSON
%   object to be generated
%
%js is a string in json format containing relevant data for a Jit
%icicle plot


%JSON fields: id, name, data{area,dim,color}, children

%set up the parent node x, to be a summary of the gene list

addpath('./json');
javaaddpath('./JSON-java.jar');
js=JSON.dump(x);
javarmpath('./JSON-java.jar');
js=strrep(js,'area','$area');%the $ prefix is needed by the Jit treemap code
js=strrep(js,'dim','$dim');
js=strrep(js,'color','$color');
%js=strrep(js,'[]','{}');%this is unlikely to be necessary
js=strrep(js,'"','\"');%the treemap code sets the data as a java string
js=['"' js '"'];