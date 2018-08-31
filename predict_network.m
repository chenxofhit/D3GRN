function edges=predict_network(scores,tf_index,varargin)

% Create list of predicted edges and optionally print them in a text file
% given a matrix of scores
%
% Syntax 1: edges=predict_network(scores,tf_index) will return the entire 
% list of edges with non-zero scores. tf_index is the set of indices 
% representing tfs.
%
% Syntax 2: edges=predict_network(scores,tf_index,'option',optionvalue) 
% same with optional parameters chosen by user.
%
% REQUIRED INPUTS: 
% - scores: a ntf*ntg matrix of scores obtained running score_edges.
% - tf_index: the set of indices representing the tfs in the scores matrix.
%
% OPTIONAL INPUTS:
% - cutoff: number of edges in the predicted network (default=all possible)
% - genenames: cell array of strings containing names of genes -- required 
%              to write network in a file. No default value.
% - name_net: path/name.ext of the file to write the network -- required to
%             to write network in a file. No default value.
% Note that the function will return a warning if any of 'name_net' or
% 'genenames' is empty. 
% 
% Example: 
% load data1
% freq=tigress(data1)
% scores=score_edges(freq,'L',3)
% edges=predict_network(scores,data1.tf_index,'cutoff',10000,'genenames',da
% ta1.genenames,'name_net','./my_net1.txt')
% 
% See also: tigress.m, score_edges.m
%
% Anne-Claire Haury, 2012


%% Parse arguments
p = inputParser;   % Create an instance of the class.
p.addRequired('scores', @isfloat);
p.addRequired('tf_index',@isfloat);
p.addParamValue('cutoff', Inf, @isfloat);
p.addParamValue('genenames',{},@iscell);
p.addParamValue('name_net','',@ischar);
p.parse(scores,tf_index,varargin{:});

%% Show which arguments were not specified in the call.
disp(' ') 
for k=1:numel(p.UsingDefaults)
   field = char(p.UsingDefaults(k));
   value = p.Results.(field);
   if isempty(value)   
       value = '[]';   
   end
   fprintf('   ''%s''    defaults to %s \n', field, value)
end


%% Initialization
ngenes=size(scores,2);
ntf=length(tf_index);
edges = zeros(ngenes*ntf-ntf,3);

%% Create list of edges
k=0;
for i=1:ntf
    for j=1:ngenes
        if tf_index(i)~=j
            k=k+1;
            edges(k,:) =[tf_index(i) j scores(tf_index(i),j)];
        end
    end
end

%% Sort the edges
[a idx] = sort(edges(:,3),'descend');
edges = edges(idx,:);

%% Remove zero scores
idx=edges(:,3)>0;
edges=edges(idx,:);
nedges=size(edges,1);

%% Threshold (nothing changes if cutoff=Inf)
if nargin > 3
    nedges = min(p.Results.cutoff,nedges);
    edges=edges(1:nedges,:);
end

%% Write the network in a text file (optional)

if ~isempty(p.Results.genenames) && ~isempty(p.Results.name_net)
    genenames=p.Results.genenames;
    name_net=p.Results.name_net;
    fid = fopen(name_net,'w');
    for i=1:nedges
         fprintf(fid,'%s\t%s\t%f\r\n',genenames{edges(i,1)},genenames{edges(i,2)},edges(i,3));
    end
    fclose(fid);
end