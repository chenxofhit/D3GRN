function scores = score_edges(freq,varargin)

% Transforms frequencies into scores
% 
% Syntax 1: scores = score_edges(freq)
% Returns a matrix of scores given a frequency matrix freq. All options are
% set to default values
%
% Syntax 2: scores = score_edges(freq,'option',option_value,...)
% Same with optional settings set by user. See below.
% 
% REQUIRED INPUT:
% - freq: the frequency matrix obtained by running tigress.
% 
% OPTIONAL INPUTS: 
% - method: either 'area' (default) or 'original' (see paper)
% - L: number of LARS steps to consider (default=max_possible)
% 
% OUTPUT:
% - scores: a ntf*ntg matrix of scores (i.e. probabilities)
% 
% Example:
% load data1
% freq=tigress(data1)
% scores=score_edges(data1,'method','area','L',5)
% 
% See also: tigress, predict_network
%
% Anne-Claire Haury, 2012 

%% Parse arguments
p = inputParser;   % Create an instance of the class.
p.addRequired('freq', @isfloat);
p.addParamValue('method','area',@(x)any(strcmpi(x,{'area','original'})));
[n nsteps ng]=size(freq);
p.addParamValue('L', nsteps, @(x)x<nsteps+1);
p.parse(freq,varargin{:})

%% Show which arguments were not specified in the call.
disp(' ') 
for k=1:numel(p.UsingDefaults)
   field = char(p.UsingDefaults(k));
   value = num2str(p.Results.(field));
   if isempty(value)   
       value = '[]';   
   end
   fprintf('   ''%s''    defaults to %s \n', field, value)
end

%% Restrain freq to the first steps we want to look at
freq=freq(:,1:p.Results.L,:);

%% Now apply any of 'original','area', or 'warea' method (can add new ones!)
switch p.Results.method
    case 'original'
        scores=reshape(max(freq,[],2),n,ng);
    case 'area' % triangles approx. of the area
        scores=reshape((zeros(n,1,ng)+freq(:,1,:))/2+sum((freq(:,2:p.Results.L,:)+freq(:,1:p.Results.L-1,:))/2,2),n,ng);
        scores=1/(p.Results.L-.5)*scores;
end