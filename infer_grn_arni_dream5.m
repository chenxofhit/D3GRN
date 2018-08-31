function infer_grn_arni(datapath, dataname, n_step, nbootstraps)

% Read input data, set default arguments and validate arguments
data2 = read_data(datapath,dataname);
gexp_data = data2.expdata;

nregs = length(data2.tf_index);

% Get data stats
[nsamples, ngenes] = size(gexp_data);

% Set regulator indices
reg_indices = 1:nregs;

% Set default parameters
if nargin < 3
    n_step = min(25, nregs-1);
    nbootstraps = 500;
elseif nargin < 4
    nbootstraps = 500;
end

samp_frac = 1;

% Number of LARS step should be at most one less than the number of
% transcription factors.
if n_step > nregs - 1
    error('Number of LARS step should be less than number of regulators minus one!');
end

%% Normalize gene expression data to zero mean and unit variance
gexp_data = bsxfun(@minus, gexp_data, mean(gexp_data));
std_data = std(gexp_data); std_data(std_data == 0) = eps;
gexp_data = bsxfun(@rdivide, gexp_data, std_data);
clear std_data;

% Initialize adjacency matrix for all bootstrap runs
idx_mat_collection = zeros(n_step, ngenes, nbootstraps);

fprintf('# Genes:%d, #tx-factor:%d, #samples: %d, #lars step:%d, #bootstraps: %d\n', ...
    ngenes,nregs,nsamples,n_step,nbootstraps);
fprintf('\n');

% Check and use if parallel computing is available
toolboxName = 'Parallel Computing Toolbox';
v = ver;
if any(strcmp(toolboxName, {v.Name}))
    parallel = true;
end

    if parallel
        %   if matlabpool('size'), matlabpool close; end
        %   if matlabpool('size') == 0
        %       matlabpool('open',4) ;
        %   end

        % if matlabpool('size'), matlabpool close; end
        % Reserve Matlab pool

         delete(gcp('nocreate'));

         parpool = gcp('nocreate'); % If no pool, do not create new one.
         fprintf('Running %d MATLABPOOL workers.\n', parpool.NumWorkers);
    end

result_dir = [datapath, '/', dataname, '_result'];
if ~isdir(result_dir)
    mkdir(result_dir);
end

last_saved_file = '';
basis_fns = '_all';

%% Run over all bootstrap (bs) runsgexp_data
for bs=1:nbootstraps
    % Start timer
    
    fprintf('Bootstrap run: %d\n',bs);
    % Bootstrap related settings. Current: samp_frac samples are selected
    % uniformly at random with replacement.
    sampk = floor(nsamples*samp_frac);
    gexp_data_bs = datasample(gexp_data, sampk,'Replace',true);
    
    idx_mat_bs = zeros(n_step, ngenes);
    
    % For each of targets, run parallel jobs
    %for t = 1:ngenes
        parfor t = 1:ngenes
        % Target data
        tic
        target = gexp_data_bs(:, t);
        
        %             % If target index is less than or equal to total number of
        %             % regulators i.e. target itself is a transcription factor then
        %             % exclude this tx-factor from potential regulators list. Else
        %             % all tx-factors are portential regulators
        %             if t<= nregs
        %                 reg_idx = [1:t-1, t+1:length(reg_indices)];
        %                 regs_indices_b = [1:bsize*(t-1), bsize*t+1:bsize*length(reg_indices)];
        %             else
        %                 reg_idx = 1:length(reg_indices);
        %                 regs_indices_b = 1:bsize*length(reg_indices);
        %             end
        %
        if t<=nregs
            reg_idx = [1:t-1, t+1:length(reg_indices)];
        else
            reg_idx = 1:length(reg_indices);
        end
        
        
        % Get expression of the potential regulators
        reg_data = gexp_data_bs(:, reg_idx);
        
        [list,vec] = reconstruct_dream5(reg_data',target',n_step);
        
        nz_indices = list(1:n_step)';
        
        idx_mat_bs(:,t) = reg_idx(:,nz_indices)';
        
        fprintf('\tFinished target :%d\n', t);
        
        toc
        
    end
    % Store the indices of selected tx-factors for each target for this
    % bootstrap run
    
    idx_mat_collection(:,:,bs) = idx_mat_bs;
    if rem(bs,1) == 0 || bs == nbootstraps
        [~,x] = system('hostname');
        hn = textscan(x,'%s','delimiter','.');
        dt = datestr(now, 'yyyymmddHHMMSS');
        cfile = [result_dir,'/','idxcol_', ...
            hn{1}{1},'_bs',num2str(bs),basis_fns,'_',dt,'.mat'];
        save(cfile,'idx_mat_collection');
        % Delete last run's output
        if last_saved_file
            delete(last_saved_file);
        end
        last_saved_file = cfile;
    end
    
end
    if parallel
        % After a hard long day, time for workers to rest now
        delete(parpool);
    end
end
