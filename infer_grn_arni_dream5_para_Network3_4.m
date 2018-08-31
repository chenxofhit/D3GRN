% To deal with Network 3 and Network4 with 100 bootstraps

%config
methodname='arni';
datapath = 'data_dream5';

dataname = 'Network3';
bootstrapmat = '/Users/chenx/Documents/Exp/bLARS-Exp/ARNI/data_dream5/Network3_result/idxcol_bio10_bs201_all_20180308193025.mat';

% dataname = 'Network4';
% bootstrapmat = '/Users/chenx/Documents/Exp/bLARS-Exp/ARNI/data_dream5/Network4_result/idxcol_bio3_bs208_all_20180308203332.mat';

% Read input data, set default arguments and validate arguments
data2 = read_data(datapath,dataname);
gexp_data = data2.expdata;

nregs = length(data2.tf_index);

% Get data stats
[nsamples, ngenes] = size(gexp_data);

% Set regulator indices
reg_indices = 1:nregs;

result_dir = [datapath, '/', dataname, '_result'];
if ~isdir(result_dir)
    mkdir(result_dir);
end

load(bootstrapmat);

score_mat = repmat(0, 2,25);
auroc_mat = repmat(0,2,25);
aupr_mat = repmat(0,2,25);

for nbootstraps=100:100:200;
    for nlars_step = 1:25;
        fprintf('# Genes:%d, #tx-factor:%d, #samples: %d, #lars step:%d, #bootstraps: %d\n', ...
            ngenes,nregs,nsamples,nlars_step,nbootstraps);
        
        subidx_mat_collection = idx_mat_collection(:,:,1:nbootstraps);
        
        % Get tx-factor selection frequnecy matrix
        freq = get_freq_mat(subidx_mat_collection, nregs);
        
       %% Score edges
        scores = score_edges(freq,'method','area','L',nlars_step);
        
       %% predicted edges
        name_net = [result_dir, '/','arni_prediction.txt'];
        
       %% Write edges
        predict_network(scores,data2.tf_index,'genenames',data2.genenames,'name_net',name_net);
        
        nw = [datapath, '/', dataname];
        goldfile = [nw, '_gold_standard.txt'];
        
        if exist(goldfile,'file')
            fprintf('Found gold-standard networks. Evaluating ... \n')
            [auroc, aupr, score] = eval_dream5_nws(datapath, dataname);
            
            score_mat(nbootstraps/100, nlars_step) = score;
            auroc_mat(nbootstraps/100, nlars_step) = auroc;
            aupr_mat(nbootstraps/100, nlars_step) = aupr;
            
        end
    end
end

xlswrite([result_dir,'/', 'idx_score.xls'], score_mat);
xlswrite([result_dir,'/','idx_auroc.xls'],  auroc_mat);
xlswrite([result_dir,'/','idx_aupr.xls'],aupr_mat);

changeplot;
t = 1:1:25;
y1 = score_mat(1,t);
y2 = score_mat(2,t);
% y3 = score_mat(3,t);
% y4 = score_mat(4,t);
% y5 = score_mat(5,t);

grid on;
hold on;
plot(t,y1);
plot(t,y2);
% plot(t,y3);
% plot(t,y4);
% plot(t,y5);
hold off;

legend('BS=100','BS=200')
title(dataname);
xlabel('Number of ARNI steps');
ylabel('Score');
saveas(gcf,[result_dir '/score_vs_S.png']);


% A helper function to compute frequency [0,1] of selection at each LARS
% steps for all tx-factors.
function freq = get_freq_mat(idx_mat_collection, nregs)
    [nlars_step, ngenes, ~] = size(idx_mat_collection);
    freq = zeros(nregs, nlars_step, ngenes);
    for i=1:ngenes
        tmp = zeros(nregs, nlars_step);
        for l=1:nlars_step
            x = tabulate(squeeze(idx_mat_collection(l,i,:)));
            tmp(x(:,1), l) = x(:,3);
        end
    freq(:,:,i) = cumsum(tmp,2)/100;
    end
end
