% parameter relationship in the alogrithm
clear all;
close all;
clc;

%% config
methodname='arni';
datapath = 'data_dream4';

subdatas = 1:5;
matdatas= {
    'idxcol_bio216-workstation_frac1_bs1000_poly_20180226231842.mat',
    'idxcol_bio216-workstation_frac1_bs1000_poly_20180227000334.mat',
    'idxcol_bio216-workstation_frac1_bs1000_poly_20180227004823.mat',
    'idxcol_bio216-workstation_frac1_bs1000_poly_20180227013320.mat',
    'idxcol_bio216-workstation_frac1_bs1000_poly_20180227021758.mat'
};

for idx=1: length(subdatas)
    subdata = int2str(subdatas(idx));
    matdata = matdatas{idx};
    
    dataname = ['insilico_size100_' subdata '_multifactorial.tsv'];
    goldenname = ['DREAM4_GoldStandard_InSilico_Size100_multifactorial_'  subdata '.tsv'];
    pdffile = [pwd,'/',datapath,'/', 'pdf_multifactorial_' subdata '.mat'];
    bootstrapmat = [pwd,'/',datapath,'/',dataname, '_result/'  matdata];
    
    
    dataset = [pwd,'/',datapath,'/',dataname];
    data2=importdata(dataset,'\t');
    
    gexp_data = data2.data;
    [nsamples, ngenes] = size(gexp_data);
    
    % Set regulator indices
    nregs = ngenes;
    reg_indices = 1:nregs;
    
    result_dir = [datapath, '/', dataname, '_result'];
    if ~isdir(result_dir)
        mkdir(result_dir);
    end
    
    load(bootstrapmat);
    
    score_mat = repmat(0,5,5);
    auroc_mat = repmat(0,5,5);
    aupr_mat = repmat(0,5,5);
    
    for nbootstraps=100:100:500;
        for n_step = 1:5;
            fprintf('# Genes:%d, #tx-factor:%d, #samples: %d, #lars step:%d, #bootstraps: %d\n', ...
                ngenes,nregs,nsamples,n_step,nbootstraps);
            
            subidx_mat_collection = idx_mat_collection(:,:,1:nbootstraps);
            
            freq = get_freq_mat(subidx_mat_collection, nregs);                      %Get tx-factor selection frequnecy matrix
            
            scores = score_edges(freq,'method','area','L',n_step);
            
            name_net = [result_dir, '/','arni_prediction.txt'];
            
            predict_network(scores,reg_indices,'genenames',strrep(data2.colheaders,'"',''),'name_net',name_net);
            
            goldfile = [datapath, '/', goldenname];
            
            %% load gold standard
            gold_data = load_dream4_network(goldfile);
            test_data = load_dream4_network(name_net);
            pdf_data = load(pdffile);
            
            
            %% calculate performance metrics
            [aupr auroc prec rec tpr fpr p_auroc p_aupr] = eval_dream4(test_data, gold_data, pdf_data);
            
            score = - 0.5*(log10(p_auroc*p_aupr) );
            
            fprintf('AUROC: %f, p-val: %e, AUPR: %f, p-val: %e, score:%f\n', ...
                auroc,p_auroc,aupr,p_aupr, score);
            
            score_mat(nbootstraps/100, n_step) = score;
            auroc_mat(nbootstraps/100, n_step) = auroc;
            aupr_mat(nbootstraps/100, n_step) = aupr;
            
            
        end
    end
    
    xlswrite([result_dir,'/', 'idx_score.xls'], score_mat);
    xlswrite([result_dir,'/','idx_auroc.xls'],  auroc_mat);
    xlswrite([result_dir,'/','idx_aupr.xls'],aupr_mat);
    
    changeplot;
    t = 1:1:5;
    
    % figure;
    % y1 = auroc_mat(1,t);
    % y2 = auroc_mat(2,t);
    % y3 = auroc_mat(3,t);
    % y4 = auroc_mat(4,t);
    % y5 = auroc_mat(5,t);
    % y6 = auroc_mat(6,t);
    %
    % grid on;
    % hold on;
    % plot(t,y1);
    % plot(t,y2);
    % plot(t,y3);
    % plot(t,y4);
    % plot(t,y5);
    % plot(t,y6);
    % hold off;
    %
    % legend('BS=50','BS=100','BS=150','BS=200','BS=250','BS=300')
    % title(dataname,'Interpreter', 'none');
    % xlabel('Number of ARNI steps');
    % ylabel('AUROC');
    % saveas(gcf,[result_dir '/auroc_vs_S.png']);
    
    figure;
    y1 = score_mat(1,t);
    y2 = score_mat(2,t);
    y3 = score_mat(3,t);
    y4 = score_mat(4,t);
    y5 = score_mat(5,t);
    % y6 = score_mat(6,t);
    % y7 = score_mat(7,t);
    % y8 = score_mat(8,t);
    % y9 = score_mat(9,t);
    % y10 = score_mat(10,t);
    
    grid on;
    hold on;
    plot(t,y1);
    plot(t,y2);
    plot(t,y3);
    plot(t,y4);
    plot(t,y5);
    % plot(t,y6);
    % plot(t,y7);
    % plot(t,y8);
    % plot(t,y9);
    % plot(t,y10);
    hold off;
    
    legend('BS=100','BS=200','BS=300','BS=400','BS=500')
    %,'BS=600','BS=700','BS=800','BS=900','BS=1000'
    title(dataname,'Interpreter', 'none');
    xlabel('Number of ARNI steps');
    ylabel('Score');
    saveas(gcf,[result_dir '/score_vs_S.png']);
end

function freq = get_freq_mat(idx_mat_collection, nregs) % A helper function to compute frequency [0,1] of selection at each ARNI
% steps for all regulators.
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