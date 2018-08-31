function [AUROC,AUPR,fpr,recall,recall_new,precision_new]=auc_from_ranks(acc_ranked,adj,n_points)

% Created by Nan Papili Gao
%            Institute for Chemical and Bioengineering 
%            ETH Zurich
%            E-mail:  nanp@ethz.ch
%
% Copyright. November 1, 2016.

maxWEIGHT=max(max(abs(acc_ranked)));
acc_ranked=acc_ranked/maxWEIGHT;
adj=adj*3;

fp=zeros(1,n_points+1);
tn=zeros(1,n_points+1);
tp=zeros(1,n_points+1);
fn=zeros(1,n_points+1);
fpr=zeros(1,n_points+1);
precision=zeros(1,n_points+1);
recall=zeros(1,n_points+1);

thr=0:1/n_points:1;
for i=1:n_points+1
    sign_info=sign(acc_ranked);
    sign_info(sign_info==0)=1;
    adj_m=(abs(acc_ranked)>=thr(i)).*sign_info;
    
    compare=abs(adj_m+adj);
    fp(1,i)=nnz(compare<3 & compare>0);
    tp(1,i)=nnz(compare==4);
    fn(1,i)=nnz(compare==3);
    tn(1,i)=nnz(compare==0);
    precision(1,i)=tp(1,i)/(tp(1,i)+fp(1,i));
    recall(1,i)=tp(1,i)/(tp(1,i)+fn(1,i));
    fpr(1,i)=fp(1,i)/(fp(1,i)+tn(1,i));
end
i=1;
precision=fliplr(precision);
recall=fliplr(recall);
precision_new=[1 precision];
recall_new=[0 recall];
fpr=fliplr(fpr);

AUROC(i)=auc(fpr(i,:),recall(i,:));
AUPR(i)=auc(recall_new(i,:),precision_new(i,:));
end