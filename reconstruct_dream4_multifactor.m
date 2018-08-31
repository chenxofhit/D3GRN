function [list,vec]=reconstruct_dream4_multifactor(Xtf,Ygene,n_step)

% BASIS='polynomial';
ORDER=5;
th=0.000001;

disp('Initiating reconstruction...');

X=Xtf;
Yg=Ygene;

% Beginning of reconstruction algorithm
disp('Performing ARNI...');

K=ORDER;
[N,M]=size(X);

Y=zeros(K+1,M,N);

for n=1:N
    for k=0:K
        Y(k+1,:,n)=X(n,:).^(k);
    end
end

nolist=1:N;
list=[];
%cost=[];
b=1;
vec=zeros(1,N);
while (~isempty(nolist)) && (b==1);
    % Composition of inferred subspaces
    Z=[];
    for n=1:length(list)
        Z=[Z;Y(:,:,list(n))];
    end
    
    % Projection on remaining composite spaces
    P=zeros(length(nolist),2);
    %cost_err=zeros(length(nolist),1);
    for n=1:length(nolist)
        % Composition of a possible space
        R=[Z;Y(:,:,nolist(n))];
        % Error of projection on possible composite space
        prjYg = Yg*pinv(R)*R;
        P(n,1)=std(Yg-prjYg);
        P(n,2)=nolist(n);
        % Fitting cost of possible composite space
        %cost_err(n,1)=1/M *norm(Yg-prjYg);
        R=[];
    end
    
    if std(P(:,1))<th || length(list) >=n_step
        b=0;
        break
    else
        % Selection of composite space which minimizes
        % projection error
        [MIN,block]=min(P(:,1));
        list=[list,P(block,2)];
        nolist(nolist==P(block,2))=[];
        vec(1,P(block,2))=MIN;
        %cost=[cost,cost_err(block,1)];
    end
end
end