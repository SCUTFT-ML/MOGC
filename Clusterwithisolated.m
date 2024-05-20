% MOGC - Multi-order Graph Clustering with Adaptive Node-level Weight Learning
% 
% *** Description ***
%
% See Ye Liu, Xuelei Lin, Yejia Chen, and Reynold Cheng, Multi-order Graph Clustering with Adaptive 
% Node-level Weight Learning,  Submitted.
%
% ****** Input ******
%       Ablk: adjacency matrix of motif-based graph
%       Wblk: indicator matrix of isolated nodes for each graph
%       alpha: regularization parameter for lambda term
%       K: number of Cluster
%       normalize: whether to normalize the eigenvector
%
% ****** Output ******
%       U: indicator matrix of cluster
%       lambda: weight matrix of each motif to each node
%
% 

function [U,lambda]=Clusterwithisolated(Ablk,Wblk,alpha,K,normalize)

    nummotif=size(Ablk,2);
    numnodes=size(Ablk{1,1},1);
    
%     initialization
    lambda=ones(numnodes,nummotif)/nummotif;
    lambda = lambda./repmat(sum(lambda,2),1,nummotif);
    Atilde=[];
    lambdatilde=[];
    eta = 1;
    for i = 1:1:nummotif
        Atmp=Ablk{1,i}*Wblk{1,i};
        Atilde=[Atilde;Atmp];
        lambdatilde=[lambdatilde,diag(lambda(:,i))];
    end
    
%     Initial adjacency matrix 
    A=0.5*(lambdatilde*Atilde+Atilde'*lambdatilde');
    
% generalized eigen-decomposition
    eig_opts = {};
    eig_opts.tol = 1e-6;
    eig_opts.isreal = true;
    eig_opts.issym = true;
    D = sum(A,2);
    L = diag(D)-A;
%     %         solve U
    [U,V] = eig(L,diag(D));
%     [U, ~] = eigs(L,K, 'sa', eig_opts); would be considered to accelerate
    U=U(:,1:K);
    if normalize
        T = spdiag(1 ./ sqrt(sum(abs(U).^2, 2))); 
        U = T * U;
    end
        
    lambda0=zeros(numnodes,nummotif);
    U0 = zeros(numnodes,K);
    error=min(norm(lambda0-lambda,'fro'),norm(U-U0,'fro'));
    M=[];
    for i=1:1:nummotif
        M=[M,eye(numnodes)];
    end
    iter = 1;
    MM = M*M';
    MMinverse = (M*M')^(-1);
    for m = 1:1:nummotif
        WA{1,m} = Wblk{1,m} * Ablk{1,m};
        diagAW{1,m} = diag(Ablk{1,m}*Wblk{1,m}*ones(numnodes,1));
        WAdiag{1,m} = WA{1,m} *diag(ones(numnodes,1));
    end
    for i = 1:1:nummotif
        nodeidx{i,1} = find(diag(Wblk{1,i})==0);
    end
    while(error>10^(-4)&&iter<100)
        U0 = U;
        
%         solve lambda
        tildeA=zeros(1,numnodes*nummotif);
        for i = 1:1:K
            tildeAtmp=[];
            for m=1:1:nummotif
                tildeAtmp=[tildeAtmp,U(:,i)'*WA{1,m}*diag(U(:,i))];
            end
            tildeA=tildeA+tildeAtmp;
        end
        
        hatA=zeros(1,numnodes*nummotif);
        P1 = [];
        for i = 1:1:K
            hatAtmp=[];
            for m=1:1:nummotif
                hatAtmp=[hatAtmp,U(:,i)'*diagAW{1,m}*diag(U(:,i))];
            end
            hatA=hatA+hatAtmp;
            P1 = [P1;hatAtmp];
        end
        
        barA=zeros(1,numnodes*nummotif);
        P2 = [];
        for i = 1:1:K
            barAtmp=[];
            for m=1:1:nummotif
                barAtmp=[barAtmp,U(:,i)'*diag(U(:,i))*WAdiag{1,m}];
            end
            barA=barA+barAtmp;
            P2 = [P2;barAtmp];
        end
        
        V=(hatA+barA)*0.5-tildeA;
        
        P = P1 + P2;
%         调用cgs得到phi1, phi2
        PPinverse = (P*P')^-1;
        Acgs1 = MM-M*P'*PPinverse*P*M';
        bcgs1 = 2*alpha*ones(numnodes,1)+M*V'-M*P'*PPinverse*(2*alpha*ones(K,1)+P*V');
        phi1 = cgs(Acgs1,bcgs1);
        Acgs2 = P*P'-P*M'*MMinverse*M*P';
        bcgs2 = 2*alpha*ones(K,1)+P*V'-P*M'*MMinverse*(2*alpha*ones(numnodes,1)+M*V');
        phi2 = cgs(Acgs2, bcgs2);
        lambdatilde = (M'* phi1 + P'*phi2 - V')/2/alpha;
        lambdatilde(find(lambdatilde<0))=0;
        
%         reshape lambda
        lambdatilde=reshape(lambdatilde,[numnodes,nummotif]);
        for i = 1:1:nummotif
            lambdatilde(nodeidx{i,1},i)=0;
        end
        lambdatilde = lambdatilde./repmat(sum(lambdatilde,2),1,nummotif);
        lambdatilde(find(isnan(lambdatilde)))=0;
        lambda0 = lambda;
        lambda=lambdatilde;
        lambdatilde=[];
        for i = 1:1:nummotif
            lambdatilde=[lambdatilde,diag(lambda(:,i))];
        end
        
%       update fused adjacency matrix A
        A=0.5*(lambdatilde*Atilde+Atilde'*lambdatilde');
%       generalized eigen-decomposition
        eig_opts = {};
        eig_opts.tol = 1e-6;
        eig_opts.isreal = true;
        eig_opts.issym = true;
        D = sum(A,2);
        L = diag(D)-A;
        %         solve U
        [U,V] = eig(L,diag(D));
        U=U(:,1:K);
        if normalize
            T = spdiag(1 ./ sqrt(sum(abs(U).^2, 2))); 
            U = T * U;
        end

        error=max(norm(lambda0-lambda,'fro'),norm(U-U0,'fro'));
        iter = iter+1;
    end
end