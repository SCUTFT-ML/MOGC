clear
clc
seed = 1;
rng(seed); 
% load the data
data = {'football','polbooks','polblogs','Cora','email','Citeseer','Pubmed'};
motif = {'m_3^3'};
for dataset = 1:1:length(data)
    path = '~/data/';
    adj_path = strcat(path,strcat('/',data{dataset},'.mat'));
    load(adj_path);
    label = gnd;
    K=length(unique(gnd));
    uA = double(A);
    [LCC, lcc_inds, ci, sizes] = LargestConnectedComponent(sparse(uA)); 
    ty_lcc=label;
    ty_lcc = ty_lcc(lcc_inds,1);
    
    %   code for  edge & M_3^3
    A1=full(LCC);
    W1=eye(size(A1,1));
    D=sum(A1);
    idx=find(D==0);
    W1(idx,idx)=0;
    
    % generate the indicator matrix W1 of the isolated nodes in the original graph

    A={};
    W={};
    A{1}=A1;
    W{1}=W1;
    
    % generate adjacency matrix for other motifs
    for motifid = 1:1:length(motif)
        if motif{motifid} == "m_3^3"
            A2 = MotifAdjacency(A1, 'm4');
            A2=full(A2);
            W2=eye(size(A2,1));
            D=sum(A2);
            idx=find(D==0);
            W2(idx,idx)=0;
            A{motifid+1}=A2;
            W{motifid+1}=W2;
        end
        if motif{motifid} == "m_3^2"
            A2 = MotifAdjacency(A1, 'm13');
            A2=full(A2);
            W2=eye(size(A2,1));
            D=sum(A2);
            idx=find(D==0);
            W2(idx,idx)=0;
            A{motifid+1}=A2;
            W{motifid+1}=W2;
        end
        if motif{motifid} == "4motif"
            load(strcat(path,'/adj_4motif.mat'));
            A2=full(A_4motif);
            W2=eye(size(A2,1));
            D=sum(A2);
            idx=find(D==0);
            W2(idx,idx)=0;
            A{motifid+1}=A2;
            W{motifid+1}=W2;
        end
        if motif{motifid} == "5motif"
            load(strcat(path,'/adj_5motif.mat'));
            A2=full(A_5motif);
            W2=eye(size(A2,1));
            D=sum(A2);
            idx=find(D==0);
            W2(idx,idx)=0;
            A{motifid+1}=A2;
            W{motifid+1}=W2;
        end
    end
    
    % alpha can be tuned
    alpha = 0.9;
    [U,lambda]=Clusterwithisolated(A,W,alpha,K);

    nmitotal=0;
    RItotal=0;
    % repeat 20 times
    for i=1:1:20
        idx=litekmeans(U, K, 'Replicates', 200);
        nmi=NMI(ty_lcc,idx);
        [~,F1,RI,~,~]=RandIndex(ty_lcc,idx);
        RItotal=RItotal+RI;
        nmitotal = nmitotal + nmi;
    end
    
    fprintf('RandIndex and NMI of MOGC for %s are %6.4f\t, %6.4f\n', strjoin(motif), RItotal/20, nmitotal/20);

end