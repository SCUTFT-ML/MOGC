function [Lambda] = solveV(V,alpha)
    [numnodes,nummotif]=size(V);
    for i = 1:1:numnodes
        v = V(i,:);
        [v1,idx1]=sort(v);
        P = 0;
        for j = nummotif:-1:1
            theta = (2*alpha+sum(v1(1:j)))/j;
            if theta - v1(j) > 0
                P = j;
                break;
            end
        end
        theta = (2*alpha+sum(v1(1:P)))/min(P,nummotif);
        lambda = (theta*ones(size(v1)) - v1)/2/alpha;
        lambda(P+1:nummotif)=0;
        Lambda(i,idx1)=lambda;
    end
end