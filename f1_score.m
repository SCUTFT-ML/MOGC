function [F1]=f1_score(label,predict)
    M=confusionmat(label,predict);
    M=M';
    precision=diag(M)./(sum(M,2)+0.0001);
    recall=diag(M)./(sum(M,1)+0.0001);
    precision=mean(precision);
    recall=mean(recall);
    F1=2*precision*recall/(precision+recall);
end