function brainNetSet = MVND(BOLD,lengW,step,lab)
%  low- and high-order FC network construction based on MVND
% Input:
% BOLD signals (volumes x ROIs)
% lengW: sliding window length, number
% step: sliding window step, number
% lab: 1 is patient while -1 is normal
% Output:
% brainNetSet:include low- and high-order FC data sets. (ROIs x ROIs x subjects)
% By Yueying Zhou and Lishan Qiao
%
nSubj=length(BOLD)
nROI=size(BOLD{1},2)
nDegree=size(BOLD{1},1)%the number of sampling
lambda=[0 10 20 30 40 50 60 70 80 90 100] % the values lies in [0,100] denoting the sparsity degree
nPar=length(lambda);
brainNetSet=cell(1,nPar);
lengS=size(BOLD{1},1);%BOLD signal length
nWindow=round((lengS-lengW)/step+1)   %%number of sliding windows;
for L=1:nPar
    brainNet=zeros(nROI,nROI,nWindow,nSubj);
    for i=1:nSubj
        for w = 1:nWindow
            datasw=BOLD{i}((step*(w-1)+1:step*(w-1)+lengW),:);%slide window
            currentNet=corrcoef(datasw);
            currentNet=currentNet-diag(diag(currentNet));% no link to oneself
            threhold=prctile(abs(currentNet(:)),lambda(L)); % fractile quantile
            currentNet(find(abs(currentNet)<=threhold))=0;
            brainNet(:,:,w,i)=currentNet;
        end
        brainNetSet{L}=brainNet;
    end
    fprintf('Done %d/%d networks!\n',L,nPar);
end
% This is the progressing by sliding windows, and down is the proposed
%high-order FC estimation
for L=1:nPar
    wNet=zeros(nROI,nROI,nSubj);
    for i=1:nSubj
        W=brainNetSet{L}(:,:,:,i);
        M=sum(W,3)/nWindow;% MLE of the mean matrix
        wNet(:,:,i)=M;
        Enew=eye(nROI);
        Eold=ones(nROI);
        while sum(abs(Eold-Enew)>= 0.00001)
            Eold=Enew;
            A=zeros(nROI);
            for w=1:nWindow
                A=A+(W(:,:,w)-M)'*inv(Eold)*(W(:,:,w)-M); %MLE of the covaiance matrix
            end
            Enew=A/(nROI*nWindow)+0.01*eye(nROI);% iteration
        end
        s=diag(Enew);
        if (any (s~=1))
            Enew=Enew./sqrt(s*s'); %transfrom the covariance matrix into correlation matrix
        end 
        CNetSet(:,:,i)=Enew;
    end
    meaSet{L}=wNet;%low-order FC
    Hig{L}=CNetSet;%high-order FC
    fprintf('Done %d/%d networks!\n',L,nPar);
end
save('mean_winMCI.mat','meaSet','lab','-v7.3');
save('cov_winMCI.mat','Hig','lab','-v7.3');
fprintf('this is %d/%d networks!\n',lengW,step);





