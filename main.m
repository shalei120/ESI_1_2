

tic
fprintf('Loading J,WT, WS...\n')
% J = load('J.txt');
% 
% %% Template classify
% WT = load('TS.txt');
% WT = WT + 1.1;
%  N = size(WT,1);
% DT=diag(WT*ones(N,1));
% KT = 4;
% %% slot classify
% WS = load('SS.txt');
% WS = WS + 1.1;
% DS=diag(WS*ones(N,1));
% KS = 4;
%  save data
WT = 4 - WT;
WS = 4 - WS;
DT=diag(WT*ones(N,1));
DS=diag(WS*ones(N,1));
fprintf('Loading J,WT, WS completed! Press any key to continue...\n')
% pause;
%% initialiuze of XT, XS
fprintf('Initializing XT,XS...\n')
XT = zeros(N,KT);
for i = 1:N
    XT(i,randi(KT))=1;
end
XS = zeros(N,KS);
for i = 1:N
    XS(i,randi(KS))=1;
end
fprintf('Initializing XT,XS completed!\n')

%% Start to iterate

fprintf('Starting iteration...\n')
XT_old = XT;
XS_old = XS;
Totaliter = 1;
lambda1 = 2;%/KT;
lambda2 = 210000000;%/KT;
%% NC
% XT = SpectralClustering(WT,DT,KT);
% XS = SpectralClustering(WS,DS,KS);
%% NC + SC
while 1
    JT = 1/2/f(XS,J)/lambda1*DT*(J*J');
    WT_ast = JT' + WT + JT;
    XT = SpectralClustering(WT_ast,DT,KT);
    fprintf('\n .............WT iteration completed\n');
    
    JS = 0.5 * f(XT,J)/lambda2*(ones(N)-J*J')*DS;
    WS_ast = JS' + WS + JS ;   
  
    XS = SpectralClustering(WS_ast,DS,KS);
    fprintf('\n .............WS iteration completed\n');
    
    R1 = (XT'*XT)\XT'*XT_old;
    R2 = (XS'*XS)\XS'*XS_old;
    residual2 = norm(R1'*R1 - eye(KT))+ norm(R2'*R2 - eye(KS));
    if residual2<0.001
        break;
    end
    XT_old = XT;
    XS_old = XS;
    fprintf('\n **************************************************************\n');
    fprintf('             Total iteration %d completed | Cost = %f',Totaliter,residual2);
    fprintf('\n **************************************************************\n');
    Totaliter = Totaliter + 1;
end
fprintf('Iteration Completed, Staring to discretize...\n')
XT = discretization(XT,KT,N);
XS = discretization(XS,KS,N);
fprintf('Discretize completed, start to write to file...\n')
%% output the result(template,slot)
% save slot.mat XT XS
fp = fopen('NCSC.txt','w');
for t = 1:KT
    ct = XT(:,t) == 1;
    for i = 1:KS
        c = XS(:,i) == 1;
        ts = ct & c;
        ts = find(ts == 1);
        fprintf(fp,'%d\n',size(ts,1));
        for j=1:size(ts,1)
            fprintf(fp,'%d  ',ts(j));
        end
        fprintf(fp,'\n\n');
    end
end

%% For statistics
% fp = fopen('D:\Grade1_2\RandomWalk\NCSC-stat.txt','w');
% for e = 1 : N
%     fprintf(fp,'%d\t',find(J(e,:)==1));
%     fprintf(fp,'%d\t',find(XT(e,:)==1));
%     fprintf(fp,'%d\t',find(XS(e,:)==1));
%     fprintf(fp,'\n');
% end
% fprintf('Every thing is done!!! Congratulations!!!\n')
% toc




% %% output image
%
% I = X_ast * [0 64 128 192 256]';
% I = reshape(I,93,93);
% I = uint8(I);
% imshow(I)

%% output the result(template)
% save template.mat X_ast
% fp = fopen('templates.txt','w');
% for i = 1:K
%     c = find(X_ast(:,i) == 1);
%     fprintf(fp,'%d\n',size(c,1));
%     for j=1:size(c,1)
%         fprintf(fp,'%d  ',c(j));
%     end
%     fprintf(fp,'\n\n');
% end

