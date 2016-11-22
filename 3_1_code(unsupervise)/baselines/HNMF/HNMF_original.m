function [  RD,F,Precision,Recall,jaccard,rat,Z_filter,pathway_ncbi ] = HNMF(g_p_4,g_p_5,M,mgi_id,standard,standard_row,K,MaxIter,InnerMaxIter,t_circle,T)

%因为第5层是父层，hnmf的效果应该是和NMF是一样的。

V1 = g_p_4;
V2 = g_p_5;
[m,n1]=size(V1);
[~,n2]=size(V2);

alphas=[0];
gamas=[0];
%lamtas1=[0.01,0.1,1,10,20];
%lamtas2=[0.01,0.1,1,10,20];
lamtas1=[10];
lamtas2=[10];

for lamta1=lamtas1
    for lamta2=lamtas2
        t=1;
        num=0;
        while 1
            num=num+1;
            W= rand(m,K);
            H1 = rand(K,n1);
            H2 = rand(K,n2);
                     
            [W_out,H_out,V_out,rmse,isCoverage] = HNMF_sparse(MaxIter,V1,V2,W,H1,H2,M,lamta1,lamta2);
            directory='../4_result/HMF';
                          if(~exist(directory,'dir'))
                              mkdir(directory);
                          end
            if(isCoverage)
                fn = ['../4_result/HMF/HNMF_result_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2)  '_t&' num2str(t)  '.mat'];
                disp([datestr(now) ':  '  fn ]);
                save(fn,'W_out','H_out','rmse');
            end
            t=t+1;
            if(t>=t_circle+1 || num>50)
                break;
            end
        end
    end
end

a = length(alphas);
b = length(gamas);
c=length(lamtas1);
d = length(lamtas2);
jaccards=zeros(c*d,t_circle);
RD = zeros(a*b*c*d,t_circle);
F = zeros(a*b*c*d,t_circle);
Z_filter = cell(a*b*c*d,t_circle);
pathway_ncbi = cell(a*b*c*d,t_circle);
Precision=zeros(a*b*c*d,t_circle);
Recall=zeros(a*b*c*d,t_circle);
jaccard=zeros(a*b*c*d,t_circle);
alpha_ratio=1;

T=3;

for j=1:c
    lamta1=lamtas1(j);
    for k=1:d
        lamta2=lamtas2(k);
        for t=1:t_circle
            fn = ['../4_result/HMF/HNMF_result_lamta1&' num2str(lamta1)  '_lamta2&' num2str(lamta2)  '_t&' num2str(t)  '.mat'];
            if(~exist(fn,'file'))
                continue;
            end
            load(fn);
            location=(1-1)*b*c*d+(1-1)*c*d+(j-1)*d+k;
            [Z_filter{location,t},pathway_ncbi{location,t}]= predicted_pathway(W_out,T,mgi_id);
            [RD(location,t),F(location,t),Precision(location,t),Recall(location,t),jaccard(location,t)]=rand_index(Z_filter{location,t}*Z_filter{location,t}',mgi_id,standard,standard_row,alpha_ratio);
            rat(t) = ratio(pathway_ncbi{location,t});
        end
    end
end

rat = 1./rat;

end

