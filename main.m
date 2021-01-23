clc;
clear;
close all;
data=generate_SNR_data(10);
%% 对数据进行列去均值处理

% [m n]=size(data);
% mn=mean(data,1);
% data=data-repmat(mn,m,1);
%% 
%% 对数据进行行去均值
% [m n]=size(data);
% mn=mean(data,2);
% data=data-repmat(mn,1,n);
%% 
fea=data';
options=[];
options.k = 5;
% options.WeightMode = 'HeatKernel';
options.We2145ightMode = 'Binary';
% options.WeightMode = 'Cosine';
options.t = 1;
W = constructW(fea,options);
lamda=10;
alpha =100;
k=4;
miu = 10^-2;
rho = 1.2; % 1<rho<1.5, smaller rho gives better solution but slower speed
X=fea';
p=0.1;
for i=1:9
    [U,V,W1] = pl21PCA(X,W,lamda,alpha,k,miu,rho,p);
    p=p+0.1;
end
toc

