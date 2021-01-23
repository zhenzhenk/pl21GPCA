function [U,V,W1]=pl21PCA(X,W,lamda,alpha,k,miu,rho,p)
%% %% pL21PCA Based on Graph Laplacian Regularization and P-Norm for Gene Selection and Clustering
%% Written by Xiangzhen Kong in Qufu Normal University 2019/11/23
% input: vector data X (size:m*n,m:orignal data dim,n:data num); 
%        graph data W(size:n*n)
%        lamda:weighting parameter for U
%        alpha:weighting parameter for 
%       
%        k:reduced dimension
%        Y:Lagrangian multipliers
% output:embedding vector V (size:n*k,k:reduction dim);
%        projection U (size:m*k)
%        auxiliary variable E (size: m*n)
% obj: J = ||X-UV^T||^p+lamda*||U||2,1+alf*trV^T*(D-W)*V    s.t. V^T*V=I


[m,n]=size(X);
epsilon=10^-3;
Y=zeros(m,n);
W1=zeros(m,n);
d=sum(W);
D=diag(d);
L=D-W;
[U,S,V]=svd(X,'econ');
% U=U(:,1:k)*S(1:k,1:k);
% V=V(:,1:k);
U=rand(m,k);
V=rand(n,k);
maxmiu=10^10;
% while norm(W1)<epsilon
flag=1;
count=0; 
new_folder=strcat('E:/simulated data/PL21GPCA/Resultp',num2str(p));
 mkdir(new_folder);
while flag==1
 %% 
    DD=zeros(m,m);
    for i=1:m
         DD(i,i)=norm(U(i,:))^-1;
    end
    I1=eye(m,m);
    A=(I1+(2*(lamda/miu)*DD))^-1;
    H=X-W1-Y/miu;
    K=-H'*A*H+(alpha/miu)*L;
    U=A*H*V;    
 %%
    [U0,V0] = eig(K);
    val = diag(V0);
    [~,ind] = sort(val);
    V = U0(:,ind(1:k));
   % U=A*H*V;
 %%
    T=X-U*V'-Y/miu;
    a=1/miu;
    pe=a*(abs(T).^(p-1));
    e=abs(T)-pe;
    W2=zeros(size(T));
    ID2=W2<e;
    W2(ID2)=e(ID2);
    W1=(W2.*T)./abs(T);
   
    
%     DD2=(e.*T)./abs(T);
% %     ID=abs(T)>pe;
% %     W1(ID)=DD2(ID);
%     ID=W1<DD2;
%     W1(ID)=DD2(ID);
 normw1=norm(W1);
 count=count+1;
 filepath=pwd;
 cd(new_folder);
 filename=strcat('GLUVW',num2str(count),'.mat');
 save(filename,'U','V','W1','count','normw1');
 cd(filepath);
    Y=Y+miu*(W1-X+U*V');
    miu=rho*miu;
%      if count==5 || norm(W1)<epsilon
         if count==5 
%        if norm(W1)<epsilon
     flag=0;
     end
end
% end
    
    