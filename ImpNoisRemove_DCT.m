function [Result1,Dict]=ImpNoisRemove_DCT(IMin0,Nimg,flag)

bb=8; 
stp=3;
K=256; 


IMin0=double(IMin0);


DCT=zeros(bb,sqrt(K));
for k=0:1:sqrt(K)-1,
    V=cos([0:1:bb-1]'*k*pi/sqrt(K));
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);



IMin=Nimg;


blkMatrixIm=Im2Patch(IMin,bb,stp);

flag1=flag;
blkMask=Im2Patch(flag1,bb,stp);

[Coeff]=OMPerrWeight(DCT,blkMatrixIm,blkMask,5); 


 Result1=RecoverImage(IMin,DCT,Coeff,stp); 

Result1=max(min(Result1,255),0);

Dict=DCT;
return;


function [A]=OMPerrWeight(D,X,Mask,errorGoal)

[n,P]=size(X);
[n,K]=size(D);
E2 = errorGoal^2*n;
maxNumCoef = 10;
 A = sparse(size(D,2),size(X,2));
for k=1:1:P,
    a=[];
    Mpos=find(Mask(:,k)); 
    Msup=Mask(Mpos,k);
    E2M=E2*length(Mpos)/n;
    x=X(Mpos,k);
    MaskMatr=kron(ones(1,K),Msup);
     Dict=D(Mpos,:);
     Dict=MaskMatr.*Dict;
     W=1./(sqrt(diag(Dict'*Dict))+eps); 
     Dict=Dict*diag(W);
    residual=x;
    indx = [];
    a = [];
    currResNorm2 = sum(residual.^2);
    j = 0;
    while currResNorm2>E2M && j < maxNumCoef,
        j = j+1;
        proj=Dict'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(Dict(:,indx(1:j)))*x;
        residual=x-Dict(:,indx(1:j))*a;
        currResNorm2=sum(residual.^2);
    end;
    if (~isempty(indx))
        A(indx,k)=a;
        A(:,k)=W.*A(:,k);
    end
end;
return;
% ========================================================
% ========================================================

%% averaging the denoised patches as an image 
function [im_out]=RecoverImage(y,D,CoefMatrix,stp)
[h1,w1]=size(y); 
n=sqrt(size(D,1)); 
im_out=zeros(h1,w1); 
im_wei=zeros(h1,w1); 
Num=size(CoefMatrix,2);
patchMatrix=zeros(n^2,Num);
s=stp;
b=8;
N        =  h1-b+1;
M        =  w1-b+1;
r        =  [1:s:N];
r        =  [r r(end)+1:N];
c        =  [1:s:M];         
c        =  [c c(end)+1:M];
N        =  length(r);   
M        =  length(c);   
for k=1:1:Num,
patch=D*CoefMatrix(:,k);
patchMatrix(:,k)=patch;
end;
kk=0;
for i  = 1:n
    for j  = 1:n
        kk    =  kk+1;
        im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( patchMatrix(kk,:)', [N M]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;       
    end
end
im_out  =  im_out./im_wei;
return;
