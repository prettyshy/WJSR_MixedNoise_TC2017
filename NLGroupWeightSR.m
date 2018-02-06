function [outpDimg,outpTemp]=NLGroupWeightSR(Nimg,Dimg,par,dict,flagNim,lambda1,lambda2)

n_im=double(Nimg);
d_im=double(Dimg);
stp=par.step;
D=dict;
p_size=par.patchsize;

 blkNimgMatrix                        =   Im2Patch( n_im,p_size,stp);
 [Arr  Wei]                           =    find_blks( d_im, par ); 
 blkNimgMatrixNL                        =   Im2Patch( n_im,p_size,1);

 % blkNimgMatrixDN                        =   Im2Patch( Dimg,p_size,1);
 % blkNimgMatrixDN2                       =   Im2Patch( Dimg,p_size,stp);
 
 blkflagMask = Im2Patch( flagNim, p_size, stp);
 blkflagMaskNL= Im2Patch( flagNim, p_size, 1);
 

blkNimgMatrix_mask=blkNimgMatrix.*blkflagMask;
patchNim_m=sum(blkNimgMatrix_mask)./(sum(blkflagMask));
patchNimMean=repmat(patchNim_m,size(blkNimgMatrix,1),1);
blkNimgMatrix=blkNimgMatrix-patchNimMean;


blkNimgMatrixNL_mask=blkNimgMatrixNL.*blkflagMaskNL;
patchNim_m=sum(blkNimgMatrixNL_mask)./(sum(blkflagMaskNL));
patchNimNLMean=repmat(patchNim_m,size(blkNimgMatrixNL,1),1);
blkNimgMatrixNL=blkNimgMatrixNL-patchNimNLMean;
 
 PNum=size(blkNimgMatrix,2);
 coefMatrix=zeros(size(D,2),PNum);
 
  h=waitbar(0,'W-SOMP on each group ...');
 for i=1:PNum
        waitbar(i/PNum);
     X=blkNimgMatrixNL(:,Arr(:,i));
     X=[blkNimgMatrix(:,i),X];   
     
     NLweit=blkflagMaskNL(:,Arr(:,i));
     NLweit=[blkflagMask(:,i),NLweit]; 
   
     [Coeff]=W_SOMP(D,X,NLweit,par.const*par.sig,10);
     coefMatrix(:,i)=Coeff(:,1); 
 end
   close(h);
 outpTemp=RecoverImage(Nimg,patchNimMean,D,coefMatrix,stp);  

 weit=(flagNim.*flagNim)+lambda1;
 G_y_x=(lambda1*(outpTemp-Nimg))./weit;
 U=lambda2*(1-flagNim)./weit;
 outpDimg=Nimg+shrink(G_y_x,U);
 
return;

function [im_out]=RecoverImage(y,p_mean,D,CoefMatrix,stp)
[h1,w1]=size(y); 
n=sqrt(size(D,1)); 
im_out=zeros(h1,w1); 
im_wei=zeros(h1,w1); 
Num=size(CoefMatrix,2);
patchMatrix=zeros(n^2,Num);
s=stp;
b=n;
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
 patchMatrix(:,k)=patch+p_mean(:,k);  
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