function [A]=W_SOMP(D,X,Mask,errorGoal,maxNumCoef)
% ========================================================
% W_SOMP: Weighted Simultaneous Orthogonal Matching Pursuit,  Weighted
% Simultaneous Sparse coding of a group of signals based on a given
% dictionary;
% error to get.

% input arguments: D - the dictionary
%                           X - the signals matrix to represent
%                           errorGoal - the maximal allowed representation error
% output arguments: A - sparse coefficient matrix which is expected to be row sparse.

% the model: || W ¡C\odot(X-DA)||_F^2, s.t. ||A||_2,1<=maxNumCoef. 
% where \odot is the Hadamard product.

% written by Licheng Liu, 2014.07.30.
% ========================================================

[n,P]=size(X);
[n,K]=size(D);
E2 = errorGoal^2*n;

A = sparse(size(D,2),size(X,2));  

    X=Mask.*X; 
    R=X;
    indx = [];
    A_s = [];
    currResNorm2 = sum(sum(R.^2))/P;
    j = 0;
    
    Dict=zeros(n,K,P);
    W=zeros(K,P);
    for i=1:P
        Msup=Mask(:,i);
        MaskMatr=kron(ones(1,K),Msup);
        dd=MaskMatr.*D;  
        
        W(:,i)=1./(sqrt(diag(dd'*dd))+eps); 
        Dict(:,:,i)=dd*diag(W(:,i)); 
    end
    
while currResNorm2>E2 && j < maxNumCoef,
  for k=1:1:P,
      dic=Dict(:,:,k);
      proj(:,k)=(dic')*R(:,k); 
  end
    j=j+1;
    projEachD=sum(abs(proj),2);
     pos=find(abs(projEachD)==max(abs(projEachD)));
     pos=pos(1);
      indx(j)=pos;
     
      A_s=zeros(j,P);
      for k=1:1:P
          x_s=X(:,k);
          Dct=Dict(:,:,k);
          D_s=Dct(:,indx(1:j));
          alpha=pinv(D_s)*x_s;
          A_s(:,k)=alpha;
          r=x_s-D_s*alpha;
          R(:,k)=r;
      end
      currResNorm2 = sum(sum(R.^2))/P;
end
        if (~isempty(indx))
          A(indx,:)=A_s;
          A(:,1)=W(:,1).*A(:,1);
        end
return;

