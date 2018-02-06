function [PSNR, MSE] = psnr(x0, y0)
% ???????PSNR??????MSE
 
% By lyqmath
 
% Dalian University of Technology

% School of Mathematical Sciences
 
% ???????PSNR??????MSE
 
% ????Y??????X???????PSNR?MSE
 
X=double(x0);
Y=double(y0);
if nargin<2
 
    D = X;
 
else
 
    if any(size(X)~=size(Y))
 
        error('The input size is not equal to each other!');
 
    end
 
    D = X-Y;
 
end
 
MSE = sum(D(:).*D(:))/prod(size(X));
 
PSNR = 10*log10(255^2/MSE);