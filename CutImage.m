function outpt=CutImage(image,n)

f=image; 
[row,col]=size(f);  
K=row*col;        
mid=zeros(row-2*n,col-2*n); 

mid=f((n+1):(row-n),(n+1):(col-n));

outpt=mid;