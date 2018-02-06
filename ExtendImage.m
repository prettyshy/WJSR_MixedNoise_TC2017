function outpt=ExtendImage(image,n)

f=image; 
[row,col]=size(f); 
K=row*col;          
mid=zeros(row+2*n,col+2*n);  

mid((n+1):(row+n),(n+1):(col+n))=f;

for i=1:n
    mid(i,1:n)=f(i,1:n); 
     mid(i,(n+1):(col+n))=f(i+1,:); 
    mid(i,(col+n+1):(col+2*n))=f(i,(col-n+1):col);
end   
for i=(row+n+1):(row+2*n) 
    mid(i,1:n)=f(i-2*n,1:n);
     mid(i,(n+1):(col+n))=f(i-2*n-1,:); 
    mid(i,(col+n+1):(col+2*n))=f(i-2*n,(col-n+1):col);
end

for j=1:n
     mid((n+1):(row+n),j)=f(:,j+1); 
end   
for j=(col+n+1):(col+2*n)
    mid((n+1):(row+n),j)=f(:,j-2*n-1); 
end
outpt=mid;