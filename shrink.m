function yshk=shrink(b,a)

b=double(b);
a=double(a);
yshk=(abs(b)>=a).*(abs(b)-a).*sign(b);
return