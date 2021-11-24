function [yfit,Yfit]=fitGaussHilbertCurves(x,b)

% ‘ункци€ √аусса, представл€юща€ пик
Gauss=@(a,x) a(1)*exp(-(x-a(2)).^2/(2*a(3)^2));

bsize=size(b);
for j=1:bsize(1)
    yfit(:,j) = Gauss(b(j,:),x);
end
Yfit=sum(yfit')';