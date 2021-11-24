function B=fitGaussHilbertCascade(x,y,lvl,varargin)

% ‘ункци€ √аусса, представл€юща€ пик
Gauss=@(a,x) a(1)*exp(-(x-a(2)).^2/(2*a(3)^2));

b{1}=fitGaussHilbert(x,y,varargin);
B0{1}=b{1};
yfit{1}=x*0;
bsize=size(b{1});
for j=1:bsize(1);
    yfit{1}=yfit{1}+Gauss([b{1}(j,1),b{1}(j,2),b{1}(j,3)],x);
end
rez{1}=y-yfit{1};

if lvl>1
    for k=2:lvl
        b{k}=fitGaussHilbert(x,rez{k-1},varargin);
        bsize=size(b{k});
        yfit{k}=x*0;
        for j=1:bsize(1);
            yfit{k}=yfit{k}+Gauss([b{k}(j,1),b{k}(j,2),b{k}(j,3)],x);
        end
        bsizeprev=size(b{k-1});
        for j=1:bsizeprev(1);
            for m=1:bsize(1)
                b{k-1}(j,1)=b{k-1}(j,1)-Gauss([b{k}(m,1),b{k}(m,2),b{k}(m,3)],b{k-1}(j,2));
            end
        end
        yfit{k-1}=x*0;
        for j=1:bsizeprev(1);
            yfit{k-1}=yfit{k-1}+Gauss([b{k-1}(j,1),b{k-1}(j,2),b{k-1}(j,3)],x);
        end
        rez{k}=y;
        for m=1:k
            rez{k}=rez{k}-yfit{m};
        end
    end
end

bsize=size(b);
n=1;
B=[];
for j=1:bsize(2);
    B=[B;b{j}];
end
[~,idx] = sort(B(:,2));
B=B(idx,:);
