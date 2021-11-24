function b=fitGaussHilbert(x,y,varargin)
% One iteration of the Cascade Hilbert-Zero Decomposition
% 
% Input variables: independent variable x and the smaple's values y in the
% corresponding points
%
% Possible parameters:
% 'LinIntrv': The length (number of points) of the interval left and right 
% from the signal's Hankel transform zero crossing, within which 
% the approximation of the signal by the Gaussian is searched for.
% 'nsigma': The maximal dispersion of the identified Gaussian compoment, 
% which is considered as a realistic peak measured in units of the length
% of the interval, where the Hilbert tranform is linear in the vicinity of
% the signal's Hankel transform zero crossing
% 'val_corr': The minimal correlation between the signal's Hankel transform
% and a straint line, which allows deciding that the former is a linear 
% function within the interval chosen
% 'threshold': The threshold for idenitified Gaussian's amplitudes;
% magnitudes less than this value are not considered

%% Параметры по умолчанию
% Максимальная длина полуинтервала (отступаемого в обе стороны от точки  
% пересечения нуля преобразованием Гильберта) на котором ищется локальная 
% аппроксимация участка входного сигнала гауссианой
LinIntrv=20;
% Отсечение слишком широких гауссиан, с шириной, в nsigma раз, большей
% длины участка линейности функции -- результата преобразования 
% Гильберта в окрестности пересечения нуля
nsigma=6;
% Пороговое значение коэффициента корреляции преобразования с прямой 
val_corr=0.98;
% Пороговое значение, амплитуды ниже которого не учитываются
threshold=0.05;
% Замена значений по умолчанию пользовательскими, если они есть
Npar=length(varargin)/2;
if Npar>0
    for k=1:Npar
        eval([genvarname(varargin{2*k-1}),'=varargin{2*k}']);
    end
end
% Функция Гаусса, представляющая пик
Gauss=@(a,x) a(1)*exp(-(x-a(2)).^2/(2*a(3)^2));
% Численное пробразование Гильберта
v=hilbert(y);
% Поиск нулей преобразования Гильберта, соответствующих максимумам функции
L=length(y);
j=2:L;
Hy=imag(v);
% Метод - найти те элементы, для которых происходит переход кривой
% преобразования Гильберта из отрицательной области в положительную;
ind0=find((Hy(j)>0)&(Hy(j-1)<0));
ind0bound=find((ind0>20)&(ind0<(L-20)));
ind0=ind0(ind0bound);
for k=1:length(ind0)
% Проверка области линейности в окрестностях точки максимума: отступ на 
% j точек от найденного максимума от 1 до LinIntrv, выборка +/-j от максимума,
% вычисление линейного фита преобразования Гильберта на этом интервале и
% проверка коэффициента корреляции этого куска с линейной функцией
for j=1:LinIntrv;
    if ((ind0(k)-j)>1)&((ind0(k)+j)<L)
        xl=x([ind0(k)-j:ind0(k)+j]);
        yl=Hy([ind0(k)-j:ind0(k)+j]);
        pl=polyfit(xl,yl,1);
        cmatr=corrcoef(yl,polyval(pl,xl));
        cc(j)=cmatr(1,2);
    end
end
% Искомая длина окрестности - та, на которой коэффициент корреляции
% превышает val_corr %
    lcf=find(cc>val_corr);
    lc=length(lcf);
% Нелинейная аппроксимация на этом участке с поиском параметров гауссианы 
% (lc>0 чтобы было достаточно точек для аппроксимации, т.к. три параметра,
% то 3 точки (центральная плюс 1 в обе стороны) - минимально возможное
% число для фиттинга
if lc>0
b{k}=nlinfit(x([ind0(k)-lc:ind0(k)+lc]),y([ind0(k)-lc:ind0(k)+lc]),Gauss,[y(ind0(k)),x(ind0(k)),x(ind0(k)+lc)-x(ind0(k))]);
else
    b{k}=[0,0,0];
end
if (b{k}(3)>(nsigma*(x(ind0(k)+lc)-x(ind0(k)))))
    b{k}=[0,0,0];
end
end
% Исключение нулей
WaveNum=[];
n=1;
for j=1:length(b);
    nonEmpty=sum(b{j});
    if nonEmpty>0
       WaveNum(n,1:3)=b{j};
       n=n+1;
    end
end
% Проверка того, не является ли матрица пустой
[sw1,sw2]=size(WaveNum);
if (sw1*sw2)>0
% Исключение слишком малых амплитуд пиком, меньших tol
    idx=find(WaveNum(:,1)>threshold);
    b=WaveNum(idx,:);
else
    b=[];
end



