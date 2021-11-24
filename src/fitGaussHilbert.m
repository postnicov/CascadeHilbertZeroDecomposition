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

%% ��������� �� ���������
% ������������ ����� ������������� (������������ � ��� ������� �� �����  
% ����������� ���� ��������������� ���������) �� ������� ������ ��������� 
% ������������� ������� �������� ������� ����������
LinIntrv=20;
% ��������� ������� ������� ��������, � �������, � nsigma ���, �������
% ����� ������� ���������� ������� -- ���������� �������������� 
% ��������� � ����������� ����������� ����
nsigma=6;
% ��������� �������� ������������ ���������� �������������� � ������ 
val_corr=0.98;
% ��������� ��������, ��������� ���� �������� �� �����������
threshold=0.05;
% ������ �������� �� ��������� �����������������, ���� ��� ����
Npar=length(varargin)/2;
if Npar>0
    for k=1:Npar
        eval([genvarname(varargin{2*k-1}),'=varargin{2*k}']);
    end
end
% ������� ������, �������������� ���
Gauss=@(a,x) a(1)*exp(-(x-a(2)).^2/(2*a(3)^2));
% ��������� ������������� ���������
v=hilbert(y);
% ����� ����� �������������� ���������, ��������������� ���������� �������
L=length(y);
j=2:L;
Hy=imag(v);
% ����� - ����� �� ��������, ��� ������� ���������� ������� ������
% �������������� ��������� �� ������������� ������� � �������������;
ind0=find((Hy(j)>0)&(Hy(j-1)<0));
ind0bound=find((ind0>20)&(ind0<(L-20)));
ind0=ind0(ind0bound);
for k=1:length(ind0)
% �������� ������� ���������� � ������������ ����� ���������: ������ �� 
% j ����� �� ���������� ��������� �� 1 �� LinIntrv, ������� +/-j �� ���������,
% ���������� ��������� ���� �������������� ��������� �� ���� ��������� �
% �������� ������������ ���������� ����� ����� � �������� ��������
for j=1:LinIntrv;
    if ((ind0(k)-j)>1)&((ind0(k)+j)<L)
        xl=x([ind0(k)-j:ind0(k)+j]);
        yl=Hy([ind0(k)-j:ind0(k)+j]);
        pl=polyfit(xl,yl,1);
        cmatr=corrcoef(yl,polyval(pl,xl));
        cc(j)=cmatr(1,2);
    end
end
% ������� ����� ����������� - ��, �� ������� ����������� ����������
% ��������� val_corr %
    lcf=find(cc>val_corr);
    lc=length(lcf);
% ���������� ������������� �� ���� ������� � ������� ���������� ��������� 
% (lc>0 ����� ���� ���������� ����� ��� �������������, �.�. ��� ���������,
% �� 3 ����� (����������� ���� 1 � ��� �������) - ���������� ���������
% ����� ��� ��������
if lc>0
b{k}=nlinfit(x([ind0(k)-lc:ind0(k)+lc]),y([ind0(k)-lc:ind0(k)+lc]),Gauss,[y(ind0(k)),x(ind0(k)),x(ind0(k)+lc)-x(ind0(k))]);
else
    b{k}=[0,0,0];
end
if (b{k}(3)>(nsigma*(x(ind0(k)+lc)-x(ind0(k)))))
    b{k}=[0,0,0];
end
end
% ���������� �����
WaveNum=[];
n=1;
for j=1:length(b);
    nonEmpty=sum(b{j});
    if nonEmpty>0
       WaveNum(n,1:3)=b{j};
       n=n+1;
    end
end
% �������� ����, �� �������� �� ������� ������
[sw1,sw2]=size(WaveNum);
if (sw1*sw2)>0
% ���������� ������� ����� �������� �����, ������� tol
    idx=find(WaveNum(:,1)>threshold);
    b=WaveNum(idx,:);
else
    b=[];
end



