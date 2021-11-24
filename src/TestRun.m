% Processing spectra with CHZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
% Printing options
printfigures=1;% 1 - "print"; 0 - "do not print";
frmt='-dpng';
res=200;%resolution in dpi
FS=19;%Font size for figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������� ������ ��� ���������
% % �������, ����� ������� ��������� � mat-file
% load dataset109
% ������� � � ������� �� ���������� �����
%��� ���������� �����:
fname='sample109';%'������_109_633_80sec_10re_��������_��_600_��_1800�����������_���������_���_���';
[x,y]=textread([fname,'.txt'],'%f %f');
% �������������� y � ��������� 100 �� �������� 1
y=y/100;
% ���� � � �������� ������� ���� �� ��������, �� ���������� ��������������  
% �� �����������
if x(2)<x(1)
   x=flipud(x);
   y=flipud(y);
end
%% ��������� ��������������
% ���������� ������� (�����) ��������������
lvls=3;
% ��������� ��������������; r������ �������� ������� B: ��������� ����, 
% ��� ��������� (����������� �� ��� �����������), ��������� (������) ����
B=fitGaussHilbertCascade(x,y,lvls);
% ��������� ���� � �� �����
[yfit,Y]=fitGaussHilbertCurves(x,B);
%% ������
figure(1)
% ������� ����� - �������� ������
plot(x,y,'-','color','green','LineWidth',1.5)
hold on
% ������ ������� ����� - �������������� ����
plot(x,yfit)
% ������ ��������� - ��� ������ ��������������
plot(x,Y,'--','color','black','LineWidth',1)
xlim([x(1) x(end)])
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman intensity')
set(gca,'FontSize',FS)

if printfigures==1
  % ���������� �������  
  figure(1)
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
  print('Figure_Gaussians',frmt,['-r',num2str(res)])
  % ���������� ������� � ����������� �����
  xlswrite(['param_Gaussians.xlsx'],B)
end
