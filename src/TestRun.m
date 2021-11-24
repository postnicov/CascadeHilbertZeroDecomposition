% Processing spectra with CHZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;
% Printing options
printfigures=1;% 1 - "print"; 0 - "do not print";
frmt='-dpng';
res=200;%resolution in dpi
FS=19;%Font size for figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Загрузка данных для обработки
% % Вариант, когда таблица сохранена в mat-file
% load dataset109
% Вариант с с чтением из текстового файла
%Имя текстового файла:
fname='sample109';%'Внелег_109_633_80sec_10re_отрезано_от_600_до_1800нормировано_коррекция_баз_лин';
[x,y]=textread([fname,'.txt'],'%f %f');
% Перенормировка y с максимума 100 на максимум 1
y=y/100;
% Если х в исходной таблице идет по убыванию, то проводится пересортировка  
% по возрастанию
if x(2)<x(1)
   x=flipud(x);
   y=flipud(y);
end
%% Каскадное преобразование
% Количество уровней (шагов) преобразования
lvls=3;
% Выполение преобразования; rолонки выходной таблицы B: амплитуда пика, 
% его положение (упорядочено по его возрастанию), дисперсия (ширина) пика
B=fitGaussHilbertCascade(x,y,lvls);
% Отдельные пики и их сумма
[yfit,Y]=fitGaussHilbertCurves(x,B);
%% График
figure(1)
% Зеленая линия - исходных спектр
plot(x,y,'-','color','green','LineWidth',1.5)
hold on
% Тонкие цветные линии - индивидуальные пики
plot(x,yfit)
% Черная штриховая - его полная аппрокисимация
plot(x,Y,'--','color','black','LineWidth',1)
xlim([x(1) x(end)])
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman intensity')
set(gca,'FontSize',FS)

if printfigures==1
  % Сохранение графика  
  figure(1)
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
  print('Figure_Gaussians',frmt,['-r',num2str(res)])
  % Сохранение таблицы с параметрами пиков
  xlswrite(['param_Gaussians.xlsx'],B)
end
