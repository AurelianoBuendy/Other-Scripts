function [beta, I_c] = SI_2019ncov_ga
    I=[291,440,571,830,1287,1975,2744,4515,5974,7711,9692,11791,14380,...
        17205,20438,24324,28018,31161,34546,37198,40171];
    t=1:21;
    N=1400000000;
    I0=I(1,1);
    S0=N-I0;
    a=S0/I0;
    options = gaoptimset('PopulationSize',200,'PopulationType',...
        'doubleVector','PlotFcns',{@gaplotbestf,@gaplotbestindiv,...
        @gaplotexpectation,@gaplotstopping});
    tmp = ga(@ga_fit, 1, [],[],[],[],[],[],[],options);
    %N = tmp(1);
    %a = tmp(2);
    beta = tmp(1);
    %t2 = t + tmp(2);
    t2 = 1:37;
    I_c = N./(a.*exp(-beta.*t2) + 1);
    %date = {'Jan_20','Jan_21','Jan_22','Jan_23','Jan_24','Jan_25',...
    %     'Jan_26','Jan_27','Jan_28','Jan_29','Jan_30','Jan_31','Feb_1',...
    %     'Feb_2','Feb_3','Feb_4','Feb_5','Feb_6','Feb_7','Feb_8',...
    %     'Feb_9','Feb_10','Feb_11','Feb_12','Feb_13','Feb_14','Feb_15',...
    %     'Feb_16','Feb_17','Feb_18','Feb_19','Feb_20','Feb_21'};
    date = [2020,1,20;2020,1,21;2020,1,22;2020,1,23;2020,1,24;...
        2020,1,25;2020,1,26;2020,1,27;2020,1,28;2020,1,29;2020,1,30;...
        2020,1,31;2020,2,1;2020,2,2;2020,2,3;2020,2,4;2020,2,5;2020,2,6;...
        2020,2,7;2020,2,8;2020,2,9;2020,2,10;2020,2,11;2020,2,12;...
        2020,2,13;2020,2,14;2020,2,15;2020,2,16;2020,2,17;2020,2,18;...
        2020,2,19;2020,2,20;2020,2,21;2020,2,22;2020,2,23;2020,2,24;...
        2020,2,25];
    t_date = datenum(date(1:21,1),date(1:21,2),date(1:21,3),0,0,0);
    t_date2 = datenum(date(1:37,1),date(1:37,2),date(1:37,3),0,0,0);
    figure;
    plot(t_date2+1, I_c);
    hold on
    plot(t_date+1, I, '*');
    datetick('x', 6);
    %set(gca,'XLim',[t_date2(1)+1, t_date2(end)+1]);%X轴的数据显示范围
    %set(gca,'XTick',t_date2(1)+1:5:t_date2(end)+1);%设置要显示坐标刻度
    set(gca,'FontSize',8.0,'Fontname', 'Times New Roman');
    %set(gca,'XTickLabel',date{1,1:21});%给坐标加标签 
    %set(gca,'YLim',[0, 100000]);%Y轴的数据显示范围
    
function y = ga_fit(x)
    I = [291,440,571,830,1287,1975,2744,4515,5974,7711,9692,11791,14380,...
    17205,20438,24324,28018,31161,34546,37198,40171];
    t = 1:21;
    N = 1400000000;
    I0 = I(1,1);
    S0 = N - I0;
    a = S0/I0;
    I_c = N./(a.*exp(-x(1).*t) + 1);
    %I_c = x(1)./(x(2).*exp(-x(3).*t) + 1);
    sum = 0;
    for loopi = 1:21
        sum = sum + (I_c(1,loopi)-I(1,loopi))^2;
    end
    y = sum;