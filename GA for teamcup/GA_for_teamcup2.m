function [Xp, LC1, LC2, Xp_f, LC1_f, LC2_f] = GA_for_teamcup2
    clf;
    popsize = 200; %poplutation size
    gene = 200; %hereditary algebra
    pm = 0.8; %variation frequence
    err = 0.001;%error
    [Xp, LC1, LC2] = teamcup_GA(popsize, gene, pm, err, 1);
    [Xp_f, LC1_f, LC2_f] = teamcup_GA(popsize, gene, pm, err, -1);
    
    

function [Xp, LC1, LC2] = teamcup_GA(popsize, gene, pm, err, flag_num)
   %% 输入参数列表
    %  gene     遗传进化迭代次数
    %  popsize     种群规模
    %  pm    变异概率
   %% 输出参数列表
    %  Xp    最优个体
    %  LC1   目标的收敛曲线
    %  LC2   平均适应度函数的收敛曲线
   %% plot convergence curve
    if flag_num == 1
        subplot(2,2,2);
        xlabel('迭代次数');
        ylabel('最大成功概率的最优值');
        title('带复活甲最大成功概率函数的收敛曲线');
        drawnow
        hold on
        subplot(2,2,1);
        xlabel('迭代次数');
        ylabel('适应度的平均值');
        title('带复活甲最大成功概率适应度的收敛曲线');
        drawnow
        hold on
    else
        subplot(2,2,4);
        xlabel('迭代次数');
        ylabel('最低成功概率的最优值');
        title('带复活甲最低成功概率函数的收敛曲线');
        drawnow
        hold on
        subplot(2,2,3);
        xlabel('迭代次数');
        ylabel('适应度的平均值');
        title('带复活甲最低成功概率适应度的收敛曲线');
        drawnow
        hold on
    end
    
    %% initial population
    Xp = randperm_num(98, 5);
    LC1 = zeros(1,gene);
    LC2 = zeros(1,gene);
    Best = inf;
    farmer = cell(1, popsize);
    k = 1;
    while k <= popsize
        tmp = randperm_num(98, 5);
        x = sort(tmp + 1);
        farmer{1, k} = x;
        k = k + 1;
    end
    %% evolution
    g = 1;
    counter = 0;
    while g <= gene
       %% cross
        k = 1;
        newfarmer = cell(1, 2 * popsize);
        while k <= popsize
            flag=randperm_num( popsize, 2);
            father = farmer{1, flag(1)}; %father
            mother = farmer{1, flag(2)}; %mother
            tmp = [father, mother];
            tmp = sort(tmp);
            tmp = unique(tmp);
            child1 = tmp(1: 5); 
            child2 = tmp(end - 4: end);
            newfarmer{1, 2 * k - 1} = child1;
            newfarmer{1, 2 * k} = child2;
            k = k + 1;
        end
        allfarmer = [farmer, newfarmer];
       %% variation
        k = 1;
        while k <= 3 * popsize
            if rand(1) <= pm
                x = allfarmer{1, k};
                label_flag = 1;
                while label_flag > 0
                    tmp = randperm_num(98, 1);
                    replace_num = tmp + 1;
                    label_flag = length(find(x == replace_num));
                end
                label = randperm_num(5, 1);
                x(label) = replace_num;
                allfarmer{1, k} = sort(x);
            end
            k = k + 1;
        end
       %% selection
        k = 1;
        freq = ones(1, 3 * popsize);
        opti_freq = zeros(1,popsize);
        while k <= 3 * popsize
            x = allfarmer{1, k};
            if flag_num == 1
                freq(1, k) = fitness_func(x);
            else
                freq(1, k) = 1 - fitness_func(x) ;
            end
            %frequece of fail or victory   
            k = k + 1;
        end     
        k = 1;
        while k <= popsize
            freqmin = min(freq);
            posfreq = find(freq == freqmin);
            pos = posfreq(1);
            farmer{k} = allfarmer{pos};
            opti_freq(k) = freq(pos);
            k = k + 1;
            freq(pos) = inf; 
        end %select the best popsize
        %choose the best of farmers
        minfreq = min(opti_freq);
        meanfreq = mean(opti_freq);
        pos = find(opti_freq == minfreq);
        LC2(g) = meanfreq;
        if minfreq < Best 
            Best = minfreq;
            Xp = farmer{pos(1)};
        end
        if flag_num == 1
            LC1(g) = 1 - Best;
            subplot(2,2,2);
            plot(g, LC1(g), '*');
            legend(['最大成功概率:' num2str(LC1(g))]);
            drawnow
            hold on
            subplot(2,2,1);
            plot(g,LC2(g),'*');
            drawnow
            hold on
        else
            LC1(g) = Best;
            subplot(2,2,4);
            plot(g, LC1(g), '*'); 
            legend(['最低成功概率:' num2str(LC1(g))]);
            drawnow
            hold on
            subplot(2,2,3);
            plot(g,LC2(g),'*');
            drawnow
            hold on
        end
        if g >= 2 && err > abs(LC1(g) - LC1(g - 1))
            counter = counter + 1;
        end
        if counter <= 20
            tmp = randperm_num(98, 5);
            x = sort(tmp + 1);
            farmer{1, randperm_num(popsize, 1)} = x;
        else
            break
        end
        g = g + 1;
    end

  
function frequence = fitness_func_markov(x)
    % fitness function
    fail = 0;
    for loopi = 1 : 1000
        steps = 1;
        while steps < 100
            dice = randperm_num(6, 1);
            steps = steps + dice;
            tmp = length(find(x == steps));
            if tmp > 0
                fail = fail + 1;
                break
            end
        end
    end
    frequence = fail/1000;

function failfreq = fitness_func_no_defense(x)
    % fitness function
    pos = x - 1;
    freq = zeros(1, 104);
    freq(6) = 1;
    if length(pos) == 1
        for loopi = 7:104
            freq(loopi) = 1/6.0 * (freq(loopi - 1) + freq(loopi - 2) +...
            freq(loopi - 3) + freq(loopi - 4) + freq(loopi - 5) +...
            freq(loopi - 6) );
            if loopi - 6 == pos
                failfreq = freq(loopi);
                break
            end
        end
    else
        k = 1;
        for loopi = 7:104
            freq(loopi) = 1/6.0 * (freq(loopi - 1) + freq(loopi - 2) +...
            freq(loopi - 3) + freq(loopi - 4) + freq(loopi - 5) +...
            freq(loopi - 6) );
            if loopi - 6 == pos(end)
                failfreq = freq(loopi) + fitness_func_no_defense(x(1: end - 1));
                break
            end
            if loopi - 6 == pos(k)
                freq(loopi) = 0;
                k = k + 1;
            end
        end
    end
    
    function failfreq = fitness_func(x)
    pos = sort(x - 1);
    failfreq = 0;
    for loopi = 1: 4
        freq = zeros(1, 104);
        freq(pos(loopi) + 6) = get_the_pos(pos(1:loopi) + 1);
        k = loopi + 1;
        for loopj = 7:104
            if loopj - 6 > pos(loopi)
                freq(loopj) = 1/6.0 * (freq(loopj - 1) + freq(loopj - 2) +...
                freq(loopj - 3) + freq(loopj - 4) + freq(loopj - 5) +...
                freq(loopj - 6) );
            end
            if loopj - 6 == pos(k)
                failfreq = failfreq + freq(loopj);
                freq(loopj) = 0;
                k = k + 1;                
            end
            if k > 5
                break
            end
        end
    end
    
    function failfreq = get_the_pos(x)
    % fitness function
    pos = x - 1;
    freq = zeros(1, 104);
    freq(6) = 1;
    if length(pos) == 1
        for loopi = 7:104
            freq(loopi) = 1/6.0 * (freq(loopi - 1) + freq(loopi - 2) +...
            freq(loopi - 3) + freq(loopi - 4) + freq(loopi - 5) +...
            freq(loopi - 6) );
            if loopi - 6 == pos
                failfreq = freq(loopi);
                break
            end
        end
    else
        k = 1;
        for loopi = 7:104
            freq(loopi) = 1/6.0 * (freq(loopi - 1) + freq(loopi - 2) +...
            freq(loopi - 3) + freq(loopi - 4) + freq(loopi - 5) +...
            freq(loopi - 6) );
            if loopi - 6 == pos(end)
                failfreq = freq(loopi);
                break
            end
            if loopi - 6 == pos(k)
                freq(loopi) = 0;
                k = k + 1;
            end
        end
    end
    
function rand_list = randperm_num(num, nums)
    tmp = randperm(num);
    rand_list = tmp(1: nums);