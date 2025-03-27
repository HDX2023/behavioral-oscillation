clear
close all
t = (-100:10:100) / 1000;
load('C:\Users\hdxst\Desktop\behave oscillation\double flash behavioral oscillation\徐哲浩-毕业论文数据-不同的声音呈现时间对视觉整合的影响研究\test\selected data ana\subdata_selected.mat', 'subdata_selected');
%load('C:\Users\hdxst\Desktop\behave oscillation\double flash behavioral oscillation\徐哲浩-毕业论文数据-不同的声音呈现时间对视觉整合的影响研究\test\subdata_new.mat','subdata_checked');
data=subdata_selected.*100;
%% Detrend of raw data
[num_subjects, num_timepoints] = size(data);

% initialize
detrended_data = zeros(num_subjects, num_timepoints);

for i = 1:num_subjects
    % current subject's data
    y = data(i, :);
    
    % time series
    %x = (1:num_timepoints)';
    x = (1:num_timepoints);
    
    % 一次多项式拟合?
    p = polyfit(x, y, 2); %多项式系数polynomial coefficient
    
    % 多项式拟合结果计算数据的趋势
    trend(i,:) = polyval(p, x); %trend_value for each time(x)
    
    % 去除趋势成分
    detrended_data(i, :) = y - trend(i,:);
end
%% plot the trend
figure;
plot(t*1000,mean(detrended_data,1),'-o');
hold on;
title('Detrended Proportion time course')
xlim([-110,110]);
%ylim([-0.08,0.08]);
%yticks(-0.05:0.05:0.05)
xlabel('Audio-visual lags (ms)')
ylabel('Proportion')
set(gca,'FontSize', 14);
% ylim([-10 10]);
% yticks(-10:5:10);

figure;
plot(t*1000,mean(trend,1),'-o');
ylim([20 60]);
yticks(20:10:60);

%% fft power spctrm for each sub
% initialize parameters
var_soa = 10; % in ms
samplerate = 1/var_soa * 1000; % in hz

t_padded = (-500:10:500) / 1000; % in second
N = length(t_padded); % data points 201
detrended_data_padded = [zeros(num_subjects,(N-length(t))/2),detrended_data,zeros(num_subjects,(N-length(t))/2)];% padding data points = length(t_padded)
%detrended_data_padded = [zeros(num_subjects,(N-length(t))/2),data,zeros(num_subjects,(N-length(t))/2)];% padding data points = length(t_padded)

T = [1:N]/samplerate; % checked data length 201

% spectral analysis
data_freq = fft(detrended_data_padded,[],2);
PS = abs(data_freq/N).^2;
singleside_PS = PS(:,1:N/2+1);
singleside_PS(:,2:end-1) = 2*singleside_PS(:,2:end-1);
faxis = samplerate/2*linspace(0,1,N/2+1);% [0*samplerate/2,1*samplerate/2]中间打了N/2+1个频率点

figure;
plot(faxis,mean(singleside_PS,1));
title('Average Group Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([5,30]);
%% permutation test
% remember to clarify the permutation vars !
PS_interest = singleside_PS;

interest_freq_index = (faxis >= 5 & faxis <= 30);
obs_statistics = PS_interest(:,interest_freq_index);
freq = faxis(:,interest_freq_index);
num_permute = 1000;
permute_stats = zeros(num_permute,size(obs_statistics,1),size(obs_statistics,2));

for ti = 1:num_permute
    
    %permute_detrend_data = detrended_data(:,randperm(size(detrended_data,2)));
    permute_data = data(:,randperm(size(data,2)));
    var_soa = 10; % in ms
    samplerate = 1/var_soa * 1000; % 100 in hz
    permute_tpad = (-500:10:500)/1000; % in second
    
    %set the detrend method
    %split-detrend_data(permute_data, t, 'split_poly2', struct()); 
    %ushape-detrend_data(permute_data, t, 'ushape_poly2', struct());
    %sg- params.w =;params.polyOrder = 2;detrended_data = detrend_data(permute_data, t, 'sgolay', params);
    %movmean - params.window_size =;detrended_data = detrend_data(permute_data, t, 'movmean', params);
    params.w =9; params.polyOrder =2;
    permute_detrend_data = detrend_data(permute_data, t, 'ushape_poly2', struct());
    
    permute_N = length(permute_tpad); % 101 data points
    permute_datapad = [zeros(num_subjects,(permute_N-length(t))/2),permute_detrend_data,zeros(num_subjects,(permute_N-length(t))/2)];% t=21;length pad 101
    
    permutedata_freq = fft(permute_datapad,[],2);
    permute_PS = abs(permutedata_freq/permute_N).^2;
    permuteside_PS = permute_PS(:,1:permute_N/2+1);
    permuteside_PS(:,2:end-1) = 2*permuteside_PS(:,2:end-1);
    permute_faxis = samplerate/2*linspace(0,1,permute_N/2+1);% [0*samplerate/2,1*samplerate/2]中间打了N/2+1个频率点
    
    permute_stat = permuteside_PS(:,interest_freq_index);
    permute_stats(ti,:,:) = permute_stat;
end

p_value = sum(squeeze(mean(permute_stats,2))>= mean(obs_statistics,1),1) / num_permute;
% plot pvalue
figure;
plot(freq, p_value, '-o');
xlabel('Frequency (Hz)');
ylabel('P-value');
title('Permutation Test P-values');
ylim([0 1]);
yline(0.05, '--r', 'alpha=0.05');
%% correction
% label 95th data
alpha = 0.05;
% 求每次置换的频率点对应的power最大值
perm_stat_avg = squeeze(sort(mean(permute_stats,2),1,'ascend'));
permute_stat_95th = perm_stat_avg(round(0.95 *num_permute),:); % 1000 * 1

% max 95th value
% sorted_max_permute = sort(max_permute_stat,'ascend');
% max_95th = sorted_max_permute(round((1-alpha)*num_permute));

%% plot
% setting defaults
set(0, 'DefaultTextColor', 'black');        % 璁剧疆榛樿鐨勬枃鏈鑹?
set(0, 'DefaultAxesXColor', 'black');       % 璁剧疆 x 杞寸殑榛樿棰滆壊
set(0, 'DefaultAxesYColor', 'black');       % 璁剧疆 y 杞寸殑榛樿棰滆壊
myfontsize = 14;
mypapersize = [4,3];
mypaperpos = [0,0,4,3];
myfontname = 'Arial';
hexStr = {'#bababa', '#e6550d','#2171b5'};
color = hex2rgb(hexStr);

% plot data (soa function)
soa = -100:10:100;
figure;
s2=shadedErrorBar(soa,mean(data,1)/100,(std(data,0,1)/100)./sqrt(size(data,1)),'lineProps',{'-','color',color(2,:)});
s2.mainLine.LineWidth = 1.5;
title('Proportion time course')
xlim([-110,110]);
ylim([0.25,0.55]);
yticks(0.30:0.10:0.50);
xlabel('Audio-visual lags(ms)')
ylabel('Proportion of 2-flash percepts')
set(gca,'FontSize', 15);
h=gcf;
box off
h.PaperUnits='inch';
h.PaperSize = mypapersize;
h.PaperPosition=mypaperpos;

%% plot trend
soa = -100:10:100;
figure;
s7=shadedErrorBar(soa,mean(trend,1)/100,(std(trend,0,1)/100)./sqrt(size(trend,1)),'lineProps',{'-','color',color(2,:)});
s7.mainLine.LineWidth = 1.5;
title('trend data time course')
xlabel('Audio-visual lags(ms)')
ylabel('Proportion of 2-flash percepts')
set(gca,'FontSize', 15);
h=gcf;
box off
h.PaperUnits='inch';
h.PaperSize = mypapersize;
h.PaperPosition=mypaperpos;

%% plot detrend data
figure;
s3 = shadedErrorBar(soa,mean(detrended_data,1),std(detrended_data,0,1) ./ sqrt(size(detrended_data,1)),'lineProps',{'-','color',color(2,:)});
s3.mainLine.LineWidth = 1.5;
title('Detrended Proportion time course')
xlim([-110,110]);
%ylim([-0.08,0.08]);
%yticks(-0.05:0.05:0.05)
xlabel('Audio-visual lags (ms)')
ylabel('Proportion')
set(gca,'FontSize', 14);
h=gcf;
box off
h.PaperUnits='inch';
h.PaperSize = mypapersize;
h.PaperPosition=mypaperpos;

%% fft
figure;
%corrected_pavlues = repmat(max(null_distribution_95th),1,size(null_distribution_95th,2));
%plot(faxis(interest_freq_index),corrected_pavlues,'--k','LineWidth', 0.5);
plot(faxis(interest_freq_index), ...
     repmat(max(permute_stat_95th), 1, sum(interest_freq_index)), ...
     '--k', 'LineWidth', 0.5);
hold on 
mean_power = mean(obs_statistics,1);
[sorted_vals, sorted_idx] = max(mean_power);
max_freqs = freq(sorted_idx); % 获取对应的频率值

text_x = 29; % right
text_y =max(permute_stat_95th)+0.01; % legend position
text(text_x,text_y,'p=0.05  corrected','FontSize',10,'FontName',myfontname,'Color', 'k','HorizontalAlignment', 'right','VerticalAlignment', 'bottom');
% %s1 = shadedErrorBar(faxis(interest_freq_index),mean(obs_statistics,1),std(obs_statistics,0,1) ./ sqrt(size(obs_statistics,1)),'lineProps',{'-','color',color(2,:)});
% %s1.mainLine.LineWidth = 1.5;
plot(faxis(interest_freq_index),mean(obs_statistics,1),'Color',color(2,:),'LineWidth', 1.5);
hold on
%plot(max_freqs, sorted_vals+0.02, '*', 'Color', color(3,:), 'MarkerSize', 5, 'LineWidth', 0.8);
title('FFT of proportion');
xlim([5 30]);
% ylim([0.05,0.22]);
% yticks(0.05:0.05:0.2);
% 确保 ylim 先更新
ylim auto;  
yl = ylim; % 获取最新的 y 轴范围
y_min = yl(1);
y_max = max(mean_power);
hold on;
line([max_freqs,max_freqs],[y_min,y_max], 'LineStyle','--','Color', 'k', 'LineWidth', 0.8);

xlabel('Frequency(Hz)');
ylabel('Power');
set(gca,'FontSize', myfontsize,'fontname',myfontname,'Box','off');
h=gcf;
box off;
h.PaperUnits='inch';
h.PaperSize = mypapersize;
h.PaperPosition=mypaperpos;