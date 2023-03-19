clear
clc
load('MarsDelta.mat');
bin = ones(1,10);
bin = [bin,bin*2,bin*3,bin*4];
delta1 = [6,8,10];%index of seleted group 1 deltas
La1 = [ 22090,27710,22560];%avulsion length of group1 delta
delta2 = [5,7,9];%index of seleted group 2 deltas
La2 = [ 17800, 8058,24510];%avulsion length of group2 delta


for i = 1:3 
    deltaNo = delta1(i);
    temp = MarsDelta.r(deltaNo);
    aveSinu2(i,:) = boxAverageSinu(temp,La1(i));
    xStep = 0.1;
    xlimit = 4;
    xNormal = (xStep : xStep : xlimit);
end

bootStat = bootstrp(100,@mean,aveSinu2);
prc90 = prctile(bootStat,90);
prc10 = prctile(bootStat,10);

aveSinu = mean(aveSinu2,'omitnan');

for i = 1:3
    sinuLb = mean(aveSinu(bin==i),'omitnan');
    tempU = mean(prc90(bin==i),'omitnan');
    tempL = mean(prc10(bin==i),'omitnan');
    figure(2);hold on
    subplot(2,2,3);hold on
    p(2) = fill([i-1,i,i,i-1],[tempU,tempU,tempL,tempL],[.9 .9 .9],'LineStyle','none');
    p(1) = plot([i-1,i],[sinuLb,sinuLb],'-r');
end

figure(2);hold on
subplot(2,2,3);hold on
p(3) = stairs(xNormal,mean(aveSinu2,'omitnan'),'color',[0,114,178]/256,'LineWidth',1);
ylabel('binned average sinuosity \it\Omega');
ylim([1, 1.6]);
xlim([0, 2]);

xticks([0 1 2 3 4]);
xlabel('\itPD/L_a');
set(gca,'yscale','log')

text(.1,.9,'c   variable sinuosity rivers (Mars)','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
[H,p_value]=Mann_Kendall(aveSinu(~isnan(aveSinu)),0.05);

clearvars aveSinu2


for i = 1:2 
    deltaNo = delta2(i);
    temp = MarsDelta.r(deltaNo);
    aveSinu2(i,:) = boxAverageSinu(temp,La2(i));
    xStep = 0.1;
    xlimit = 4;
    xNormal = (xStep : xStep : xlimit);
end

bootStat = bootstrp(100,@mean,aveSinu2);
prc90 = prctile(bootStat,90);
prc10 = prctile(bootStat,10);
aveSinu = mean(aveSinu2,'omitnan');

for i = 1:3
    sinuLb = mean(aveSinu(bin==i),'omitnan');
    tempU = mean(prc90(bin==i),'omitnan');
    tempL = mean(prc10(bin==i),'omitnan');
    figure(2);hold on
    subplot(2,2,4);hold on
    p(2) = fill([i-1,i,i,i-1],[tempU,tempU,tempL,tempL],[.9 .9 .9],'LineStyle','none');
    p(1) = plot([i-1,i],[sinuLb,sinuLb],'-r');
end

figure(2);hold on
subplot(2,2,4);hold on
p(3) = stairs(xNormal,mean(aveSinu2,'omitnan'),'color',[0,158,115]/256,'LineWidth',1);
ylabel('binned average sinuosity \it\Omega');
ylim([1, 1.6]);
xlim([0, 2]);

xticks([0 1 2 3 4]);
xlabel('\itPD/L_a');
set(gca,'yscale','log')
text(.1,.9,'d   constant sinuosity rivers (Mars)','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
[H,p_value]=Mann_Kendall(aveSinu(~isnan(aveSinu)),0.05);
set(gcf,'unit','inch','position',[1 1 7 7]);

clearvars aveSinu2

function aveSinu2 = boxAverageSinu(delta,La)    
    %sinu2
    delta.xOver = delta.xOver(delta.sinu2>1.001);
    delta.sinu2 = delta.sinu2(delta.sinu2>1.001);
    tempxOver = delta.xOver/La;
    tempxOver = tempxOver*10;
    binLocation = floor(tempxOver);
    binLocation = binLocation/10;
    binLocation = binLocation(2:end);% resize, cause there are n+1 RKOver and n sinu
    m = 1;
    xlimit = 4;
    xStep = 0.1;
    aveSinu2 = nan(size(xStep:xStep:xlimit));
    upperLimit = nan(size(xStep:xStep:xlimit));
    lowerLimit = nan(size(xStep:xStep:xlimit));
    for i = xStep : xStep : xlimit
        temp = delta.sinu2(binLocation>i-0.1 & binLocation<=i);
        aveSinu2(m) = mean(temp);
        upperLimit(m) = prctile(temp,75) - aveSinu2(m);
        lowerLimit(m) = aveSinu2(m) - prctile(temp,25);
        m = m + 1;
    end
end