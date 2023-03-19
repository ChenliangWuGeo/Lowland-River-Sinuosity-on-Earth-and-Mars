clear 
clc
load('modernRivers.mat')
riverGroup1 = [mississippi,trinity,rhine,meuse,danube,nile,don,kobuk,yana,mackenzie,tombigbee,alabama,apalachicola];
riverGroup2 = [brazos,rioGrande,colorado,indus,chaoPhraya,sabine,neches,suwanee];

bin = ones(1,10);
bin = [bin,bin*2,bin*3,bin*4];

%plot group1
for i = 1:13 %6 is number of rivers in group 1   
    aveSinu2(i,:) = boxAverageSinu(riverGroup1(i));
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
    subplot(2,2,1);hold on
    p(2) = fill([i-1,i,i,i-1],[tempU,tempU,tempL,tempL],[.9 .9 .9],'LineStyle','none');
    p(1) = plot([i-1,i],[sinuLb,sinuLb],'-r');
end

figure(2);hold on
subplot(2,2,1);hold on
p(3) = stairs(xNormal,mean(aveSinu2,'omitnan'),'color',[0,114,178]/256,'LineWidth',1);
ylabel('binned average sinuosity \it\Omega');
ylim([1, 1.6]);
xlim([0, 3]);

xticks([0 1 2 3 4]);
xlabel('\itRK/L_a');

legend([p(3),p(1),p(2)],{'bin size 0.1\it{L_a}','bin size \it{L_a}','10-90%'},...
    'box','off','location','northeast');
text(.1,.9,'a','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);
text(1.5,1.1,'variable sinuosity rivers (Earth)','VerticalAlignment','middle','HorizontalAlignment','center','fontsize',8);
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
set(gca,'yscale','log')

%plot group2
[H,p_value]=Mann_Kendall(aveSinu(1:30),0.05);
clearvars aveSinu2

for i = 1:8 %5 is number of rivers in group 1      
    aveSinu2(i,:) = boxAverageSinu(riverGroup2(i));
    
    xStep = 0.1;
    xlimit = 4;
    xNormal = (xStep : xStep : xlimit);
end

bootStat = bootstrp(100,@mean,aveSinu2);
prc90 = prctile(bootStat,90);
prc10 = prctile(bootStat,10);

aveSinu = mean(aveSinu2,'omitnan');

for i = 1:3
    sinuLb = mean(aveSinu(bin==i));
    tempU = mean(prc90(bin==i));
    tempL = mean(prc10(bin==i));
    figure(2);hold on
    subplot(2,2,2);hold on
    fill([i-1,i,i,i-1],[tempU,tempU,tempL,tempL],[.9 .9 .9],'LineStyle','none');
    
    plot([i-1,i],[sinuLb,sinuLb],'-r');
end

figure(2);hold on
subplot(2,2,2);hold on
p(4) = stairs(xNormal,mean(aveSinu2,'omitnan'),'color',[0,158,115]/256,'LineWidth',1);
ylabel('binned average sinuosity \it\Omega');
ylim([1, 1.6]);
xlim([0, 3]);
xticks([0 1 2 3 4]);
xlabel('\itRK/L_a');
legend([p(4)],{'bin size 0.1\it{L_a}'},...
    'box','off','location','northeast');
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
text(.1,.9,'b','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);
text(1.5,1.1,'constant sinuosity rivers (Earth)','VerticalAlignment','middle','HorizontalAlignment','center','fontsize',8);
set(gca,'yscale','log')

% entire figure
set(gcf,'unit','inch','position',[1 1 7 7]);

[H,p_value]=Mann_Kendall(aveSinu(1:30),0.05);

function aveSinu2 = boxAverageSinu(river)    
    %sinu2
    tempRKOver = river.Data.RKOver/river.Lb/1e3;
    tempRKOver = tempRKOver*10;
    binLocation = floor(tempRKOver);
    binLocation = binLocation/10;
    binLocation = binLocation(2:end);% resize, cause there are n+1 RKOver and n sinu
    m = 1;
    xlimit = 4;
    xStep = 0.1;
    aveSinu2 = nan(size(xStep:xStep:xlimit));
    upperLimit = nan(size(xStep:xStep:xlimit));
    lowerLimit = nan(size(xStep:xStep:xlimit));
    for i = xStep : xStep : xlimit
        temp = river.Data.sinu2(binLocation>i-0.1 & binLocation<=i);
        aveSinu2(m) = mean(temp);
        upperLimit(m) = prctile(temp,75) - aveSinu2(m);
        lowerLimit(m) = aveSinu2(m) - prctile(temp,25);
        m = m + 1;
    end
end
