load('modernRivers.mat');%load modern river data
load('Result1000.mat');%load simulation results
noSimulation  = 1000;%number of simulation
tSpan = (0.1:0.01:10);%normalized tspan
[~,temp] = size(tSpan);
sinu = nan(noSimulation,temp);

for i = 1:noSimulation
    sinu(i,:) = Result(i).r.aveSinu;
end

aveSinu = mean(sinu,'omitnan');
sinu005 = prctile(sinu,.5);%.5 percentile of simulated sinuosity
sinu995 = prctile(sinu,99.5);%99.5 percentile of simulated sinuosity
Ti = interp1(aveSinu(tSpan<1),tSpan(tSpan<1),(1.001));%interpolate for ti when sinuosity is 1.001
tSpan = (tSpan-Ti)/(1-Ti);%rescale tSpan

aveSinu = aveSinu(tSpan>0);
sinu005 = sinu005(tSpan>0);
sinu995 = sinu995(tSpan>0);
tSpan = tSpan(tSpan>0);

aveSinu = aveSinu(~isnan(aveSinu));%filter nan
sinu005 = sinu005(~isnan(aveSinu));
sinu995 = sinu995(~isnan(aveSinu));
tSpan = tSpan(~isnan(aveSinu));

figure(3);hold on
p1 = fill([tSpan,fliplr(tSpan)],[sinu005,fliplr(sinu995)],[.8 .8 .8],'LineStyle','none');
p2 = plot(tSpan,aveSinu,'-k');
p3 = plot([1 1],[1 2.05],'k:');
h3 = text(1,2.2,'t^*_1','VerticalAlignment','top','HorizontalAlignment','center','fontsize',9);
p4 = plot([2.84 2.84],[1 2.05],'k:');
h4 = text(2.84,2.2,'t^*_s','VerticalAlignment','top','HorizontalAlignment','center','fontsize',9);

ylim([1 2.2]);
set(gca,'xscale','log');  
xlim([0.2 15]);
xticks([.2 .4 .6 .8 1 2 4 6]);
xticklabels({'0.2','0.4','0.6','0.8','1','2','4','6'});

xlabel('\itt*, t_a^*','fontsize',9);
ylabel('average sinuoisty \Omega','fontsize',9);
set(gca,'fontsize',9)

%plot group 1
% mississippi = rmfield(mississippi,{'dist','cumdist','Lrs','norm_toLb','sinu','slope'});
% trinity = rmfield(trinity,{'dist','cumdist','Lrs','norm_toLb','sinu','slope'});
% danube = rmfield(danube,{'dist','cumdist','Lrs','norm_toLb','sinu','slope'});
names = ["Mississippi","Trinity","Rhine","Meuse","Danube"];
grp1 = [mississippi; trinity; rhine; meuse; danube];
T1 = [744.380416130343,647.883880700573,566.563233487428;...
   572.864457337056,436.176512808747,358.312498897035;...
    [423.686763701032,365.968659885271,322.647111038905];...
    2220,2220,2220;...
    1524,1338,1190];%calculated t1 for variable sinuosity rivers

Ta = [1000,934,1000,1760,1570];%age of variable sinuosity rivers
for i = 1:5
    river = grp1(i);
    tempSinu = getSinu1(river,i);%mean, max,min, std
    ta = Ta(i);%
    t1 = T1(i,:);%max, mean, and min
    ta = ta./t1;
    taU = ta(3) - ta(2);
    taL = ta(2) - ta(1);
    sinuU = tempSinu(3) - tempSinu(1);
    sinuL = tempSinu(1) - tempSinu(2);
    if i == 5
        tempSinu = 1.9;%for the case of Danube
    elseif i == 3
        tempSinu = 1.38; %for the case of Rhine, use apex sinuosity
    end
    p5 = errorbar(ta(2),tempSinu(1),[],[],taL,taU,...
            'LineStyle','none','LineWidth',1,'color','k','Marker','o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0,114,178]/256,'MarkerSize',8,'CapSize',1);   
    txt = sprintf('%s',names(i));
    h = text(ta(2)+0.1,tempSinu(1),txt,'VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8);
    temp = h.Position;
    set(h,'Position',temp + [-0.05 0 0]);
    sinuG1t(i) = tempSinu(1);
    tG1t(i) = ta(2);
end

% plot lafourche
ta = 750./[1027.06154958642,836.427239119121,697.849601403565];
taU = ta(3) - ta(2);
taL = ta(2) - ta(1);

errorbar(750/836.427239119121,1.48,[],[],taL,taU,...
            'LineStyle','none','LineWidth',1,'color','k','Marker','o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0,114,178]/256,'MarkerSize',8,'CapSize',1); 
text(ta(2)+0.1,1.48,'Lafourche','VerticalAlignment','bottom','HorizontalAlignment','left','fontsize',8);

tG1 = [tG1t,ta(2)];%include lafourche
sinuG1 = [sinuG1t,1.4818];%include lafourche
tG1 = [tG1(1:3),tG1(5:6)];
sinuG1 = [sinuG1(1:3),sinuG1(5:6)];

load('sinuG1sim.mat')
sinuG1sim = [sinuG1sim(1:3),sinuG1sim(5:6)];

idx = [1,3,5];
a = tG1(idx);
b = sinuG1(idx);
c = sinuG1sim(idx);
% sinuG1sim(4) = 1.8364;%lafourch
arrow([a',b'+0.02],[a',c'-0.03],...
    'Length',7)
p7 = plot(tG1,sinuG1sim,'o','MarkerEdgeColor',[0,114,178]/256,...
            'MarkerFaceColor','none','MarkerSize',8,'linewidth',1.5);

%plot group 2
grp2 = [rioGrande; indus; brazos;colorado;chaoPhraya];%need to work on colorado
names = ["Rio Grande","Indus","Brazos","Colorado","Chao Phraya"];
T1 = [197.439934313268,165.758619937278,152.622308099659;...
    155.468320101498,128.358790977132,108.789099252873;...
    [215.904345140132,185.607961899178,168.280619776027];
    532.323719010047,452.655820237516,392.202201126703;...
    536,371,278];%t1 for constant sinuosity rivers
Ta = [1000,1000,2500,2500,2000];%age for constant sinuosity rivers
xtext = [5 0 0 0 0];%text location adjustment (x)
ytext = [0.05,0 0 0 0];%text location adjustment (y)
for i = 1:5
    river = grp2(i);
    tempSinu = getSinu2(river);
    ta = Ta(i);
    t1 = T1(i,:);%max, mean, and min
    ta = ta./t1;
    taU = ta(3) - ta(2);
    taL = ta(2) - ta(1);
    sinuU = tempSinu(3) - tempSinu(1);
    sinuL = tempSinu(1) - tempSinu(2);
    p6 = errorbar(ta(2),tempSinu(1),[],[],taL,taU,...
            'LineStyle','none','LineWidth',1,'color','k','Marker','d',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0,158,115]/256,'MarkerSize',9,'CapSize',1);
    txt = sprintf('%s',names(i));
    h = text(ta(2)+ xtext(i),tempSinu(1)+ytext(i),txt,'VerticalAlignment','top','HorizontalAlignment','right','fontsize',8);
    temp = h.Position;
    set(h,'Position',temp + [-0.1 0 0]);
    sinuG2(i) = tempSinu(1);
    tG2(i) = ta(2);
end

load('sinuG2sim.mat')
arrow([tG2(2:5)',sinuG2(2:5)'+0.03],[tG2(2:5)',sinuG2sim(2:5)'-0.03],...
    'Length',7)
p8 = plot(tG2,sinuG2sim,'d','MarkerEdgeColor',[0,158,115]/256,...
            'MarkerFaceColor','none','MarkerSize',9,'linewidth',1.5);
set(gca,'xscale','log');  
xlim([0.2 15]);
xticks([.2 .4 .6 .8 1 2 4 6]);
xticklabels({'0.2','0.4','0.6','0.8','1','2','4','6'});
ylabel('average sinuoisty \Omega','fontsize',9);
set(gca,'fontsize',9)
set(gca,'layer','top','XColor',[0 0 0],'YColor',[0 0 0],'LineWidth',1);
legend([p1,p2,p5,p6],{'.5-99.5%','average','variable \Omega','constant \Omega'},...
    'location','northwest');

