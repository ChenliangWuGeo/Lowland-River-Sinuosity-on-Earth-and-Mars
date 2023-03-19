load batch1DischargeData.mat
load batch2DischargeData.mat
load batch2LMR.mat
%organize discharge data
typeA = [dischargeData1(1:6,:);dischargeData2(1:7,:)];
typeB = [dischargeData1(7:11,:);dischargeData2(8:10,:)];
discharge = [typeA;typeB];
%get migration rate
batch2B = [410,523,975,988,280,165,240,120,95,240];
Qp15a = [9.32E+00 9.8E+00 3.3 1.91E+00 10.4856,nan];
Qp15b = [3.11E+00 2.7537 4.254 10 6.11E+00];
migRateA = [9.8,3.067163,3.3,1.0324,8.18,nan];
migRateB = [5.951325, 9.0, 2.037274, 7.8e+01, 9.3];%
widthA = [1300,200,500,150,1250,nan];
widthB = [200,125,200,780,360];
migRateG1 = [migRateA,migRateB]./([widthA, widthB] .* [Qp15a,Qp15b]);
migMaxG1 = [13.92708355	3.994140855	4.2	nan	nan	nan 7.409461318	11.91003815	2.617478011	146.670545 11.3];
migMinG1 = [7.635198967	1.423298441	1.4	nan	nan	nan 4.094669859	5.645771564	1.457419709	37.53476435 6.3];
migRateG1max = migMaxG1./([widthA, widthB] .* [Qp15a,Qp15b]);
migRateG1min = migMinG1./([widthA, widthB] .* [Qp15a,Qp15b]);


migRateG2 = batch2LMR(:,3)'./batch2B./dischargeData2(:,5)';
migRateG2min = batch2LMR(:,1)'./batch2B./dischargeData2(:,5)';
migRateG2max = batch2LMR(:,2)'./batch2B./dischargeData2(:,5)';

%organize qsqw data
mig1 = [migRateG1(1:6),migRateG2(1:7)];
mig2 = [migRateG1(7:11),migRateG2(8:10)];
migR = [mig1,mig2];
qsqw1 = [3.02E-04	9.05E-05	2.54E-05	1.02E-04	4.32E-05...
    5.11E-04	1.47E-05	2.48E-04	3.24E-05	5.47E-05...
    3.12E-05	3.35E-06    3.20E-05];
qsqw2 = [6.81E-04	9.20E-03	7.62E-04	4.72E-01...
    6.27E-04	3.84E-05	9.70E-04    nan];
qsqw = [qsqw1,qsqw2];

figure(4);hold on
ytext = ["\itDV","\itQ_I"];
labelText = ["c","d"];
for i = 1:2
    subplot(1,4,i+2);hold on
    plotPara = 2*i;
    g = [zeros(length(typeA(:,plotPara)), 1);ones(length(typeB(:,plotPara)), 1)];
    vs = violinplot([typeA(:,plotPara);typeB(:,plotPara)],g,...
        'ShowMean', true,...
        'ViolinColor', [0,114,178;0,158,115]/256,...
        'ViolinAlpha',0.4);
    xticklabels({'variable \it\Omega','constant \it\Omega'});
    ylabel(ytext(i));
    text(.1, .9,labelText(i),'unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);
end

%qsqw
subplot(1,4,2)
g = [zeros(length(qsqw1), 1);ones(length(qsqw2), 1)];
vs = violinplot(log([(qsqw1)';(qsqw2)']),g,...
    'ShowMean', true,...
    'ViolinColor', [0,114,178;0,158,115]/256,...
    'ViolinAlpha',0.4);
xticklabels({'variable \it\Omega','constant \it\Omega'});
ylabel('log({\itQ_s/Q_w})');

xtext = ["\itDV","\itQ_I"];
text(.1, .9,'b','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);

%migration rate
subplot(1,4,1)
g = [zeros(length([migRateG1(1:6),migRateG2(1:7)]), 1);...
    ones(length([migRateG1(7:11),migRateG2(8:10)]), 1)];
vs = violinplot([log([migRateG1(1:6),migRateG2(1:7)]');...
    log([migRateG1(7:11),migRateG2(8:10)]')],g,...
    'ShowMean', true,...
    'ViolinColor', [0,114,178;0,158,115]/256,...
    'ViolinAlpha',0.4);
xticklabels({'variable \it\Omega','constant \it\Omega'});
ylabel('log({\itr*})');
xtext = ["\itDV","\itQ_I"];
text(.1, .9,'a','unit','normalized','VerticalAlignment','middle','HorizontalAlignment','left','fontsize',8);

set(gcf,'unit','inch','position',[1 1 12 3]);








