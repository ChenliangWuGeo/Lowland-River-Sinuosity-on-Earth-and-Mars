%calculate sinuosity as ratio between alongstream distance and cartisian
%distance between two succeeding cross over points
function output = calculateSinu2(curvature_t,x_t,y_t,accuS_t)

    [pks1, loca1]=findpeaks(curvature_t);
    [pks2, loca2]=findpeaks(-curvature_t);

    locaAll = [loca1,loca2];
    locaAll = sort(locaAll);

    [~,lgth] = size(locaAll);
    if mod(lgth,2) ~= 0 %remainder after divided by 2, in case of odd number
        locaAll = locaAll(1:end-1);
        lgth = lgth - 1;
    end

    xOver = nan(1,lgth);
    yOver = nan(1,lgth);
    accuSOver = nan(1,lgth);

    intrinsicLength = nan(1,lgth-1);
    directLength = nan(1,lgth-1);

    for i = 1  : lgth-1
        locaTemp = locaAll(i):locaAll(i+1);
        xTemp = x_t(locaTemp);
        yTemp = y_t(locaTemp);
        accuSTemp = accuS_t(locaTemp);
        curvTemp = curvature_t(locaTemp);
        xOver(i) = interp1(curvTemp,xTemp,[0]);
        yOver(i) = interp1(curvTemp,yTemp,[0]);
        accuSOver(i) = interp1(curvTemp,accuSTemp,[0]);
    end

    xOver = xOver(~isnan(xOver));
    yOver = yOver(~isnan(yOver));
    accuSOver = accuSOver(~isnan(accuSOver));

    intrinsicLength = accuSOver(2:end)-accuSOver(1:end-1);
    directLength = sqrt(...
        (xOver(2:end)-xOver(1:end-1)).^2 +...
        (yOver(2:end)-yOver(1:end-1)).^2);

    sinu2 = intrinsicLength./directLength;
    
    output = sinu2;
end

