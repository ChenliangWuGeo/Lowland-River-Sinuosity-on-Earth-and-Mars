function output = simuFunction(hInput,bInput,uInput,CfInput,EInput,sInput)
    Cf = CfInput;
    rng('default')
    year = 365*24*60*60;
    noNode = 200;%number of modeling node
    sLog = -log10(sInput);
    
    noYear = round (1000 + 4000*(sLog-3)/3);%longer run time for lower slope channels
    tf = 1;%time scaling factor
    noTime = noYear/tf;%modeling time steps
    T1 = nan(1,noTime);%first cut off time
    xStep = bInput*3;
    xlimit = noNode*xStep;
    x = xStep:xStep:xlimit;
    y = zeros(1,noNode) + rand(1,noNode)*bInput/10;
    xNor = x/x(end)*25*pi;
    xi = x;

    X = nan(noYear,1000);
    Y = nan(noYear,1000);
    SINU = nan(noYear,1000);

    us0 = uInput; %normal flow velocity
    h0 = hInput; %normal flow depth
    b = bInput;%half channel width
    g = 9.81;%gravitational acceleration
    A = 10; %scouring factor
    alfa = 0.077;
    Chi1 = alfa/sqrt(Cf);
    Chi = Chi1 - 1/3;
    As = 181 * (h0/b)^2 * 1/Chi1 * (2*Chi^2 + 0.8*Chi + 1/15);
    A = A + As -1;
    E = EInput;
    k = 80; %searching index range for cutoff point

    curvature = nan(1,noNode);
    dCurve = nan(1,noNode);
    s = nan(1,noNode);
    accuS = nan(1,noNode);
    ds = nan(1,noNode);
    dx1 = nan(1,noNode);
    dx2 = nan(1,noNode);
    dy1 = nan(1,noNode);
    dy2 = nan(1,noNode);
    usb = zeros(1,noNode);
    migrRate = zeros(1,noNode);
    maxSinu = zeros(1,noYear);
    aveSinu = zeros(1,noYear);
    maxSinu2 = zeros(1,noYear);
    aveSinu2 = zeros(1,noYear);
    aveMigRate = zeros(1,noYear);
    maxMigRate = zeros(1,noYear);
    aveCurve = zeros(1,noYear);
    maxCurve = zeros(1,noYear);

    sinu25 = zeros(1,noYear);
    sinu75 = zeros(1,noYear);


    for t=1:noTime
    ds(2:end) = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2);
    ds(1) = ds(2);

    dx1(2:end) =  (x(2:end)-x(1:end-1))./ds(2:end);
    dx1(1) = dx1(2);
    dx2(2:end) = (dx1(2:end)-dx1(1:end-1))./ds(2:end);
    dx2(1) = dx2(2);

    dy1(2:end) = (y(2:end)-y(1:end-1))./ds(2:end);
    dy1(1) = dy1(2);
    dy2(2:end) = (dy1(2:end)-dy1(1:end-1))./ds(2:end);
    dy2(1) = dy2(2);

    curvature = (dx1.*dy2-dy1.*dx2)./(dx1.^2+dy1.^2).^(3/2);
    dCurve(2:end) = (curvature(2:end)-curvature(1:end-1))./ds(2:end);
    dCurve(1) = dCurve(2);
    dnx = (-dy1)./sqrt(dx1.^2+dy1.^2);
    dny = dx1./sqrt(dx1.^2+dy1.^2);

    usb(2:end) = b./(us0./ds(2:end) + 2*us0./h0*Cf) .*...
        (-us0.^2.*dCurve(2:end) + ...
        Cf*curvature(2:end) .* ((us0.^4/g./h0.^2) + A*us0.^2./h0)...
        + us0./ds(2:end) .* usb(1:end-1)/b);
    usb(1) = 0;

    migrRate = E.*usb*year*tf;
    L = sqrt(dnx.^2+dny.^2);
    dx = -dnx.*migrRate;
    dy = -dny.*migrRate;

    dx(end) = 0;
    dy(end) = 0;

    x = x + dx;
    y = y + dy;

    %save x, curvature and s from previous time step,
    x_t = x;
    y_t = y;
    curvature_t = curvature;

    accuS = streamLineDistance(ds);
    accuS_t = accuS;%
    maxCurve(t) = max(abs(curvature_t));   
    aveCurve(t) = mean(abs(curvature_t));
    aveMigRate(t) = mean(abs(migrRate));
    maxMigRate(t) = max(abs(migrRate));

    %calculate cutoff
    [x,y,T1,noNode]=cutOff(noNode,k,x,y,b,t,T1);
    %calculate cutoff for the most downstream section
    [x,y,T1,noNode]=cutOff2(noNode,k,x,y,b,t,T1);

    [~,noNode] = size(x);
    ds = nan(1,noNode);
    ds(2:end) = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2);
    ds(1) = 0;

    accuS = streamLineDistance(ds);
    s_reMesh = accuS(1):xStep:accuS(end);
    x_reMesh = interp1(accuS,x,s_reMesh,'spline');%'makima'
    y_reMesh = interp1(accuS,y,s_reMesh,'spline');
    x_reMesh = [x_reMesh, xlimit];%add the end point, so model is stable
    y_reMesh = [y_reMesh, 0];
    x = x_reMesh;
    y = y_reMesh;


    y(1) = y(1) + (rand(1,1)-0.5)*1;% small perturbation at the upstream, without this, upstream become straight    
    y(end) = y(end-1) + (rand(1,1)-0.5)*1;

    [~,noNode] = size(x);
    us0_t = us0;
    usb_t = usb;
    migrRate_t = migrRate;
    curvature = nan(1,noNode);
    dCurve = nan(1,noNode);
    s = nan(1,noNode);
    accuS = nan(1,noNode);
    ds = nan(1,noNode);
    dx1 = nan(1,noNode);
    dx2 = nan(1,noNode);
    dy1 = nan(1,noNode);
    dy2 = nan(1,noNode);
    usb = zeros(1,noNode);
    migrRate = zeros(1,noNode);

    sinu = calculateSinu(x_t,y_t,xStep);
    maxSinu(t) = max(sinu);
    aveSinu(t) = mean(sinu);
    [~,L] = size(sinu);

    %calculate intrinsic sinuosity
    sinu2 = calculateSinu2(curvature_t,x_t,y_t,accuS_t);
    maxSinu2(t) = max(sinu2);
    aveSinu2(t) = mean(sinu2);
    sinu25(t) =   prctile(sinu2,25);
    sinu75(t) =   prctile(sinu2,75);

    %store sinuosity through time
    X(t,1:L) = x_t(1:L);
    Y(t,1:L) = y_t(1:L);
    SINU(t,1:L) = sinu;
    end

    T1 = T1(~isnan(T1));
    T1 = T1(1);
    tSpan = (0.1:0.01:10);
    tNormal = (1:noYear)/T1;
    aveSinuInterp = interp1(tNormal,aveSinu2,tSpan);
    sinu25Interp = interp1(tNormal,sinu25,tSpan);
    sinu75Interp = interp1(tNormal,sinu75,tSpan);

    result.T1 = T1;
    result.aveSinu = aveSinuInterp;
    result.sinu25 = sinu25Interp;
    result.sinu75 = sinu75Interp;

    output = result;

end

