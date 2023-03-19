    function output = calculateSinu(x_t,y_t,xStep)
    n = 20;
    Lpath = n * xStep;
    Ldist = sqrt( (x_t(n+1:end) - x_t(1:end-n)).^2+...
        (y_t(n+1:end)- y_t(1:end-n)).^2);
    sinu = Lpath./Ldist;
    output = sinu;
    end