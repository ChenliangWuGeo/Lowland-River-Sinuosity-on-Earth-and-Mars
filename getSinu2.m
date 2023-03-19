function output = getSinu2(river)
    sinu = river.Data.sinu2;
    Lb = river.Data.Lb;
    RKOver = river.Data.RKOver;
    RKOver = RKOver/1e3/Lb;
    RKOver = RKOver(2:end);
    output(1) = mean(sinu(RKOver<1));
%     output(2) = output(1) - std(sinu(RKOver<1));
%     output(3) = output(1) + std(sinu(RKOver<1));
    output(2) = prctile(sinu(RKOver<1),25);
    output(3) = prctile(sinu(RKOver<1),75);



end