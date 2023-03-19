function output = getSinu1(river,i)
    sinu = river.Data.sinu2;
    Lb = river.Data.Lb;
    RKOver = river.Data.RKOver;
    RKOver = RKOver/1e3/Lb;
    RKOver = RKOver(2:end);
%     if i == 3
%         temp = sinu(RKOver>1.5 & RKOver<2);
%     else
        temp = sinu(RKOver>0.67 & RKOver<1);
%     end
%     temp = temp(temp>1.05);
    output(1) = mean(temp);
%     output(2) = output(1) - std(sinu(RKOver>0.67 & RKOver<1));
%     output(3) = output(1) + std(sinu(RKOver>0.67 & RKOver<1));
    output(2) = prctile(temp,25);
    output(3) = prctile(temp,75);
    output(4) = mean(bootstrp(1000,@(x)(std(x)),temp));  
%     histogram(bootstrp(1000,@(x)(std(x)),temp));
%     keyboard
end