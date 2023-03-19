% calculate cutoff for most x domain, except for the last (noNode-k-3)
% points
function [x,y,T1,noNode] = cutOff(noNode,k,x,y,b,t,T1)
    for i = 1:(noNode-k-3)
        if i>noNode-k-3
            break
        end
        T1(t)=T1(t);

        p = [x(i+3:i+3+k);y(i+3:i+3+k)]';
        pq = [x(i),y(i)];

        [idx,dist] = dsearchn(p,pq);
        if dist < 2*b
            x = [x(1:i),x(i+4+idx:end)];
            y = [y(1:i),y(i+4+idx:end)];
            T1(t) = t;   
            clearvars p
        end

        [~,noNode] = size(x);

    end
end