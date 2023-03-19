% calculate cutoff for the last (noNode-k-5)points
function [x,y,T1,noNode] = cutOff2(noNode,k,x,y,b,t,T1)
    for i = noNode-k-5:(noNode-5)
        if i>noNode-3
            break
        end
        T1(t)=T1(t);

        p = [x(i+3:end);y(i+3:end)]';
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