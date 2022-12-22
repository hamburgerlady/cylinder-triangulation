function R = getRfromV(V)

bestmin = inf;
Ve = zeros(3);
T = [1     0     0;...
     0     0     1;...
     0     1     0];

     
 
for iii = 1:2
    Ve(:,1) = (-1)^iii*V(:,1);
    for jjj = 1:2
        Ve(:,2) = (-1)^jjj*V(:,2);
        for kkk = 1:2
           Ve(:,3) = (-1)^kkk*V(:,3);
           Re = T*Ve';
           if det(Re)>0
               mini = norm(Re-eye(3));
               if mini<bestmin
                   bestmin = mini;
                   R = Re;
               end
           end
        end
    end
end



for iii = 1:2
    Ve(:,1) = (-1)^iii*V(:,2);
    for jjj = 1:2
        Ve(:,2) = (-1)^jjj*V(:,1);
        for kkk = 1:2
           Ve(:,3) = (-1)^kkk*V(:,3);
           Re = T*Ve';
           if det(Re)>0
               mini = norm(Re-eye(3));
               if mini<bestmin
                   bestmin = mini;
                   R = Re;
               end
           end
        end
    end
end


