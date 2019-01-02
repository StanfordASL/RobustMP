function path_red = triangle_reduce(path,obstacles)

n_wp = size(path,1);
path_red = path(1,:);

i_red = 1; %added nodes so far
i = 1; %base node to expand
while (i < n_wp) 
    next = i+1;
    if (next == n_wp)
        %auto add last node
        path_red(i_red+1,:) = path(next,:);
        break;
    end
    is_feas = 1;
    while (next <= n_wp) && (is_feas)
        if checkCollision(path(i,:),path(next,:),obstacles)
            if (next == n_wp) %no more further nodes
               path_red(i_red+1,:) = path(next,:);
               i = n_wp;
               break;
            else
                next = next + 1;
            end
        else
            is_feas = 0;
            %add furthest collision free node
            path_red(i_red+1,:) = path(next-1,:); 
            %augment count
            i_red = i_red + 1;
            %update expand node
            i = next-1;
        end
    end
end            

end