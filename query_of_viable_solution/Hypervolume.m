function Vol = Hypervolume(PopObj)
    % Calculate the Hypervolume of each solution in the same front
    
        [N,M]    = size(PopObj);
        Vol = ones(1,N);
        for i = 1 : M
            [~,rank] = sortrows(PopObj(:,i));
            Vol(rank(1))   = 0;
            Vol(rank(end)) = 0;
            for j = 2 : N
                Vol(rank(j)) = Vol(rank(j))*(PopObj(rank(j),i)-PopObj(rank(j-1),i));
            end
        end
  end         