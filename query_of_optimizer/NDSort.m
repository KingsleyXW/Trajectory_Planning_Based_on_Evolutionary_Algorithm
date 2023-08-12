function [FrontNO,MaxFNO] = NDSort(PopObj,nSort)
    %NDSort - Do non-dominated sorting on the population
    %
    %   FrontNO = NDSort(A,s) does non-dominated sorting on A, where A is a
    %   matrix which stores the objective values of all the individuals in the
    %   population, and s is the number of individuals being sorted at least.
    %   FrontNO(i) means the number of front of the i-th individual.
    %
    %   [FrontNO,K] = NDSort(...) also returns the maximum number of fronts,
    %   except for the value of inf.
    %
    %   In particular, s = 1 stands for find only the first non-dominated
    %   front, s = size(A,1)/2 stands for sort only half of the population
    %   (which is often used in the algorithm), and s = inf stands for sort the
    %   whole population.
    %
    %   Example:
    %   [FrontNO,MaxFNO] = NDSort(PopObj,1)
    
        [N,M] = size(PopObj);
        
        FrontNO = inf(1,N);
        MaxFNO  = 0;
        [PopObj,rank] = sortrows(PopObj);
        while sum(FrontNO<inf) < min(nSort,N)
            MaxFNO = MaxFNO + 1;
            for i = 1 : N
                if FrontNO(i) == inf
                    Dominated = false;
                    for j = i-1 : -1 : 1
                        if FrontNO(j) == MaxFNO
                            m = 2;
                            while m <= M && PopObj(i,m) >= PopObj(j,m)
                                m = m + 1;
                            end
                            Dominated = m > M;
                            if Dominated || M == 2
                                break;
                            end
                        end
                    end
                    if ~Dominated
                        FrontNO(i) = MaxFNO;
                    end
                end
            end
        end
        FrontNO(rank) = FrontNO;
    end
    