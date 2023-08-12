% AMOACO to enhance the road trajectories of cutting head
% to form high quality Section in Coal Mine Roadway
% .....................................by Kingsley W

clear all
close all
clc
tic
% create the gripmap
a = ones(16); % changeable 
% a(3,5:6)=0;
% a(4,5:6)=0;
% a(3,7:8)=0;
% a(4,7:8)=0;
% a(5,5:8)=0;
% a(6,5:8)=0;
%  
% a(10:15,13:16)=0;
% a(11:14,5:8)=0;

% a(7:10,7:10)=0;

% a(7:10,1:4)=0;
% a(7:10,13:16)=0;

a(7:10,1:2)=0;
a(7:10,15:16)=0;
a(7:10,7:10)=0;
a(13:16,7:10)=0;

% a(6:8,5:8)=0;
% a(10:15,13:16)=0;
% a(11:14,5:8)=0;

% a(17:20,9:12)=0;
% a(9:12,1:4)=0;
% a(9:12,17:20)=0;
% a(7:10,9:12)=0;

b=a;
b(end+1,end+1)=0;
a=b
colormap([0 0 0;1 1 1]);  % use colormap to create the wanted color
X_grid=[1:size(a,2)].*ones(size(a,1),1);
Y_grid=X_grid';

%create the decision variable
Nx=size(a,2)-1;
Ny=size(a,1)-1;
n=Nx*Ny/4;
n1=Nx/2;
n2=Ny/2;
XY=zeros(n2,n1);
coordinatesXY=zeros(n2*n1,2);
S=zeros(1,n2*n1); % coordinates sequence
R=0.3; % cutting radium
ratio=1/R; % cutting width unit/head diameter

x_grid=(X_grid(1,:)-1)/ratio; 
y_grid=(Y_grid(:,1)-1)/ratio;  
x=x_grid(2:2:[length(x_grid)-1]); % x coordinates when ratio=2 x=[0.5:n1-0.5]; 
y=y_grid(2:2:[length(y_grid)-1]); % y coordinates when ratio=2 y=[0.5:n2-0.5]; 
pcolor(x_grid,y_grid,a);
set(gca,'XTick',x_grid,'YTick',y_grid);  % axis setting
xlabel('Width');
ylabel('Height');
title('Coal Mine Section Surface')
axis image xy; 
hold on


% self version coordinates index to variable XY index S 
k=1;
for i=1:n2
     for j=1:n1
         XY(i,j)=k;
         S(k)=k;
         coordinatesXY(k,:)=[x(j),y(i)];
         plot(x(j),y(i),'ro')
         labl=num2str(k);
         labl=[' ',labl];
         text(x(j),y(i),labl)
         hold on
         k=k+1;
     end
end
hold off
Edge=[XY(1,:)',XY(:,1),XY(:,length(XY)),XY(length(XY),:)']; % define edge matrix
Edge(end,1)=0; Edge(1,4)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % objective function parameter setting
    Hw=6;% Specific energy consumption of cutting with different hardness of coal rock
    A=R*(size(a,2)-1); % road width
    B=R*(size(a,1)-1); % road height
    h=0.2; % average cutting depth
    R;
    %S2=9.25; % block area
    S0=(R^2-pi*R^2/4)*2; % unit left area
    S1=4*R^2;
    K=0.6; %Kjol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % aco algorithm parameter
    Ite = 250; % max iteration
    Ant_num = 100;
    Alpha = 1;
    Rho = 0.15;
    Q = 100;
    Beta = 5;

    City=coordinatesXY;
    City_num = size(City,1);%n indicate the number of the nodes(cities)
    Distance = zeros(City_num,City_num);% adjunct matrix
   
    % caculate the adjunct distance method 1
    for i=1:n
        for j=i+1:n
            Distance(i,j)=sqrt(sum((City(i,:)-City(j,:)).^2));
            Distance(j,i)=Distance(i,j);
        end
    end
%     for i=1:City_num % caculate the adjunct distance method 2
%         for j=1:City_num
%             if i~=j
%                 Distance(i,j)= sqrt((City(i,1)-City(j,1))^2 + (City(i,2)-City(j,2))^2);
%             else
%                 Distance(i,j) = eps;
%             end
%         end
%     end
%     Distance = abs(Distance+Distance');

    % presetting, initiating the algorithm parameter
    Eta = 1./Distance; 
    Tau = ones(City_num,City_num);   
    c=20; % constant hermone distribution value
    for i=1:City_num 
        if ismember([i],[length(x)*[1:length(y)]])==1
            route=[i-1];
        elseif ismember([i],[length(x)*[1:length(y)-1]+1])==1
            route=[i+1];
        else
            route=[i-1,i+1];
        end
        for j=1:City_num 
            if a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))~=0 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))~=0 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))~=0 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1),1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))~=0
               if ismember([j],route)==1
                  Tau(i,j)=c/Distance(i,j);
               end
            end
        end
    end
    for i=1:City_num % presetting hormone Tau
        if ismember([i],[length(x)*[2:length(y)-1]])==1
           judge=[i-1,0
                  i+length(x),i-length(x)];
                  %0,i+length(x)-1
                  %0,i-length(x)-1];
        elseif ismember([i],[length(x)*[1:length(y)]+1])==1
           judge=[0,i+1
           i+length(x),i-length(x)];
           %i+length(x)+1,0
           %i-length(x)+1,0];
        elseif ismember([i],[1])==1
            judge=[i+1,0
           i+length(x),0];
           %i+length(x)+1,0
           %0,0];
        elseif ismember([i],[length(x)])==1
            judge=[i-1,0
           i+length(x),0];
           %0,i+length(x)-1
           %0,0];
           
        elseif ismember([i],[(length(x)-1)*length(y)+1])==1
            judge=[i+1,0
           0,i-length(x)];
           %i-length(x)+1,0
          % 0,0];
        elseif ismember([i],[length(x)*length(y)])==1
            judge=[0,i-1
           0,i-length(x)];
          % 0,i-length(x)-1
          % 0,0];
        else
            judge=[i-1,i+1
            i+length(x),i-length(x)];
            %i+length(x)+1,i-length(x)-1
            %i-length(x)+1,i+length(x)-1];
        end
        for j=1:City_num
            if ismember([i,j],Edge)==1 & ismember([j],judge)==1
                Tau(i,j)=1.8*c/Distance(i,j);
            elseif a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==0 || a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1),1),X_grid(1,ratio*coordinatesXY(j,1)+1))==0 
                Tau(i,j)=eps;
            elseif a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==0 || a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1),1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==0
                Tau(i,j)=eps;
            elseif ismember([j],judge)==1 
                Tau(i,j)=0.5*c/Distance(i,j);
            else
                Tau(i,j)=eps/Distance(i,j);
            end
        end
    end

    

    % remove the city that is in the barricades
    City_num_temp=[1:City_num];
    City_num=City_num_temp;
    for j=1:length(City_num)
        if a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==0 || a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==0 
             City_num_temp(1,j)=0; 
        elseif a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==0 || a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==0
             City_num_temp(1,j)=0;
        end
    end
     [r0,c0]=find(City_num_temp==0);
     block_num=length(c0);
     S2=block_num*S1; % block area
     City_num_temp(:,c0)=[]; 
     City_num=City_num_temp;
     block_num=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACO progress parameter defining 
Ite_num = 1; % iteration time (generation time)
R_rec = zeros(Ant_num,length(City_num));
R_rec_temp = zeros(Ant_num,length(City_num));
R_best = zeros(Ite,length(City_num)); % best road of every generation 
L_best = inf.*ones(Ite,1); % best length of best road of each generation 
L_ave = zeros(Ite,1); % best length of best road of each generation   
archived_solution=zeros(Ite,2); 
L_best1=zeros(Ite,1);
 % iteration parameter preallocation
archived_phermone=zeros(Ant_num,2);
best_routes=zeros(Ant_num,length(City_num));
bad_ant=0;
Index=0;

%ACO begins
while Ite_num<=Ite
    % put 'ant_num' ants into 'city_num' cities randomly
        Init_ant_R_rec = []; % A temporary matrix that randomly assigns the ants in the initial state to each city 
        for i = 1:(ceil(Ant_num/length(City_num)))
            Init_ant_R_rec = [Init_ant_R_rec,ones(1,(length(City_num)))]; % The first city of each ants
        end

        R_rec(:,1) = (Init_ant_R_rec(1,1:Ant_num))'; % After the matrix is transposed, 
        R_rec_temp= R_rec(:,1);                         % it becomes an initial matrix of Ant_num row-1 column, each row represents an ant, 
                                                        % and each value represents the city code of the current ant. 
    % choose next city using probability function
    for j = 2:length(City_num) % city loop, current city not included i=city_num
        for i = 1:Ant_num % for each ant
            City_visited = R_rec(i,1:(j-1)); % tatu table
            City_unvisited = zeros(1,(length(City_num)-j+1)); 
            %P = eps*(ones(1,(length(City_num)-j+1))); % prealloacting probability function
            count = 1;
            
            city_current=City_visited(end);
            if ismember([city_current],length(x)*[2:length(y)-1])==1
                judge=[city_current-1,0
                       city_current+length(x),city_current-length(x)
                       0,city_current+length(x)-1
                       0,city_current-length(x)-1];
            elseif ismember([city_current],[length(x)*[1:length(y)]+1])==1
                judge=[city_current+1,0
                       city_current+length(x),city_current-length(x)
                       city_current+length(x)+1,0
                       city_current-length(x)+1,0];
            elseif ismember([city_current],[1])==1
                judge=[city_current+1,0
                      city_current+length(x),0
                      city_current+length(x)+1,0
                      0,0];
            elseif ismember([city_current],[length(x)])==1
                judge=[city_current-1,0
                      city_current+length(x),0
                      0,city_current+length(x)-1
                      0,0];
               
            elseif ismember([city_current],[(length(x)-1)*length(y)+1])==1
                judge=[city_current+1,0
                      0,city_current-length(x)
                      city_current-length(x)+1,0
                      0,0];
            elseif ismember([city_current],[length(x)*length(y)])==1
                judge=[0,city_current-1
                       0,city_current-length(x)
                       0,city_current-length(x)-1
                       0,0];
            else
                 judge=[city_current-1,city_current+1
                 city_current+length(x),city_current-length(x)
                 city_current+length(x)+1,city_current+length(x)-1
                 city_current-length(x)+1,city_current-length(x)-1];
            end
            % choose the city unvisited
            for k =City_num(2:end)
                %if isempty(find(City_visited == k,1)) % avoid going repeated city
                if ismember([k],judge)==1 && isempty(find(City_visited == k,1)) ...
                && a(Y_grid(ratio*coordinatesXY(k,2)+1,1),X_grid(1,ratio*coordinatesXY(k,1)+1))==1 && a(Y_grid(round(ratio*(coordinatesXY(k,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(k,1)+1))==1 ... 
                && a(Y_grid(ratio*coordinatesXY(k,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(k,1)-1/ratio))+1))==1 && a(Y_grid(round(ratio*(coordinatesXY(k,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(k,1)-1/ratio)+1)))==1
                   City_unvisited(count) = k;
                   count = count+1;
                end
            end
            [r,c]=find(City_unvisited==0);
            City_unvisited(:,c)=[];         
            if length(City_unvisited)==0 && length(City_num)-length(City_visited)<=block_num && City_visited(end)-City_visited(end-1)~=0  % need to be reviesed
               City_unvisited = zeros(1,(length(City_num)-j+1));
               for kb =c0(1:end) 
                   City_unvisited(count)=kb; % start cutting block rocks
                   count = count+1;
               end
               [r,c]=find(City_unvisited==0);
               City_unvisited(:,c)=[]; 
            elseif length(City_unvisited)==0 & (length(City_num)-length(City_visited)>block_num | City_visited(end)-City_visited(end-1)==0)
                   bad_ant=i; % record the the ant unable to cover all city
                   City_unvisited=city_current;
            end
            % caculate the probility of city to be visited
            P = eps*(ones(1,length(City_unvisited)));
            for kk = 1:length(City_unvisited)
                P(kk) = (Tau(City_visited(end),City_unvisited(kk))^Alpha)*...
                (Eta(City_visited(end),City_unvisited(kk))^Beta);           
                if P(kk)==0
                   P(kk)=eps;
                end
            end
            P=P/sum(P);

            % wheel selection to choose next city
            P_cum=cumsum(P);
            Select = find(P_cum>=rand);
            if length(City_unvisited)==1
               To_visit=City_unvisited;
            else
                To_visit = City_unvisited(Select(1));
            end

            R_rec(i,j) = To_visit; % The i-th ant will go to the city where the j-th step will go. After the cycle is over, 

            if i==bad_ant
            R_rec_temp(i,j)=0;
            else
            R_rec_temp(i,j)=To_visit;   
            end                      % the path of each ant will be obtained
        end
    end
    
   % fix solution methods
   nn=1;
   propor=6;
   while nn>=1
    fix_range=ceil(length(City_num)/propor);
    for i = 1:Ant_num
        [r,c]=find(R_rec_temp(i,:)==0);
        rec_temp=R_rec_temp(i,:);
        rec_temp(c)=[];
        if length(rec_temp)>=length(City_num)-fix_range
           [logic,lable]=ismember(City_num,rec_temp);
           [rr,cc]=find(lable==0);
           City_left=City_num(cc);
           for j=City_left(1:length(cc))
               if ismember(j,Edge)==0 & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio))+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==1
                 City_adjunct=[j-1,j+1,j+Nx/2,j-Nx/2];
               elseif ismember(j,Edge(:,1))==1 & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio))+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==1
                 City_adjunct=[j-1,j+1,j+Nx/2];
               elseif ismember(j,Edge(:,2))==1 & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio))+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==1
                 City_adjunct=[j+1,j+Nx/2,j-Nx/2]; 
               elseif ismember(j,Edge(:,3))==1 & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio))+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==1
                 City_adjunct=[j-1,j+Nx/2,j-Nx/2];
               elseif ismember(j,Edge(:,4))==1 & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,ratio*coordinatesXY(j,1)+1))==1 ... 
               & a(Y_grid(ratio*coordinatesXY(j,2)+1,1),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio))+1))==1 & a(Y_grid(round(ratio*(coordinatesXY(j,2)-1/ratio)+1,1)),X_grid(1,round(ratio*(coordinatesXY(j,1)-1/ratio)+1)))==1
                 City_adjunct=[j-1,j+1,j-Nx/2];
               end
               for kk=1:length(City_adjunct)
                   [rrr,ccc]=find(rec_temp==City_adjunct(kk));
                   if length(ccc)~=0 & (City_adjunct(kk)~=1 | (City_adjunct(kk)==1 & length(rec_temp(1:ccc(1)-1))~=0))
                      rec_temp=[rec_temp(1:ccc(1)-1),j,rec_temp(ccc(1):end)];
                      break
                   end 
               end
            end
        end
        R_rec_temp(i,:)=[rec_temp,zeros(1,length(City_num)-length(rec_temp))];
        R_rec(i,:)=[rec_temp,rec_temp(1,end)*ones(1,(length(City_num)-length(rec_temp)))];
    end
    nn=nn-1;
   end   
    Len = zeros(Ant_num,1);
    for i=1:Ant_num
           R_temp = R_rec(i,:);
           for j = 1:(length(City_num)-1)
               Len(i) = Len(i) + Distance(R_temp(j),R_temp(j+1));  
           end
           %Len(i)=Len(i)+Distance(R_temp(1),R_temp(City_num)); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multi-objective function defining 
    N=zeros(Ant_num,length(City_num));
    N_all=zeros(Ant_num,1);
    f=zeros(Ant_num,2); % multi-objective function
    w1=0.1;w2=0.9; %function weight
    u=zeros(Ant_num,2); % satisfaction rate matrix
    U=zeros(Ant_num,1);
    for i = 1 : Ant_num
        for j=2:length(City_num)-1
            if ismember([R_rec(i,j)],Edge)==1 & R_rec(i,j+1)~= R_rec(i,j)
               if (abs(R_rec(i,j+1)-R_rec(i,j))==1 & abs(R_rec(i,j)-R_rec(i,j-1))==1) | (abs(R_rec(i,j+1)-R_rec(i,j))==Nx/2 & abs(R_rec(i,j)-R_rec(i,j-1))==Nx/2) ...
                | (abs(R_rec(i,j+1)-R_rec(i,j))==1 & abs(R_rec(i,j)-R_rec(i,j-1))==Nx/2) | (abs(R_rec(i,j+1)-R_rec(i,j))==Nx/2 & abs(R_rec(i,j)-R_rec(i,j-1))==1)
                   N(i,j)=0;
               else
                   N(i,j)=1/2;
               end
            elseif ismember([R_rec(i,length(City_num))],Edge)==1
                N(i,end)=1;
            else
                N(i,j)=0;
            end
        end
        [r,c]=find(R_rec_temp(i,:)==0);
        rec_temp=R_rec_temp(i,:);
        rec_temp(c)=[];
        N_all(i,1)=sum(N(i,:));
        SX=S0*(N_all(i)+2)+(length(City_num)-length(rec_temp))*S1;% total left area
        Q_1= Hw*(A*B-SX-S2)*h;
        Q_2= Len(i)*K;              
        f(i,1)= Q_1+Q_2;
        f(i,2)= SX;
    end
    f1_min=min(f(:,1));
    f2_min=min(f(:,2));
    u(:,1)=(f(:,1)-f1_min)/f1_min;
    u(:,2)=(f(:,2)-f2_min)/f2_min;
    for i=1:Ant_num
      U(i)=u(i,1)*w1+u(i,2)*w2; %satisfaction rate for each ants
    end
    %%%%%%%%%%%%%%%%%%%% 
    %update the hermone on whole area
    R0=10; %self-adaption value
    R1=R0*(1-exp(Ite_num-Ite));  % self-adaption function
    %R1=R0*0.5^(Ite);  % self-adaption function
    %R1=R0;
    Delta_Tau = zeros(length(S), length(S)); 
    for i = 1:Ant_num
        [r,c]=find(R_rec_temp(i,:)==0);
        rec_temp=R_rec_temp(i,:);
        rec_temp(c)=[];
        if length(rec_temp)==(length(City_num)-block_num) %if constraints not satsified, do not update the hermone
           for j = 1:(length(City_num)-1)
               Delta_Tau(R_rec(i,j), R_rec(i,j+1)) = Q/(R1*(f(i,1)+f(i,2))+U(i))+Delta_Tau(R_rec(i,j), R_rec(i,j+1)); 
           end
           Rho=0.15;
        else
            Rho=1;
        end
        %Delta_Tau(R_rec(i,length(City_num)), R_rec(i,1)) = Delta_Tau(R_rec(i,length(City_num)), R_rec(i,1))+ Q/(R1*(f(i,1)+f(i,2))+U(i));
    end
    R_rec_temp;
    Tau = Rho.*Tau +(1-Rho)*Delta_Tau; % update, whole area

    %%%%%%%%%%%%%%%%%%%%%%% 
    % Non-dominated sorting and Archiving
    if Ite_num >= 2
        R_rec= [R_rec;best_routes]; % best routes
        f=[f;archived_phermone];
    end
    [FrontNO,MaxFNO] = NDSort(f,inf);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % record the best rout of this iteration
    Index = find(FrontNO==1);
    CrowdDis=CrowdingDistance(f(Index(1:end),:));
    CrowdDis=CrowdDis';
    [Y,IN]=sort(CrowdDis,'descend');
    Index=Index(IN);
    [r,c]=find(Y==inf);
    if length(r)>2
       Y(r,:)=[];
       IN(r,:)=[];
       Index(IN(r,:))=[];
    end
  
    %archived_phermone=zeros(ceil(length(Index)/2),2);
    %best_routes=zeros(ceil(length(Index)/2),length(City_num));
    %archived_phermone=f(Index(1:ceil(length(Index)/2)),:);
    %best_routes=R_rec(Index(1:ceil(length(Index)/2)),:); % best routes(solutions)

    archived_phermone=zeros(length(Index),2); % preallocation for each iteration
    best_routes=zeros(length(Index),length(City_num)); % preallocation for each iteration
    archived_phermone=f(Index(1:end),:);
    best_routes=R_rec(Index(1:end),:); % best routes(solutions)
    % delete the routes invalid
    

    for i=1:size(best_routes,1)
       for j=1:length(City_num)-1
           if best_routes(i,j+1)-best_routes(i,j)==0
              best_routes(i,:)=zeros(1,length(City_num));
              break
           end
       end
    end
   [r0,c0]=find(best_routes==0);
   best_routes(r0,:)=[];
   archived_phermone(r0,:)=[];
   if size(archived_phermone,1)==0
      [L_best(Ite_num), index] = min(f(Index(1:end),2)); % find the best solution of each iteration by minimum of left area
      %[L_best1(Ite_num), index1] = min(f(Index(1:end),1)); % find the best solution of each iteration by minimum of energy consumption
      R_best(Ite_num,:) = R_rec(Index(index),:); % best route
      archived_solution(Ite_num,:)=f(Index(index),:);
      best_routes=R_best(Ite_num,:);
      archived_phermone=archived_solution(Ite_num,:);
   else
      [L_best(Ite_num), index] = min(archived_phermone(1:end,2)); % find the best solution of each iteration by minimum of left area
      %[L_best1(Ite_num), index1] = min(archived_phermone(1:end,1)); % find the best solution of each iteration by minimum of energy consumption
      archived_solution(Ite_num,:)=archived_phermone(index,:); % best results by minimum of left area
      R_best(Ite_num,:) = best_routes(index,:); % best route(solution) by minimum of left area
   end
   % print the outcome of each iteration
%    figure(1);
%    subplot(2,2,1)
%    PlotCosts(f(Index(1:end),:));
%    title(['Pareto Front of Iteration ',num2str(Ite_num)])
%    
%    subplot(2,2,3)
%    PlotCosts(archived_phermone(1:end,:));
%    title(['Pareto valid solution front of Iteration ',num2str(Ite_num)])
% 
%    subplot(2,2,1)
%    PlotCosts(f(Index(1:end),:));
%    title(['Pareto Front of Iteration ',num2str(Ite_num)])
%    
%    subplot(2,2,4)
%    PlotCosts(f); 
%    title(['Function Costs of Iteration ',num2str(Ite_num)])
% 
%    subplot(2,2,2)
%    It_BestR=R_best(Ite_num,:);
%    pcolor(x_grid,y_grid,a);
%    hold on
%    scatter(City(:,1),City(:,2));
%    for ii=2:length(It_BestR)
%        plot([City(It_BestR(ii-1),1),City(It_BestR(ii),1)],[City(It_BestR(ii-1),2),City(It_BestR(ii),2)],'r','LineWidth',2)
%        hold on
%    end
%    title(['Best Route of Iteration ',num2str(Ite_num)]);
%    drawnow;
%    hold off
%    
   if Ite_num==Ite
      f=f(Index(1:ceil(length(Index)/2)),:); % to plot the pareto fronts
   end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clear vivted city table
    R_rec=zeros(Ant_num,length(City_num));
    R_rec_temp=R_rec;
    disp(['Iteration ' num2str(Ite_num)])
    Ite_num=Ite_num+1;
end % ACO ends
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print the outcome

Pos = find(archived_phermone(:,2)==min(archived_phermone(:,2))); % based on best solution of whole algorithm
Shortest_Route=best_routes(Pos(1),:);
%Pos = find(L_best==min(L_best)); % based on best solution of each interation
%Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));

N=length(Shortest_Route);
RS=Shortest_Route;
figure(2)
pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy
hold on
scatter(City(:,1),City(:,2));
hold on
for ii=2:N
    plot([City(RS(ii-1),1),City(RS(ii),1)],[City(RS(ii-1),2),City(RS(ii),2)],'r','LineWidth',2)
    pause(0.2) 
end
title('Best Route Regarding Left Area')
xlabel('Width'); ylabel('Height');

hold off

figure(3) % plot the valid solution from pareto front 1 
PlotCosts(archived_phermone);

figure(4) % plot the pareto front 1
PlotCosts(f);

% judge the quality of the algorithm
vol_0=Hypervolume(archived_phermone);
Vol_Pareto_solution=sum(vol_0)
vol_1=Hypervolume(f);
Vol_Pareto_front=sum(vol_1)

figure(5) % plot the convergence graph with regard to left area
x=linspace(0,Ite,Ite);
y=L_best(:,1);
plot(x,y,'b-o');
xlabel('Iteration'); ylabel('Best Route Regarding Left Area');

% figure(6)   %plot the convergence graph with regard to energy concumption
% x=linspace(0,Ite,Ite);
% y=L_best1(:,1);
% plot(x,y,'g-o');
% xlabel('Iteration'); ylabel('Best Route Regarding energy consumption');

toc