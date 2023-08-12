clc
clear
hermone_case_1_Ln1=load('viableS_case_1_In3.mat')
hermone_case_2_Ln1=load('viableS_case_2_In3.mat')
hermone_case_3_Ln1=load('viableS_case_3_In3.mat')
hermone_case_4_Ln1=load('viableS_case_4_In3.mat')
% figure(3) % plot the valid solution from pareto front 1 
% PlotCosts(hermone_case_1_Ln1.archived_phermone);
% hold on
% PlotCosts1(hermone_case_2_Ln1.archived_phermone);
% PlotCosts2(hermone_case_3_Ln1.archived_phermone);
% hold off

figure(4) % plot the pareto front 1
subplot(1,2,1)
PlotCosts(hermone_case_1_Ln1.f);
hold on
PlotCosts1(hermone_case_2_Ln1.f);
PlotCosts2(hermone_case_3_Ln1.f);
PlotCosts3(hermone_case_4_Ln1.f);
% judge the quality of the algorithm
subplot(1,2,2)% plot the convergence graph with regard to left area
plot(hermone_case_1_Ln1.x,hermone_case_1_Ln1.y,'r-o');
hold on
plot(hermone_case_2_Ln1.x,hermone_case_2_Ln1.y,'g-o');
plot(hermone_case_3_Ln1.x,hermone_case_3_Ln1.y,'b-o');
plot(hermone_case_4_Ln1.x,hermone_case_4_Ln1.y,'y-o');
xlabel('Iteration'); ylabel('Best Route Regarding Left Area');
title('Convergence Regarding Left Area')


x_grid=hermone_case_1_Ln1.x_grid;
y_grid=hermone_case_1_Ln1.y_grid;
N=hermone_case_1_Ln1.N;
a=hermone_case_1_Ln1.a;
City=hermone_case_1_Ln1.City;
RS1=hermone_case_1_Ln1.RS;
RS2=hermone_case_2_Ln1.RS;
RS3=hermone_case_3_Ln1.RS;
RS4=hermone_case_4_Ln1.RS;
figure(2)
subplot(1,4,1)
pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy
hold on
scatter(City(:,1),City(:,2));
hold on
for ii=2:N
    plot([City(RS1(ii-1),1),City(RS1(ii),1)],[City(RS1(ii-1),2),City(RS1(ii),2)],'r','LineWidth',2)
end
title('Best Route Regarding Left Area')
xlabel('Width'); ylabel('Height');
hold off
subplot(1,4,2)
pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy
hold on
scatter(City(:,1),City(:,2));
hold on
for ii=2:N
    plot([City(RS2(ii-1),1),City(RS2(ii),1)],[City(RS2(ii-1),2),City(RS2(ii),2)],'r','LineWidth',2)
end
title('Best Route Regarding Left Area')
xlabel('Width'); ylabel('Height');
hold off
subplot(1,4,3)
pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy
hold on
scatter(City(:,1),City(:,2));
hold on
for ii=2:N
    plot([City(RS3(ii-1),1),City(RS3(ii),1)],[City(RS3(ii-1),2),City(RS3(ii),2)],'r','LineWidth',2)
end
title('Best Route Regarding Left Area')
xlabel('Width'); ylabel('Height');
hold off
subplot(1,4,4)
pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy
hold on
scatter(City(:,1),City(:,2));
hold on
for ii=2:N
    plot([City(RS4(ii-1),1),City(RS4(ii),1)],[City(RS4(ii-1),2),City(RS4(ii),2)],'r','LineWidth',2)
end
title('Best Route Regarding Left Area')
xlabel('Width'); ylabel('Height');
hold off