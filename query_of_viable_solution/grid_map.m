clear all
close all
clc
% create the gripmap
a = ones(40);
a(3,3:7)=0;
a(4,3:7)=0;
a(5,3:7)=0;
a(16:17,13:17)=0;
a(10:17,13:14)=0;
a(16,6:7)=0
a(17,6:7)=0
a(15,7)=0
a(15,15)=0
a(10:15,20)=0
a(14,15)=0;
a(21:22,25:36)=0
a(10:14,21)=0
a(28,18:24)=0;
a(27,18:26)=0;
a(26,17:30)=0;
a(36:37,15:30)=0;
a(35,13:25)=0;
a(6:11,30:35)=0;
a(7:12,36)=0;
a(12,31:35)=0;

b=a;
b(end+1,end+1)=0;
a=b
colormap([0 0 0;1 1 1]);  % use colormap to create the wanted color
X_grid=[1:size(a,2)].*ones(size(a,1),1) 
Y_grid=X_grid'

%create the decision variable
Nx=size(a,2)-1;
Ny=size(a,1)-1;
n=Nx*Ny/4
n1=Nx/2
n2=Ny/2
XY=zeros(n2,n1);
coordinatesXY=zeros(n2*n1,2)
S=zeros(1,n2*n1); % coordinates sequence
ratio=10 % cutting width unit/head diameter

x_grid=(X_grid(1,:)-1)/ratio; 
y_grid=(Y_grid(:,1)-1)/ratio;  
x=x_grid(2:2:[length(x_grid)-1]) % x coordinates when ratio=2 x=[0.5:n1-0.5]; 
y=y_grid(2:2:[length(y_grid)-1]) % y coordinates when ratio=2 y=[0.5:n2-0.5]; 

pcolor(x_grid,y_grid,a);
set(gca,'XTick',1:size(a,2),'YTick',1:size(a,1));  % axis setting
axis image xy; 
hold on

% self version coordinates index to variable XY index S 
k=1
for i=1:n2
     for j=1:n1
         XY(i,j)=k;
         S(k)=k;
         coordinatesXY(k,:)=[x(j),y(i)];
         plot(x(j),y(i),'ro')
         hold on
         k=k+1;
     end
end

% aco function
%[R_best,L_best,L_ave,Shortest_Route,Shortest_Length]=ACATSP(coordinatesXY,Ite,Ant_num,Alpha,Beta,Rho,Q)
