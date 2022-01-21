 %% define  world
%defining 3 obstacles
obstacle = [struct('center',[-4;5],'radius',3,'influence',1)];
%defining workspace
workspace = struct('center',[0;0],'radius',15,'influence',0.5 );

%define
robot = [-9;5];
Goal = [-1;12];

close all
figure
for i = 1:size(robot,2)
color = ['r','b','c','k','m'];
plot_contourwithpath(robot(:,i),obstacle,workspace,Goal,color(:,i));
hold on
end

%plot_totpotential(Goal,obstacle,workspace)
%plot_reppotential(Goal,obstacle,workspace)
%plot_attpotential(Goal,obstacle,workspace)

%plot_potential(Goal,obstacle,workspace);

%% potential function
function [Uatt,Urep] = potential_function(q,Goal,obstacle,workspace)
%attractive function
zeta = 0.2;
eta = 3;
dis_goal = norm(q-Goal);
dist_star = 10;

if dis_goal <= dist_star
    Uatt = 0.5*zeta*dis_goal^2;
else
    Uatt = dist_star*zeta*dis_goal - 0.5*zeta*dist_star^2;
end
%repulsive potential
Urep = 0;
dist_workspace = workspace.radius - norm(q-workspace.center);

if dist_workspace < 0.00005
    Urep = 1/0;
elseif abs(dist_workspace) <= workspace.influence
    Urep = 0.5*eta*((1/dist_workspace)-(1/workspace.influence))^2;
end

for i = 1:size(obstacle,1)
    dist_obs = -obstacle(i).radius + norm(q-obstacle(i).center);
    if dist_obs < 0.00005
       U_rep = 1/0;
       Urep = Urep + U_rep;
    elseif dist_obs <= obstacle(i).influence
        U_rep = 0.5*eta*((1/dist_obs)-(1/obstacle(i).influence))^2;
        Urep = Urep + U_rep;

    end
end
end

%% gradient calculation

function [grad_U] = grad_potential(q,Goal,obstacle)
% gradient attractive function
dis_goal = norm(q-Goal);
dist_star = 10;
zeta = 0.2;
eta = 3;
if dis_goal <= dist_star
    grad_Uatt =zeta*(q-Goal);
else
    grad_Uatt = dist_star*zeta*((q-Goal)/dis_goal);
end
% gradientrepulsive function
grad_Urep = zeros(2,1);

for i = 1:size(obstacle,1)
    dist_obs = -obstacle(i).radius + norm(q-obstacle(i).center);
    if dist_obs <= obstacle(i).influence
       grad_U_rep = -eta*((1/dist_obs)-(1/obstacle(i).influence))*(1/dist_obs)^2*((q-obstacle(i).center)/dist_obs);
       grad_Urep = grad_Urep + grad_U_rep;
    end
    
end

grad_U = grad_Uatt + grad_Urep;
end 

%path calculation
function [path] = path_calc(q,Goal,obstacle)
epsilon  = 0.05;
max_step = 9999;
path = zeros(2,max_step);
path(:,1) = q;

for step = 2:max_step
    pt = path(:,step-1);
    grad = grad_potential(pt,Goal,obstacle);
    path(:,step) = pt - epsilon * grad;
   
    if norm(pt-Goal) < 0.005
        path = path(:,1:step);
        break;
    end
end
end

function [zz_att,zz_rep] = grid(xx,yy,Goal,obstacle,workspace)
zz_att = [];
zz_rep = [];
for i = 1:numel(xx)
    qx = xx(i);
    qy = yy(i);
    q = [qx; qy];
    [Uatt,Urep] = potential_function(q,Goal,obstacle,workspace);
    zz_att(i) = Uatt;
    zz_rep(i) = Urep;
end

zz_att = reshape(zz_att, size(xx));
zz_rep = reshape(zz_rep, size(xx));
end

function plot_totpotential(Goal,obstacle,workspace)
figure
x = linspace(-16, 16);
y =linspace(-16, 16);
[xx, yy] = meshgrid(x, y);
[zz_att,zz_rep] = grid(xx,yy,Goal,obstacle,workspace); 
threshold = max(zz_att,[],'all');
zz = zz_att + zz_rep;
zz(zz > threshold) = threshold;
surf(xx, yy, zz)
end

function plot_attpotential(Goal,obstacle,workspace)
figure
x = linspace(-16, 16);
y =linspace(-16, 16);
[xx, yy] = meshgrid(x, y);
[zz_att] = grid(xx,yy,Goal,obstacle,workspace); 
surf(xx, yy, zz_att)
end

function plot_reppotential(Goal,obstacle,workspace)
figure
x = linspace(-16, 16);
y =linspace(-16, 16);
[xx, yy] = meshgrid(x, y);
[zz_att,zz_rep] = grid(xx,yy,Goal,obstacle,workspace); 
threshold = max(zz_att,[],'all');
zz = zz_rep;
zz(zz > threshold) = threshold; 
surf(xx, yy, zz)
end

function [xx,yy] = pot_meshgrid(workspace) 

x_min = workspace.center(1) - workspace.radius;
x_max = workspace.center(1) + workspace.radius;
y_min = workspace.center(2) - workspace.radius;
y_max = workspace.center(2) + workspace.radius;
x = linspace(x_min - 0.1,x_max + 0.1);
y = linspace(y_min - 0.1, y_max + 0.1);

[X,Y] = meshgrid(x, y);
xx = X;
yy = Y;
end

function plot_contourwithpath(q,obstacle,workspace,Goal,color)
[xx, yy] = pot_meshgrid(workspace);
[zz_att,zz_rep] = grid(xx,yy,Goal,obstacle,workspace); 
threshold = max(zz_att,[],'all');
zz = zz_att + zz_rep;
zz(zz > threshold) = threshold;
contour(xx, yy,zz,50)
colormap(winter(10));

hold on

path = path_calc(q,Goal,obstacle);
plot(path(1,:), path(2,:),color,'LineWidth',2); %plots path
plot(Goal(1),Goal(2),"rx")
end

