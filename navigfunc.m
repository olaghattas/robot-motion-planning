%% define  world
%defining 3 obstacles
obstacle = [struct('center',[-4;5],'radius',3,'influence',1); struct('center',[-3;-7],'radius',2,'influence',1.2); struct('center',[4;4],'radius',4,'influence',1)];
%defining workspace
workspace = struct('center',[0;0],'radius',15,'influence',0.5 );

%define
robot = [-9,9,-2,-4,11;5,0,-1,-10,-5];
Goal = [-1;12];

% positive parameter K
K= 9;
close all;
%plot_navcontour(true,K,robot,Goal,obstacle,workspace);
plot_navsurface(K,Goal,obstacle,workspace)
figure
for i = 1:size(robot,2)
color = ['r','b','c','k','m'];
plot_navcontour(true,K,robot(:,i),Goal,obstacle,workspace,color(:,i));
hold on
end

%path calculation
function [path] = path_calc(robot,K,Goal,obstacle,workspace)
scale_magnitude = 0.1;
max_step = 10000;
epsilon = 1; % learning rate
starting_distance = dist2goal(robot,Goal);
path = zeros(2,max_step);
path(:,1) = robot;

for step = 2:max_step
    q = path(:,step-1);
    grad = get_grad(q,K,Goal,obstacle,workspace);
    q_goal = min(norm(q - Goal), starting_distance); %distance to goal
    if  ~isnan (scale_magnitude)
        eta = (q_goal/ starting_distance)^2;
        grad = eta * scale_magnitude * grad / norm(grad);
    end
    path(:,step) = q - epsilon * grad;
   
    
    if q_goal < 0.005
        path = path(:,1:step);
        break;
    end
end
end

%% Plotting

function plot_navcontour(withPath,K,robot,Goal,obstacle,workspace,color)
x = linspace(-20, 20);
y =linspace(-20, 20);
[xx, yy] = meshgrid(x, y);
zz_nav = grid(xx,yy,K,Goal,obstacle,workspace);

contour(xx, yy, zz_nav,'FaceColor','b')
hold on

for i=1:numel(obstacle)
    plot(obstacle(i).center(1,1), obstacle(i).center(2,1),'b');
    hold on
end

plot(Goal(1),Goal(2),"rx");
hold on

if withPath
    plot_path(K,robot,Goal,obstacle,workspace,color)
end
end

function plot_path(K,robot,Goal,obstacle,workspace,color) 
% robot is the start point
if size(robot,2) > 1
    for i = 1:size(robot,2)
        start = robot(:,i);
        path = path_calc(start,K,Goal,obstacle,workspace);
        plot(path(1,:),path(2,:),'-');
        hold on
    end
else
    path = path_calc(robot,K,Goal,obstacle,workspace);
    plot(path(1,:),path(2,:),color,'LineWidth',3);
end
end

function plot_navsurface(K,Goal,obstacle,workspace)
figure
x = linspace(-16, 16);
y =linspace(-16, 16);
[xx, yy] = meshgrid(x, y);
zz_nav = grid(xx,yy,K,Goal,obstacle,workspace);

surf(xx, yy, zz_nav)
end

function [zz_nav] = grid(xx,yy,K,Goal,obstacle,workspace)
zz_nav = [];
for i = 1:numel(xx)
    qx = xx(i);
    qy = yy(i);
    q = [qx; qy];
    [~,~,~,~,phi] = phi_values(q,K,Goal,obstacle,workspace);
    zz_nav(i) = phi;
end
zz_nav =reshape(zz_nav, size(xx));
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

function plot_contourwithpath(robot,obstacle,workspace,K,Goal,color)

[xx, yy] = pot_meshgrid(workspace);
zz_nav = grid(xx,yy,K,Goal,obstacle,workspace);
contour(xx, yy, zz_nav)
mycolors = [1 0 0; 1 1 0; 0 0 1];
colormap(mycolors);

hold on
path = path_calc(robot,K,Goal,obstacle,workspace);
plot(path(1,:), path(2, :),color,'LineWidth',2) %plots path

end


%% functions
function [beta_0,beta_q,gamma_q,alpha_q,phi] = phi_values(robot,K,Goal,obstacle,workspace)

dist_obs = 1;
for i_row = 1:size(obstacle,1)
    distanceobts= distance_Obs(robot,obstacle(i_row).center,obstacle(i_row).radius);
    dist_obs = distanceobts * dist_obs ;
end
beta_0 = -norm(robot-workspace.center)^2+ workspace.radius^2;
beta_q = dist_obs * beta_0; 


[dist_goal] = dist2goal(robot,Goal);

gamma_q = dist_goal.^(2*K);
alpha_q = gamma_q/beta_q;

%navigation function
if alpha_q < 0
    phi = 1;
else
    phi = (alpha_q / (1 + alpha_q)).^(1/K);
end
end

function [grad_phi]=get_grad(robot,K,Goal,obstacle,workspace)
%call values needed
[beta_0,beta_q,gamma_q,alpha_q] = phi_values(robot,K,Goal,obstacle,workspace);
%grad gamma
[dist_goal] = dist2goal(robot,Goal);
grad_gamma = 2*K* dist_goal^(2*K-1) * (robot-Goal)/ dist_goal;
%display(grad_gamma)

%grad beta
grad_beta_0 = -2 * (robot-workspace.center);

betas = ones(1,size(obstacle,1)+1);
betas(1,1) = beta_0;

grad_betas = ones(2,size(obstacle,1)+1);
grad_betas(:,1) = grad_beta_0;


for i_row = 1:size(obstacle,1)
    betas(1,i_row) = distance_Obs(robot,obstacle(i_row).center,obstacle(i_row).radius);
    grad_betas(:,i_row) = 2 * (robot - obstacle(i_row).center);
end

grad_beta = zeros(2,1);
n = numel(betas);
for i =1:n
    product = grad_betas(i);
    for j = 1:n
        if j ~= i
            product = product *betas(j);
        end
    end
    grad_beta = grad_beta + product;
end
%grad alpha
grad_alpha = (grad_gamma * beta_q - gamma_q * grad_beta) / beta_q^2;
%grad phi
grad_phi = (1 / K) * (alpha_q / (1 + alpha_q))^((1 - K)/ K) * (1 / (1 + alpha_q)^2)* grad_alpha;
end

%distance calc
function [dist_goal] = dist2goal(robot,Goal)
dist_goal = norm(robot-Goal);
end
function [dist_obs] = distance_Obs(robot,center,radius)
dist_obs = norm(robot-center)^2 - radius^2;
end