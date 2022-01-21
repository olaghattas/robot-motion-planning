%% Define Variables

row_num =25;
col_num = 25;

obst1 = [5,5,6,7,8,9,9;5,6,5,5,5,5,6];
obst2 = [13,13,14,14,15,16,17,18,19,20,21,15,16,17,18,19,20,21;5,6,5,6,5,5,5,5,5,5,5,6,6,6,6,6,6,6];
obst3 = [3,4,5,6,7,8,16,17,18,19,20,21,22,23,3,4,5,6,7,8,16,17,18,19,20,21,22,23;16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17];
obst4 = [12,12;16,17];
out_edge1 = [ones(1,col_num); 1:row_num];
out_edge2 = [2:col_num-1;ones(1,row_num-2)];
out_edge3 = [2:col_num-1;ones(1,row_num-2)*row_num];
out_edge4 = [ones(1,col_num)*col_num;1:row_num];

obsts = [obst1 obst2 obst3 obst4  out_edge1 out_edge2 out_edge3 out_edge4];
obst_row = [struct('obst',obst1);struct('obst',obst2);struct('obst',obst3);struct('obst',obst4);struct('obst',out_edge1);struct('obst',out_edge2);struct('obst',out_edge3);struct('obst',out_edge4)];
brushfireGrid = zeros(row_num,col_num);

close all;
[distanceMap, Gener_vd] = dist_map(brushfireGrid,obsts,obst_row,row_num,col_num);
brushfireGrid = map_grid(distanceMap,col_num,row_num);
figure
plot_gr(brushfireGrid,Gener_vd,col_num,obsts,row_num);

hold on
% path = get_path_to_goal(brushfireGrid, Gener_vd, [6;18], [17;18],row_num,col_num);
% plot(path(2,:) - 0.5,path(1,:)-0.5,'k-','linewidth',2);
% hold on
% path = get_path_to_goal(brushfireGrid, Gener_vd, [20;7], [15,23],row_num,col_num);
% plot(path(2,:) - 0.5,path(1,:)-0.5,'r-','linewidth',2);
% hold on
path = get_path_to_goal(brushfireGrid, Gener_vd, [10;9], [20;2],row_num,col_num);
plot(path(2,:) - 0.5,path(1,:)-0.5,'g-','linewidth',2);
hold on
%path = get_path_to_goal(brushfireGrid, Gener_vd, [5;2], [5;24],row_num,col_num);
%plot(path(2,:) - 0.5,path(1,:)-0.5,'c-','linewidth',2);
%% Plot Grid
function plot_gr(brushfire_grid, Gener_vd,col_num,obsts,row_num)
plot([0:row_num; 0:row_num], [0:col_num; 0:col_num], 'k')
hold on
hax = gca;
hax.XTick = 0:row_num;
hax.YTick = 0:col_num;
hax.XTickLabel = [];
hax.YTickLabel = [];
grid
axis square
% fill grid with values
for iRow = 1:row_num
    for iCol = 1:col_num
        value  = brushfire_grid(iRow ,iCol);
        value = value + 1;
        hold on
        text((iCol - 1)+0.5,(iRow -1)+0.5,num2str(value),'FontWeight','bold');
        if value == 1
            fill([0 1 1 0]+(iCol - 1), [0 0 1 1]+(iRow-1), 'b')
        end
    end
end


for i = 1:size(Gener_vd,2)
    point = Gener_vd(:,i);
    [manhattan_neighbors, diagonal_neighbors] = man_daig_neig(point(1),point(2),row_num,col_num);
    allneighbors = horzcat(manhattan_neighbors, diagonal_neighbors);
    for j=1:size(allneighbors,2)
        neighbor = allneighbors(:,j);
        bool = is_cell_on_voronoi_boundary(Gener_vd, neighbor);
        if bool
            plot([point(2) neighbor(2)] - 0.5, [point(1) neighbor(1)] - 0.5,'m.-','linewidth',2);
            hold on
        end
    end
end
end

function [bool] = is_cell_on_voronoi_boundary(Gener_vd, point)
for i = 1:size(Gener_vd,2)
    if isequal(point, Gener_vd(:,i))
        bool = true;
        return;
    end
end
bool = false;
end

function [D, Gener_vd] = dist_map(brushfireGrid,obsts,obst_row,row_num,col_num)
% fill all obsts with 1
for i = 1:size(obsts,2)
    obst = obsts(:,i);
    brushfireGrid(obst(1),obst(2)) = 1;
end

M = grid_map(brushfireGrid,row_num,col_num);
D = containers.Map();
obst = containers.Map();
open = priority_prepare();
Gener_vd = [];

% Initialize data structures
for iCol = 1:col_num
    for iRow = 1:row_num
        key = iRow + "-" + iCol;
        if M(key) == 1
            D(key) = 0;
            obst(key) = key;
            [open,~] = priority_insert(open, key, 0);
        else
            D(key) = 1/0;
        end
    end
end
% Update Grid
while size(open,1) > 0
    [open,key,~] = priority_minExtract(open);
    coordinates = split(key,'-');
    s_x = str2double(coordinates(1));
    s_y = str2double(coordinates(2));
    obst_s_value = obst(key);
    obst_s = split(obst_s_value,'-');
    obst_s = str2double(obst_s);
    % lower
    [manhattan_neighbors, diagonal_neighbors] = man_daig_neig(s_x,s_y,row_num,col_num);
    allneighbors = horzcat(manhattan_neighbors, diagonal_neighbors);
    for i = 1:size(allneighbors,2)
        neighbor = allneighbors(:,i);
        neighbor_key = neighbor(1) + "-" + neighbor(2);
        d =  int8(norm(obst_s - neighbor));
        if d < D(neighbor_key)
            D(neighbor_key) = d;
            obst(neighbor_key) = obst(key);
            [open, ~] = priority_insert(open, neighbor_key,d);
        elseif d == D(neighbor_key)
            obst_n_value = obst(neighbor_key);
            obst_n = split(obst_n_value,'-');
            obst_n = str2double(obst_n);
            if ~Same_obst(obst_n, obst_s,obst_row)
                Gener_vd = add_cell_to_voronoi_boundary(Gener_vd, neighbor);
            end
        end
    end
end
end

function [Gener_vd] = add_cell_to_voronoi_boundary(Gener_vd, cell)
for i = 1:size(Gener_vd,2)
    if isequal(cell, Gener_vd(:,i))
        return;
    end
end
Gener_vd = [Gener_vd cell];
end

function [gridMap] = grid_map(grid,row_num,col_num)
gridMap = containers.Map();
for iCol = 1:row_num
    for iRow = 1:col_num
        key = iCol + "-" + iRow;
        gridMap(key) = grid(iCol,iRow);
    end
end
end

function [grid] = map_grid(map,row_num,col_num)
grid = zeros(row_num,col_num);
for iCol = 1:row_num
    for iRow = 1:col_num
        key = iCol + "-" + iRow;
        grid(iCol,iRow) = map(key);
    end
end
end

%% Get manhattan and diagonal neighbors of a grid cell
function [manhattan_neighbors, diagonal_neighbors] = man_daig_neig(x,y,row_num,col_num)
manhattan_neighbors = [];
diagonal_neighbors = [];
max_cols = col_num;
max_rows = row_num;
% Manh
if (x + 1 <= max_cols)
    manhattan_neighbors = horzcat(manhattan_neighbors,[x + 1; y]);
end
if (x - 1 >= 1)
    manhattan_neighbors = horzcat(manhattan_neighbors,[x - 1; y]);
end
if (y + 1 <= max_rows)
    manhattan_neighbors = horzcat(manhattan_neighbors,[x; y + 1]);
end
if (y - 1 >= 1)
    manhattan_neighbors = horzcat(manhattan_neighbors,[x; y - 1]);
end
% Diag
if (x + 1 <= max_cols && y + 1 <= max_rows)
    diagonal_neighbors = horzcat(diagonal_neighbors,[x + 1; y + 1]);
end
if (x + 1 <= max_cols && y - 1 >= 1)
    diagonal_neighbors = horzcat(diagonal_neighbors,[x + 1; y - 1]);
end
if (x - 1 >= 1 && y + 1 <= max_rows)
    diagonal_neighbors = horzcat(diagonal_neighbors,[x - 1; y + 1]);
end
if (x - 1 >= 1 && y - 1 >= 1)
    diagonal_neighbors = horzcat(diagonal_neighbors,[x - 1; y - 1]);
end
end

function [bool] = Same_obst(point1, point2,obst_row)

if isequal(point1, point2)
    bool = true;
    return;
end
bool = false;
for i=1:size(obst_row,1)
    
    obst = obst_row(i,1).obst;
    bool1 = false;
    bool2 = false;
    for j = 1:size(obst,2)
        if isequal(obst(:,j), point1)
            bool1 = true;
        end
        if isequal(obst(:,j), point2)
            bool2 = true;
        end
    end
    if bool1 && bool2
        bool = true;
        return
    end
end
end
%% from robot motion planning pQueue
function [pQueue, cellUpdated]=priority_insert(pQueue,key,distance)
for i=1:size(pQueue,1)
    if pQueue(i).key == key
        cellUpdated = true;
        pQueue(i).distance = distance;
        return;
    end
end
cellUpdated = false;
pQueue=[pQueue;struct('key',key,'distance',distance)];
end
function [pQueue,key,dist]=priority_minExtract(pQueue)
elemNum=numel(pQueue);%get number of records
if elemNum==0%if the given queue is empty, return an empty queue.
    key=[];%set key as empty.
    dist=[];%set cost as empty.
    return;
end
[~,idxMin]=min([pQueue.distance]);
idxMin=idxMin(1);
key=pQueue(idxMin).key;
dist=pQueue(idxMin).distance;
pQueue(idxMin)=[];
if isempty(pQueue)
    pQueue=priority_prepare();
end
end
function [pQueue]=priority_prepare()
pQueue=repmat(struct('key',[],'distance',[]),0,1);
end

%% path planning
function [path] = get_path_to_goal(grid,Gener_vd,start,goal,row_num,col_num)
up_path = move_up(grid,Gener_vd,start,row_num,col_num);
down_path = flip(move_up(grid,Gener_vd,goal,row_num,col_num),2);
Gener_vd_path = get_voronoi_path(grid,up_path(:,end), down_path(:,1),row_num,col_num);

up_path(:,end) = [];
down_path(:,1) = [];
path = [up_path Gener_vd_path down_path];
end

function [path] = get_voronoi_path(grid,start, goal,row_num,col_num)
path = [];
queue = priority_prepare();
queue =  priority_insert(queue,start,0);

came_from = containers.Map();
cost_so_far = containers.Map();

start_key = start(1)+"-"+start(2);

came_from(start_key) = NaN;
cost_so_far(start_key) = 0;

while size(queue,1) > 0
    [queue,cell,~] = priority_minExtract(queue);
    x = cell(1);
    y = cell(2);
    cell_key = x+"-"+y;
    
    if isequal(cell,goal)
        path = [];
        path(:,1) = goal;
        current = goal;
        current_key = current(1)+"-"+current(2);
        while isKey(came_from,current_key)
            if isnan(came_from(current_key))
                break;
            end
            current = came_from(current_key);
            current_key = current(1)+"-"+current(2);
            path = [path current];
        end
        path = [path current];
        path = flip(path,2);
        
        
        return
    end
    
    [manhattan_neighbors, diagonal_neighbors] = man_daig_neig(x,y,row_num,col_num);
    allneighbors = [manhattan_neighbors diagonal_neighbors];
    
    for i = 1:size(allneighbors,2)
        neighbor = allneighbors(:,i);
        v = grid(neighbor(1),neighbor(2));
        if v == 1
            continue
        end
        neighbor_key = neighbor(1)+"-"+neighbor(2);
        candidate_cost = cost_so_far(cell_key) + 1;
        
        expand = false;
        if isKey(cost_so_far,neighbor_key)
            if candidate_cost < cost_so_far(neighbor_key)
                expand = true;
            end
        else
            expand = true;
        end
        
        if expand
            came_from(neighbor_key) = cell;
            cost_so_far(neighbor_key) = candidate_cost;
            [queue] = priority_insert(queue,neighbor,candidate_cost);
        end
    end
end
end

function [path] = move_up(grid,Gener_vd,start,row_num,col_num)
path = [];
path(:,1) = start;

while true
    cell = path(:,end);
    if is_cell_on_voronoi_boundary(Gener_vd, cell)
        break;
    end
    
    [manhattan_neighbors, diagonal_neighbors] =  man_daig_neig(cell(1), cell(2),row_num,col_num);
    max_value = 0;
    max_neighbor = NaN;
    allneighbors = [manhattan_neighbors diagonal_neighbors];
    
    for i = 1:size(allneighbors,2)
        neighbor = allneighbors(:,i);
        if grid(neighbor(1), neighbor(2)) == 1
            continue
        elseif grid(neighbor(1), neighbor(2)) > max_value
            max_value = grid(neighbor(1), neighbor(2));
            max_neighbor = neighbor;
        end
    end
    
    if isnan(max_neighbor)
        break;
    else
        path =[path max_neighbor];
    end
end
end

