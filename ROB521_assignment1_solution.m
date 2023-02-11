% ======
% ROB521_assignment1.m
% ======
%
% This assignment will introduce you to the idea of motion planning for  
% holonomic robots that can move in any direction and change direction of 
% motion instantaneously.  Although unrealistic, it can work quite well for
% complex large scale planning.  You will generate mazes to plan through 
% and employ the PRM algorithm presented in lecture as well as any 
% variations you can invent in the later sections.
% 
% There are three questions to complete (5 marks each):
%
%    Question 1: implement the PRM algorithm to construct a graph
%    connecting start to finish nodes.
%    Question 2: find the shortest path over the graph by implementing the
%    Dijkstra's or A* algorithm.
%    Question 3: identify sampling, connection or collision checking 
%    strategies that can reduce runtime for mazes.
%
% Fill in the required sections of this script with your code, run it to
% generate the requested plots, then paste the plots into a short report
% that includes a few comments about what you've observed.  Append your
% version of this script to the report.  Hand in the report as a PDF file.
%
% requires: basic Matlab, 
%
% S L Waslander, January 2022
%
clear; close all; clc;

% set random seed for repeatability if desired
rng('default');

% ==========================
% Maze Generation
% ==========================
%
% The maze function returns a map object with all of the edges in the maze.
% Each row of the map structure draws a single line of the maze.  The
% function returns the lines with coordinates [x1 y1 x2 y2].
% Bottom left corner of maze is [0.5 0.5], 
% Top right corner is [col+0.5 row+0.5]
%

row = 5; % Maze rows
col = 7; % Maze columns
map = maze(row,col); % Creates the maze
start = [0.5, 1.0]; % Start at the bottom left
finish = [col+0.5, row]; % Finish at the top right

h = figure(1);clf; hold on;
plot(start(1), start(2),'go')
plot(finish(1), finish(2),'rx')
show_maze(map,row,col,h); % Draws the maze
drawnow;

% ======================================================
% Question 1: construct a PRM connecting start and finish
% ======================================================
%
% Using 500 samples, construct a PRM graph whose milestones stay at least 
% 0.1 units away from all walls, using the MinDist2Edges function provided for 
% collision detection.  Use a nearest neighbour connection strategy and the 
% CheckCollision function provided for collision checking, and find an 
% appropriate number of connections to ensure a connection from  start to 
% finish with high probability.


% variables to store PRM components
nS = 500;  % number of samples to try for milestone creation
milestones = [start; finish];  % each row is a point [x y] in feasible space
edges = [];  % each row is should be an edge of the form [x1 y1 x2 y2]

disp("Time to create PRM graph")
tic;
% ------insert your PRM generation code here-------

% uniform random sampling
r1 = rand(500,1);
r2 = rand(500,1);

mins = min(map);
maxs = max(map);
xmin = min([mins(1) mins(3)]);
ymin = min([mins(2) mins(4)]);
xmax = max([maxs(1) maxs(3)]);
ymax = max([maxs(2) maxs(4)]);

pt_x = xmin + r1 * col;
pt_y = ymin + r2 * row;
samples = [pt_x pt_y];

% collision checking for the samples
% add to milestones list if collision-free
for i = 1:nS
    sampled_pt = samples(i,:);
    min_d_pt = MinDist2Edges(sampled_pt, map);
    if min_d_pt > 0.1
        milestones = [milestones; sampled_pt];
    end
end

% for each milestone, find its cloest 8 neighbours
% connect with an edge if the edge is collision-free
k_neighbours = 8;
n_milestones = length(milestones);

dist_milestones = zeros(n_milestones);
adj_mat = zeros(n_milestones);
for i = 1:n_milestones
    i_milestone = milestones(i, :);
    for j = i:n_milestones
        j_milestone = milestones(j, :);
        ij_dist = norm(i_milestone - j_milestone);
        dist_milestones(i, j) = ij_dist;
        dist_milestones(j, i) = ij_dist;
    end

    [dists, idxs] = mink(dist_milestones(i, :), k_neighbours + 1);
    for j = 2:k_neighbours + 1
        edge = [i_milestone milestones(idxs(j), :)];
        edge_idx = [i idxs(j)];
        [inCollision, tmp] = CheckCollision(i_milestone, milestones(idxs(j), :), map);
        if ~inCollision
            edges = [edges; edge];
            adj_mat(i, idxs(j)) = 1;
            adj_mat(idxs(j), i) = 1;
        end
    end
end

% build an adjacency matrix for the graph for running part 2
% each entry is the length of the edge or 0 if no edge exists
dist_milestones = dist_milestones .* adj_mat

% ------end of your PRM generation code -------
toc;

figure(1);
plot(milestones(:,1),milestones(:,2),'m.');
if (~isempty(edges))
    line(edges(:,1:2:3)', edges(:,2:2:4)','Color','magenta') % line uses [x1 x2 y1 y2]
end
str = sprintf('Q1 - %d X %d Maze PRM', row, col);
title(str);
drawnow;

print -dpng assignment1_q1.png


% =================================================================
% Question 2: Find the shortest path over the PRM graph
% =================================================================
%
% Using an optimal graph search method (Dijkstra's or A*) , find the 
% shortest path across the graph generated.  Please code your own 
% implementation instead of using any built in functions.

disp('Time to find shortest path');
tic;

% Variable to store shortest path
spath = []; % shortest path, stored as a milestone row index sequence


% ------insert your shortest path finding algorithm here-------

% Dijkstra's algorithm 
dist_nodes = Inf * ones(n_milestones,1);
parent = -1 * ones(n_milestones, 1);
dist_nodes(1) = 0;

queue = 1:n_milestones;
success = false;

while length(queue) > 0
    
    % find min node
    min_val = Inf;
    min_idx = -1;
    min_q_idx = -1;
    for j = 1:length(queue)
        node_j = queue(j);
        if dist_nodes(node_j) < min_val
            min_val = dist_nodes(node_j);
            min_q_idx = j;
            min_idx = node_j;
        end
    end
    if min_q_idx == -1
        break
    end
    queue = [queue(1:min_q_idx-1) queue(min_q_idx+1:end)];
    
    % update parent
    for j = 1:n_milestones
       if dist_milestones(min_idx, j) > 0
           if any(queue == j)
               if dist_nodes(min_idx) + dist_milestones(min_idx, j) < dist_nodes(j)
                   dist_nodes(j) = dist_nodes(min_idx) + dist_milestones(min_idx, j);
                   parent(j) = min_idx;
               end
           end
       end
    end
    
    if all(milestones(min_idx, :) == finish)
        success = true;
        break
    end
    
    
end

% store path
if success
    path_reverse = [2];
    i = 2;
    while true
       if all(milestones(i, :) == start)
           break
       end
       i = parent(i);
       path_reverse = [path_reverse i];
    end
    spath = flip(path_reverse);
end

    
    
% ------end of shortest path finding algorithm------- 
toc;    

% plot the shortest path
figure(1);
for i=1:length(spath)-1
    plot(milestones(spath(i:i+1),1),milestones(spath(i:i+1),2), 'go-', 'LineWidth',3);
end
str = sprintf('Q2 - %d X %d Maze Shortest Path', row, col);
title(str);
drawnow;

print -dpng assingment1_q2.png


%% ================================================================
% Question 3: find a faster way
% ================================================================
%
% Modify your milestone generation, edge connection, collision detection 
% and/or shortest path methods to reduce runtime.  What is the largest maze 
% for which you can find a shortest path from start to goal in under 20 
% seconds on your computer? (Anything larger than 40x40 will suffice for 
% full marks)

clear; close all; clc;
row = 80;
col = 80;
map = maze(row,col);
start = [0.5, 1.0];
finish = [col+0.5, row];
milestones = [start; finish];  % each row is a point [x y] in feasible space
edges = [];  % each row is should be an edge of the form [x1 y1 x2 y2]

h = figure(2);clf; hold on;
plot(start(1), start(2),'go')
plot(finish(1), finish(2),'rx')
show_maze(map,row,col,h); % Draws the maze
drawnow;

fprintf("Attempting large %d X %d maze... \n", row, col);
tic;        
% ------insert your optimized algorithm here------

% workspace-guided sampling
% sample points at the middle of every adjacent wall pairs
xset = unique(sort([map(:, 1) map(:, 3)]));
yset = unique(sort([map(:, 2) map(:, 4)]));
pt_x = (xset(1:end-1) + xset(2:end)) / 2;
pt_y = (yset(1:end-1) + yset(2:end)) / 2;

spath = [];
for i = 1:length(pt_x)
    for j = 1:length(pt_y)
        milestones = [milestones; [pt_x(i) pt_y(j)]];
    end
end

% locate closest 4 neighbours
% connect with edge if collision-free
k_neighbours = 4;
n_milestones = length(milestones);


dist_milestones = zeros(n_milestones);
adj_mat = zeros(n_milestones);
for i = 1:n_milestones
    i_milestone = milestones(i, :);
    for j = i:n_milestones
        j_milestone = milestones(j, :);
        ij_dist = norm(i_milestone - j_milestone);
        dist_milestones(i, j) = ij_dist;
        dist_milestones(j, i) = ij_dist;
    end

    [dists, idxs] = mink(dist_milestones(i, :), k_neighbours + 1);
    for j = 2:k_neighbours + 1
        edge = [i_milestone milestones(idxs(j), :)];
        edge_idx = [i idxs(j)];
        [inCollision, tmp] = CheckCollision(i_milestone, milestones(idxs(j), :), map);
        if ~inCollision
            edges = [edges; edge];
            adj_mat(i, idxs(j)) = 1;
            adj_mat(idxs(j), i) = 1;
        end
    end
end
dist_milestones = dist_milestones .* adj_mat;

% Dijkstra's algorithm 
dist_nodes = Inf * ones(n_milestones,1);
parent = -1 * ones(n_milestones, 1);
dist_nodes(1) = 0;

queue = 1:n_milestones;
success = false;

while length(queue) > 0
    
    % find min node
    min_val = Inf;
    min_idx = -1;
    min_q_idx = -1;
    for j = 1:length(queue)
        node_j = queue(j);
        if dist_nodes(node_j) < min_val
            min_val = dist_nodes(node_j);
            min_q_idx = j;
            min_idx = node_j;
        end
    end
    if min_q_idx == -1
        break
    end
    queue = [queue(1:min_q_idx-1) queue(min_q_idx+1:end)];
    
    % update parent
    for j = 1:n_milestones
       if dist_milestones(min_idx, j) > 0
           if any(queue == j)
               if dist_nodes(min_idx) + dist_milestones(min_idx, j) < dist_nodes(j)
                   dist_nodes(j) = dist_nodes(min_idx) + dist_milestones(min_idx, j);
                   parent(j) = min_idx;
               end
           end
       end
    end
    
    if all(milestones(min_idx, :) == finish)
        success = true;
        break
    end
    
    
end

% store path
spath = [];
if success
    path_reverse = [2];
    i = 2;
    while true
       if all(milestones(i, :) == start)
           break
       end
       i = parent(i);
       path_reverse = [path_reverse i];
    end
    spath = flip(path_reverse);
end



% ------end of your optimized algorithm-------
dt = toc;

figure(2); hold on;
plot(milestones(:,1),milestones(:,2),'m.');
if (~isempty(edges))
    line(edges(:,1:2:3)', edges(:,2:2:4)','Color','magenta')
end
if (~isempty(spath))
    for i=1:length(spath)-1
        plot(milestones(spath(i:i+1),1),milestones(spath(i:i+1),2), 'go-', 'LineWidth',3);
    end
end
str = sprintf('Q3 - %d X %d Maze solved in %f seconds', row, col, dt);
title(str);

print -dpng assignment1_q3.png

