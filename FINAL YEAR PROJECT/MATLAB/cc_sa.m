 clc;clear;close all;
% N =/ 1-63,64-67,68-74,75-82
trials=5;
sd=0;
N_v=50;
cost_comparison=zeros(N_v,4,trials);
for tr=1:trials
tic;
disp("Trial: "+tr);

for N=N_v:N_v

disp(N);
sd=tr;
nodes=generateNodes(N,tr);

origin=[0.5 0.5];

nodes=[origin;nodes];



%Simulated Annealing
path_sa=sa(nodes);
cost_comparison(N,4,tr)=getCost(nodes,path_sa);
end
toc;
end


% i = 1:N_v;
% mc1=mean(cost_comparison(i,1,:),3);
% mc2=mean(cost_comparison(i,2,:),3);
% mc3=mean(cost_comparison(i,3,:),3);
% mc4=mean(cost_comparison(i,4,:),3);
% plot(i, mc1, 'LineStyle', '-',Color='#ffcc00');
% 
% hold on
% plot(i,mc2,'LineStyle', '-',Color='green');
% plot(i,mc3,'LineStyle', '-',Color='#A020F0');
% plot(i,mc4,'LineStyle', '-',Color='red');
% 
%  ylim([0,ceil(cost_comparison(end,1,1))]);
% set(gca,'xtick', 1:N_v);set(gca,'ytick', 1:ceil(cost_comparison(end,1,1)));
toc;% OPTIMISATION FUNCTIONS

function path_sa=sa(nodes)
    
min_cost_sa=inf;

% Set initial temperature and other parameters
T0=1;
coolingRate = 0.99; % cooling rate = 91.25%
numIterations = 6; % number of iterations
numAttempts = 10; % number of attempts to find new solution at each temperature




sa_analysis=NaN(numIterations*numAttempts,3);
% Loop through iterations
for i = 1:numIterations
    % Cooling schedule
    T = T0*(coolingRate^i);
    
 
 
    currentPath=randperm(length(nodes))';
    currentPath(currentPath==1)=[];
    currentPath = [1;currentPath;1]; % choose random path
    currentDist = getCost(nodes,currentPath); 
    
    for j = 1:numAttempts
           sa_analysis((i-1)*numAttempts+j,1)=T;
          
          
            
          

        % Generate a new solution
         path_sa = two_Opt(nodes,currentPath);
         newDist = getCost(nodes, path_sa);
         sa_analysis((i-1)*numAttempts+j,2)=newDist;
         deltaE = newDist - currentDist;

          acceptanceProbability = exp(-deltaE/T); % If new solution is worse, accept it with probability e^(-deltaE/T)
             
            if rand() > acceptanceProbability
                 path_sa=random_swap(path_sa);
 
             end
        if deltaE <= 0
            currentPath = path_sa;
            currentDist = newDist;
       
        else
            acceptanceProbability = exp(-deltaE/T); % If new solution is worse, accept it with probability e^(-deltaE/T)
            if rand() < acceptanceProbability
                currentPath = path_sa;
                currentDist = newDist;
            end
        end
     end
    if(getCost(nodes,path_sa)<min_cost_sa)
    min_cost_sa=getCost(nodes,path_sa);
    sa_analysis((i-1)*numAttempts+j,3)=min_cost_sa;
    min_path=path_sa;
    end
end
sa_analysis(end,end)=min_cost_sa;
peak=sa_analysis(:,3);
peak_idx=isnan(peak);
peak=peak(~isnan(peak));
disp(length(peak)-1);
figure;
plot(sa_analysis(:,2))
hold on
scatter(1:length(sa_analysis),sa_analysis(:,3),'red')


stairs(find(peak_idx==0),peak,'LineStyle','-',Color='green');
path_sa=min_path;
end

function path_eg=christofides(nodes)

    distmat_mst=getDistmat(nodes);
G=graph(distmat_mst);
T=minspantree(G);
    degrees=degree(T);
odd_nodes=find(mod(degrees, 2) == 1); %finds all the nodes which have an odd number of vertices in the MST
distmat_match=getDistmat(nodes(odd_nodes,:));

distmat_match(distmat_match==0)=inf;
matching = min_perfect_matching(distmat_match);

matching=[odd_nodes(1:length(matching)),odd_nodes(matching)];

for i=1:length(matching)
    for j=1:length(matching)
        if (matching(i,1)==matching(j,2) && matching(i,2)==matching(j,1))
            matching(j,:)=[0 0];
        end
    end
end
matching(all(~matching,2), : ) = [];


eulerian_graph=[T.Edges.EndNodes];
eulerian_graph=[eulerian_graph,ones(length(eulerian_graph),1)];


matching=[matching,ones(size(matching,1),1)];

for i=1:size(matching,1)

    for j=1:size(eulerian_graph,1)
        if (matching(i,1)==eulerian_graph(j,1) && matching(i,2)==eulerian_graph(j,2)) || (matching(i,1)==eulerian_graph(j,2) && matching(i,2)==eulerian_graph(j,1))
        eulerian_graph(j,3)=eulerian_graph(j,3)+1;
        matching(i,3)=0;
        end
    end
    
end

for i=1:size(matching,1)
    if(matching(i,3)==1)
        eulerian_graph=[eulerian_graph;matching(i,:)];
    end
end




distmat_eg=zeros(length(nodes));
    
for i=1:length(eulerian_graph)
    distmat_eg(eulerian_graph(i,1),eulerian_graph(i,2))=eulerian_graph(i,3);
    distmat_eg(eulerian_graph(i,2),eulerian_graph(i,1))=eulerian_graph(i,3);
end

path_eg=get_eulerian_path(distmat_eg,nodes);
path_eg=unique(path_eg,"stable");
path_eg=[path_eg,path_eg(1)];
path_eg=path_eg';
end



function path = get_eulerian_path(adj_matrix,nodes)
    % Create adjacency matrix from distance matrix
%     adj_matrix = adj_matrix;
%     
     % Find any vertex with odd degree
%      odd_vertices = find(mod(sum(adj_matrix), 2));
%      start_vertex = odd_vertices(1);
%      if isempty(start_vertex)
        start_vertex = 1; % If no odd degree vertices, start anywhere
%      end
    
    % Initialize stack and path
    stack = start_vertex;
    path = [];
    
    while ~isempty(stack)
        % Take the last vertex on the stack
        current_vertex = stack(end);
        
        % Find all unvisited edges from the current vertex
        unvisited_edges = find(adj_matrix(current_vertex, :) > 0);
        visited_edges = find(adj_matrix(current_vertex, :) == 0);
        
        if ~isempty(unvisited_edges)
            % Visit the first unvisited edge
            next_vertex = nearest_unvisited_node(current_vertex,unvisited_edges,nodes);%unvisited_edges(1);
            
            % Mark the edge as visited
            adj_matrix(current_vertex, next_vertex) = adj_matrix(current_vertex, next_vertex)-1;
            adj_matrix(next_vertex, current_vertex) =adj_matrix(next_vertex, current_vertex)-1 ;
            
            % Push the next vertex onto the stack
            stack(end+1) = next_vertex;
        else
            % Pop the current vertex off the stack
            stack(end) = [];
            
            % Add the current vertex to the path
            path(end+1) = current_vertex;
        end
    end
    
    % Reverse the path to get the correct order
    path = fliplr(path);
end


function [MST] = prims_algorithm(nodes, distmat)

N=length(nodes);

visited = false(1, N); % initialize all nodes as unvisited
MST = zeros(N-1, 2); % initialize the minimum spanning tree
total_cost = 0; % initialize the total cost

% choose the first node as the starting point
current_node = 1;
visited(current_node) = true;

for i = 1:N-1 % repeat N-1 times
    % find the minimum edge that connects the current node to an unvisited node
    min_dist = inf;
    min_node = 0;
    for j = 1:N
        if visited(j)
            % check all neighbors of the current node
            for k = 1:N
                if ~visited(k) && distmat(j,k) < min_dist
                    min_dist = distmat(j,k);
                    current_node=j;
                    min_node = k;
                end
            end
        end
    end

    % add the minimum edge to the minimum spanning tree
    MST(i,:) = [current_node, min_node];
    %greenLine(nodes(current_node,:),nodes(min_node,:));
    % total_cost = total_cost + min_dist;
    visited(min_node) = true;
    current_node = min_node;
end

end

function min_path=two_Opt(nodes,path)
cost_2opt=getCost(nodes,path);
min_path=path;

for i=1:length(path)-3 % first edge
    for j=i+2:length(path)-1    % second edge will always be after the first edge

        path=swapEdges(path,i,j);

        c=getCost(nodes,path);
        if(c<cost_2opt)
            min_path=path;
            cost_2opt=c;


        else

            path=swapEdges(path,i,j);
        end

    end

end

end
% FUNCTIONS
function min_node=nearest_unvisited_node(current_node,unvisited_edges,nodes)
    min_dist=inf;
    min_node=1;
    for i=1:length(unvisited_edges)
        
            if(dist(nodes(current_node,:),nodes(unvisited_edges(i),:))<min_dist)
                min_dist=dist(nodes(current_node,:),nodes(unvisited_edges(i),:));
                min_node=unvisited_edges(i);
            end
    end
    
end

function distmat=MST(nodes)

distmat=zeros(length(nodes),length(nodes));
end

function plotPathDirection(nodes,path)

for i=1:length(path)-1
    p1=nodes(path(i),:);
    p2=nodes(path(i+1),:);  dp=p2-p1;
    %dp=p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),"off",'black','filled','ShowArrowHead','on');
    %quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),"off","black");

end
end

function distmat=getDistmat(nodes)
distmat=zeros(length(nodes),length(nodes));
for i=1:length(nodes) % matrix of distances between two nodes
    p1=nodes(i,:);
    for j=1:i
        p2=nodes(j,:);
        distmat(i,j)=dist(p1,p2);
        distmat(j,i)=dist(p2,p1);
        if(dist(p1,p2)==0)
            distmat(i,j)=inf;
        end
    end
end


end

function path=swapEdges(path,p1,p2)
path=path(1:end-1);


path(p1+1:p2)=flip(path(p1+1:p2));
path(end+1)=path(1);


end

function cost=getCost(nodes,path)
cost=0;
for i=1:length(path)-1
    cost=cost+dist(nodes(path(i),:),nodes(path(i+1),:));
end
end

function out=isInArray(arr,a)
out=false
for i=1:length(arr)
    if(arr(i)==a)
        out=true;
    end
end
end

function out=areIndirectlyConnected(distmat,p1,p2)
out=false;
path=getPath(distmat,p1);
if(path(end)==p2)
    out=true;
end
end

function path=getPath(distmat,start)
distmat1=distmat;
p1=start;
path=p1;

for k=1:length(distmat1(1,:))
    for i=1:length(distmat1(:,1))

        if(distmat1(p1,i)==-1)
            if(length(path)>1 && i==path(end-1))

            else
                p1=i;
                distmat1(p1,i)=0;
                distmat1(i,p1)=0;
                path=[path;i];
                if(i==start)

                    return;
                end
                break;

            end




        end
    end
end

if(~path(end)==path(1))

    path(end+1)=path(1);
end
end

function out=isIncomplete(traversed,x)
out=true;
temp=traversed;


if(temp(x(1))==1 && temp(x(2))==1)
    temp(x(1))=temp(x(1))+1;
    temp(x(2))=temp(x(2))+1;

    for i=1:length(temp)
        if(temp(i)~=2)
            out=false;

            break;
        else

        end
    end
end
end

function greenLine(p1,p2,c)
deleteLine(p1,p2)
X = [p1(:,1) p2(:,1)] ;
Y = [p1(:,2) p2(:,2)] ;
Z = [p1(:,3) p2(:,3)] ;
plot3(X',Y',Z','green');
hold on
plot3(X',Y',Z','.',Color='black');
end

function redLine(p1,p2,c)
deleteLine(p1,p2)
X = [p1(:,1) p2(:,1)] ;
Y = [p1(:,2) p2(:,2)] ;
Z = [p1(:,3) p2(:,3)] ;
plot3(X',Y',Z','red');
hold on
plot3(X',Y',Z','.',Color='black');
end
function line(p1,p2)
X = [p1(:,1) p2(:,1)] ;
Y = [p1(:,2) p2(:,2)] ;
Z = [p1(:,3) p2(:,3)] ;
plot3(X',Y',Z','black');
hold on
plot3(X',Y',Z','.',Color='black');
end

function deleteLine(p1,p2)
plot([ p1(1) p2(1)], [p1(2) p2(2)],'LineStyle','none');
end
function f=min_nz(arr) % Return the index and value minimum non-zero value in a matrix
m=inf; f=[1 1 0]; % [1 1] is the default value
for i=1:length(arr(:,1))
    for j=1:length(arr(1,:))
        if(arr(i,j)<m && arr(i,j)>0)
            m=arr(i,j);
            f=[i,j,m];

        end
    end
end
if (m==inf)
    f=[i,j,0];
end

end




function nodes=generateNodes(N,r)
rng(r); % Could be any random natural number
seed=rng;
nodes=rand(N,2);
rng(seed);
end
function dis=dist(a,b) % calculates the distance between the two points a and b on a 2d plane


dis=sqrt((b(1)-a(1))^2+(b(2)-a(2))^2);


end

    function swapped_path = random_swap(path)
% This function randomly swaps two elements in the input path array, except for the first and last element.

% Generate two random indices within the bounds of the array, but not the first or last element
idx1 = randi([2,length(path)-1]);
idx2 = randi([2,length(path)-1]);

% Ensure the two indices are different
while idx2 == idx1
    idx2 = randi([2,length(path)-1]);
end

% Swap the values at the two indices
temp = path(idx1);
path(idx1) = path(idx2);
path(idx2) = temp;

% Return the swapped path
swapped_path = path;

end