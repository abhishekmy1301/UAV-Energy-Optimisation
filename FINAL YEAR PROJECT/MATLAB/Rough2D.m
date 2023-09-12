clc;clear;close all;
tic

N=49; %Number of nodes excluding the origin
rng(17) % Could be any random natural number
seed=rng;
nodes=rand(N,2); % Coordinates of nodes
rng(seed);
%Initalizing coordinates

xvals=nodes(:,1);
yvals=nodes(:,2);
origin=[0.5 0.5];
nodes=[origin ;nodes ]; % adding origin coordinates to the front and the back of the array to complete the path
% now length of nodes is N+2


st=60.*ones(N+1,1); % service time of UAV at each node in seconds
st(1)=0; % service time at depot is 0s



figure('Name','Minimum Spanning Tree - Prim''s Algorithm','NumberTitle','off');
plotNodes(nodes);
distmat_mst=getDistmat(nodes);
prims_algorithm(nodes,distmat_mst);  

% Original Path
figure('Name','Random Order of Nodes','NumberTitle','off');
plotNodes(nodes);
path_rnd=1:N+1;
path_rnd(end+1)=path_rnd(1);
cost_rnd1=getCost(nodes,path_rnd);
fprintf('Cost of Random Path: %f\n',cost_rnd1);
plotPathDirection(nodes,path_rnd);
% Original Path - 2Opt
figure('Name','Random Order - 2-Opt Optimal','NumberTitle','off');

path_rnd=two_Opt(nodes,path_rnd);
clf;plotNodes(nodes);plotPathDirection(nodes,path_rnd);
cost_rnd2=getCost(nodes,path_rnd);
while(~(cost_rnd2==cost_rnd1))
    cost_rnd2=getCost(nodes,path_rnd);
    path_rnd=two_Opt(nodes,path_rnd);
    cost_rnd1=getCost(nodes,path_rnd);
    clf;plotNodes(nodes);plotPathDirection(nodes,path_rnd);
end
    

plotNodes(nodes);

plotPathDirection(nodes,path_rnd);
fprintf('Cost After 2-Opt Optimization: %f\n',getCost(nodes,path_rnd));

%Nearest Neighbour Heuristic
figure('Name','Nearest Neighbour Heuristic','NumberTitle','off');

% nodes(end,:) = [];
distmat_nnh=getDistmat(nodes);
plotNodes(nodes,origin);
p1=origin;
current=1;
for i=1:N
    x=min_nz(distmat_nnh(current,:));
    x(1)=current;
    p2=nodes(x(2),:);
    dp=p2-p1;
    quiver(p1(1),p1(2),dp(1),dp(2),0,'black');

    p=distmat_nnh(current,:);
    p(p>0)=0;
    distmat_nnh(current,:)=p;
    p=distmat_nnh(:,current);
    p(p>0)=0;
    distmat_nnh(:,current)=p;
    distmat_nnh(current,x(2))=-1;
    distmat_nnh(x(2),current)=-1;

    p1=p2;
    current=x(2);

end
current=1;
p=distmat_nnh(current,:);
p(p>0)=0;
distmat_nnh(current,:)=p;
p=distmat_nnh(:,current);
p(p>0)=0;
distmat_nnh(:,current)=p;
distmat_nnh(current,x(2))=-1;
distmat_nnh(x(2),current)=-1;


p1=p2;p2=origin; dp=p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),0,'black'); % Adding last edge manually
path_nnh=getPath(distmat_nnh,1);

fprintf('Cost of NNH: %f\n',getCost(nodes,path_nnh));






% 2-Opt NNH
figure('Name','Nearest Neighbour Heuristic - 2-Opt Optimal','NumberTitle','off');
plotNodes(nodes,origin);
cost_nnh1=getCost(nodes,path_nnh);
path_nnh_2opt=two_Opt(nodes,path_nnh);

 plotNodes(nodes);

plotPathDirection(nodes,path_nnh_2opt);


fprintf('Cost After 2-Opt Optimization: %f\n',getCost(nodes,path_nnh_2opt));
plotPathDirection(nodes,path_nnh_2opt);





% Greedy Heuristic Approach
figure('Name','Greedy Heuristic','NumberTitle','off');


plotNodes(nodes,origin);


distmat_gh=getDistmat(nodes);



traversed=zeros(N+1,1);
path_gh=0;

for j=1:(N+1)^2
    x=min_nz(distmat_gh);
    if(x(3)==0)
        break;
    end

    if(traversed(x(1))<2 && traversed(x(2))<2 && ~areIndirectlyConnected(distmat_gh,x(1),x(2)) )




        p1=nodes(x(1),:);p2=nodes(x(2),:);

        line(p1,p2);
        path_gh=path_gh+dist(p1,p2);

        traversed(x(1))=traversed(x(1))+1;
        traversed(x(2))=traversed(x(2))+1;

        distmat_gh(x(1),x(2))=-1;
        distmat_gh(x(2),x(1))=-1;

    else
        distmat_gh(x(1),x(2))=0;
        distmat_gh(x(2),x(1))=0;
    end





end


a=1; % Adding last edge manually
for a=1:length(traversed)
    if(traversed(a)==1)
        p1=nodes(a,:);
        break;
    end
end
for i=(a+1):length(traversed)
    if(traversed(i)==1)
        p2=nodes(i,:);
        distmat_gh(a,i)=-1;
        distmat_gh(i,a)=-1;
    end
end
line(p1,p2);
path_gh=path_gh+dist(p1,p2);

fprintf('Cost of GH: %f\n',path_gh);

hold off




figure('Name','Greedy Heuristic - 2-Opt Optimal','NumberTitle','off');
path_gh_2opt=getPath(distmat_gh,1); %  greedy heuristic


 path_gh_2opt=two_Opt(nodes,path_gh_2opt); 



plotNodes(nodes,origin);

plotPathDirection(nodes,path_gh_2opt);
 fprintf('Cost After 2-Opt Optimization: %f\n',getCost(nodes,path_gh_2opt));




%CHRISTOFIDES ALGORITHM
figure('Name','Minimum Spanning Tree - Prim''s Algorithm','NumberTitle','off');
plotNodes(nodes);
distmat_mst=getDistmat(nodes);
prims_algorithm(nodes,distmat_mst);  




path_eg=christofides(nodes);

figure('Name','Christofides Algorithm','NumberTitle','off');
plotNodes(nodes);
plotPathDirection(nodes,path_eg);
fprintf('Cost of CA: %f\n',getCost(nodes,path_eg));
path_eg=two_Opt(nodes,path_eg);
figure('Name','Christofides Algorithm - 2-Opt Optimal','NumberTitle','off');
plotNodes(nodes);
plotPathDirection(nodes,path_eg);
fprintf('Cost after 2-Opt: %f\n',getCost(nodes,path_eg));




%Simulated Annealing
path_sa=sa(nodes);
figure('Name','Simulated Annealing','NumberTitle','off');

plotNodes(nodes);plotPathDirection(nodes,path_sa);
fprintf('Cost of SA: %f\n',getCost(nodes,path_sa));







% POWER CONSUMPTION MODEL

%distance=velocity*time energy= power*(distance/velocity)=power*(m/(m/s)/3600)Wh
st=60.*ones(N,1); % service time of UAV at each node in seconds








% CLUSTERING
capacity=320; %Wh
max_cluster_cost=inf; %Wh
k=1;
while max_cluster_cost>capacity  
    str="K-Means - Number of clusters: "+num2str(k);
figure('Name',str,NumberTitle='off');
 
nodes_wo=nodes(2:end,:);
[idx,centres]=kmeans(nodes_wo,k);
plotNodes(nodes);
    colors = hsv(k); % generate k unique colors
max_cluster_cost=0;
final_cost=0;
    for i = 1:k
        
        scatter(nodes_wo((idx == i) , 1), nodes_wo((idx == i), 2),[], colors(i, :), '.');
        hold on;
        cluster_path=[1;find(idx==i)+1;1];
         cluster_path_order=[sa(nodes(cluster_path,:))];
         cluster_path_order=cluster_path_order(1:end-1);
        cluster_path = [cluster_path(cluster_path_order);1];
        cluster_cost=getCost(nodes,cluster_path);
        cluster_cost=cluster_cost+sum(st(idx==i))*1;
        final_cost=final_cost+cluster_cost;
        if(cluster_cost>max_cluster_cost)
            max_cluster_cost=cost_wh(getCost(nodes,cluster_path),2);
        end
        plotPathDirection(nodes,cluster_path);
        
    end
    
k=k+1;
end

fprintf("Final Cost in Wh = %f\n",cost_wh(final_cost,2));
toc


% FUNCTIONS

function min_path = three_Opt(nodes, path)
    cost_3opt = getCost(nodes, path);
    min_path = path;

    for i = 1:length(path)-5   % first edge
        for j = i+2:length(path)-3 % second edge will always be after the first edge
            for k = j+2:length(path)-1 % third edge will always be after the second edge

                path_1 = swapEdges(path, i, j-1);
                path_2 = swapEdges(path, j, k-1);
                path_3 = swapEdges(path, i, k-1);

                c_1 = getCost(nodes, path_1);
                c_2 = getCost(nodes, path_2);
                c_3 = getCost(nodes, path_3);

                [c_min, ind] = min([c_1, c_2, c_3]);

                if(c_min < cost_3opt)
                    cost_3opt = c_min;
                    if(ind == 1)
                        min_path = path_1;
                    elseif(ind == 2)
                        min_path = path_2;
                    else
                        min_path = path_3;
                    end
                end
            end
        end
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
    dp=p2-p1;quiver(p1(1),p1(2),dp(1),dp(2),"off",'black','filled','ShowArrowHead','on');

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

function plotPath(nodes,path)
clf;
plotNodes(nodes,[0.5 0.5]);
for i=1:(length(path)-1)
    p1=nodes(path(i),:);
    p2=nodes(path(i+1),:);

    line(p1,p2);

end
end

function greenLine(p1,p2)
deleteLine(p1,p2)
plot([ p1(1) p2(1)], [p1(2) p2(2)],'green');
end

function line(p1,p2)
plot([ p1(1) p2(1)], [p1(2) p2(2)],'black');
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

function redLine(p1,p2)
deleteLine(p1,p2)
plot([ p1(1) p2(1)], [p1(2) p2(2)],'red');
end

function plotNodes(nodes,origin)
clf;
scatter(nodes(:,1),nodes(:,2),10,'blue','filled');
hold on
scatter(nodes(1,1),nodes(1,2),250,'black','x');
text(nodes(1,1),nodes(1,2),'Depot','HorizontalAlignment','left','FontSize',10,'VerticalAlignment','top');
for i=2:(length(nodes(:,1)))
    text(nodes(i,1),nodes(i,2),num2str(i),'HorizontalAlignment','left','FontSize',10,'VerticalAlignment','top');

end
xlim([0 1]); ylim([0,1]);
end

function plotClusters(nodes,idx,k)

    colors = hsv(k); % generate k unique colors
for i = 1:k
    scatter(nodes((idx == i) , 1), nodes((idx == i), 2),[], colors(i, :), '.');
    hold on;
    cluster_path=[1;find(idx==i);1];
    plotPathDirection(nodes,two_Opt(nodes,cluster_path));
   
end

end

function out=isTraversed(traversed,a)% Checks if a node has already been traversed in a path
out = false;
for i=1:length(traversed(:,1))
    if(a==traversed(i,:))
        out = true;
    end
end
end

function dis=dist(a,b) % calculates the distance between the two points a and b on a 2d plane

dis=sqrt((b(1)-a(1))^2+(b(2)-a(2))^2);
end

function newPath = generateNewPath(path)
    % Generates a new solution by swapping two random cities
    N = length(path);
    i = randi(N);
    j = randi(N);
    while i == j
        j = randi(N);
    end
    newPath = path;
    newPath([i j]) = newPath([j i]);
end

function P_T=power_uav(V)
U_tip=200;%m/s
phi=1.225; %kg/m^3
A=0.79; %rotor disk area m^2
d_o=0.3;%Fuselage drag ratio

s_hat=0.05;%rotor solidity
% Aircraft forward sped in m/s


small_sigma=0.012;%profile drag coefficient
s=0.05;%rotor solidity
omega=400;%Blade angular velocity in radians/second
R=0.5;%Rotor radius in m
m=3.29;%mass kg
g=9.81; %gravity m/s^2
k=2;%Thrust to wight ratio T/W
v_o=sqrt((m*g)/(2*phi*A));%motor induced velocity m/s
P_o=(small_sigma/8)*phi*s*A*(omega^3)*(R^3);

P_i=(1+k)*(m*g)^1.5/sqrt(2*phi*A);
P_T=P_o*(1+3*(V^2)/U_tip^2)+ P_i*sqrt((sqrt(1+V^4/(4*v_o^4)))-V^2/(2*v_o^2))+0.5*d_o*phi*s*A*V^3;
end

function cost_wh=cost_wh(cp_cost,V)
%    cp=path; %chosen path
% cp_cost=getCost(nodes,cp);%

P_T=power_uav(V); % power to travel at V m/s

total_st=(cp_cost*1000/3)/3600;
cost_wh=P_T*(total_st);
end

% OPTIMISATION FUNCTIONS



function path_sa=sa(nodes)
    
min_cost_sa=inf;
% Set initial temperature and other parameters
T0=1;
coolingRate = 0.95; % cooling rate
numIterations = 10; % number of iterations
numAttempts = 10; % number of attempts to find new solution at each temperature

% Generate initial solution
currentPath = [1:length(nodes),1]'; % choose original path

% Evaluate objective function
currentDist = getCost(nodes,currentPath); % calculate the total tour distance

% Loop through iterations
for i = 1:numIterations
    % Cooling schedule
    T = T0 * (coolingRate ^ i);
    
    % Try numAttempts times to find a new solution
    for j = 1:numAttempts
        % Generate a new solution
        path_sa = two_Opt(nodes,currentPath);
        
        % Evaluate objective function
        newDist = getCost(nodes, path_sa);
        
        % Calculate delta E
        deltaE = newDist - currentDist;
        
        % If new solution is better, accept it
        if deltaE < 0
            currentPath = path_sa;
            currentDist = newDist;
        % If new solution is worse, accept it with probability e^(-deltaE/T)
        else
            acceptanceProbability = exp(-deltaE/T);
            if rand() < acceptanceProbability
                currentPath = path_sa;
                currentDist = newDist;
            end
        end
    end
end

% find(newPath==1);
% newPath=circshift(newPath,N+2-(find(newPath==1)));
if(getCost(nodes,path_sa)<min_cost_sa)
    min_cost_sa=getCost(nodes,path_sa);
    
end
path_sa=[path_sa];
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
    redLine(nodes(matching(i,1),:),nodes(matching(i,2),:));
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




function [MST] = prims_algorithm(nodes, distmat)

N=length(nodes);
nodes=nodes;
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
    greenLine(nodes(current_node,:),nodes(min_node,:));
    % total_cost = total_cost + min_dist;
    visited(min_node) = true;
    current_node = min_node;
end

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