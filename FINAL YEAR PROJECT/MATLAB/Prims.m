N=20;

rng(144); % Could be any random natural number
    seed=rng;
    nodes=rand(N,2); 
    rng(seed);
origin=[0.5 0.5 0.5];

% nodes=[origin;nodes];

distances = getDistmat(nodes);

prims_algorithm(nodes,distances)

function [MST, total_cost] = prims_algorithm(nodes, distances)
% Prim's algorithm for minimum spanning tree
% nodes: an Nx2 matrix containing the x and y coordinates of the nodes
% distances: an NxN matrix containing the pairwise distances between the nodes

N = size(nodes, 1); % number of nodes
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
                if ~visited(k) && distances(j,k) < min_dist
                    min_dist = distances(j,k);
                    min_node = k;
                end
            end
        end
    end
    
    % add the minimum edge to the minimum spanning tree
    MST(i,:) = [current_node, min_node];
    total_cost = total_cost + min_dist;
    visited(min_node) = true;
    current_node = min_node;
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


function dis=dist(a,b) % calculates the distance between the two points a and b on a 3d plane


dis=sqrt((b(1)-a(1))^2+(b(2)-a(2))^2+(b(3)-a(3))^2);

% if(b(3)<a(3))
%     dis=sqrt((b(1)-a(1))^2+(b(2)-a(2))^2);
% end

end