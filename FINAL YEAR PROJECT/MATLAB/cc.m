close all;
% N =/ 1-63,64-67,68-74,75-82
trials=40;

N_v=100 ;
cost_comparison=zeros(N_v,4,trials);
for tr=1:trials
tic;
disp("Trial: "+tr);

for N=2 :N_v

disp(N);
origin=[0.5 0.5];
nodes=[origin;generateNodes(N,tr+50)];


 








%Nearest Neighbour Heuristic


distmat_nnh=getDistmat(nodes);

p1=origin;
current=1;
for i=1:N
    
    x=min_nz(distmat_nnh(current,:));
    x(1)=current;
    p2=nodes(x(2),:);
    

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
p1=p2;p2=origin; dp=p2-p1;
path_nnh=getPath(distmat_nnh,1);


cost_nnh=getCost(nodes,path_nnh);

cost_comparison(N,1,tr)=getCost(nodes,path_nnh);




% Greedy Heuristic Approach




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

       % line(p1,p2);
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

path_gh=path_gh+dist(p1,p2);
cost_gh=path_gh;
path_gh=getPath(distmat_gh,1);





cost_comparison(N,2,tr)=cost_gh;

% path_gh1=two_Opt(nodes,path_gh);
% while(~(path_gh==path_gh1))
%     path_gh=two_Opt();
% end


%CHRISTOFIDES ALGORITHM


path_eg= christofides(nodes);

cost_comparison(N,3,tr)=getCost(nodes,path_eg);


%Simulated Annealing
path_sa=sa(nodes);

cost_comparison(N,4,tr)=getCost(nodes,path_sa);
end
toc;
end


i = 1:N_v;
mc1=mean(cost_comparison(i,1,:),3);
mc2=mean(cost_comparison(i,2,:),3);
mc3=mean(cost_comparison(i,3,:),3);
mc4=mean(cost_comparison(i,4,:),3);
%plot(i, mc1, 'LineStyle', '-',Color='#ffcc00',LineWidth=2.5);
plot(i,mc4,'LineStyle', '-',Color='red');
hold on


N_v = 100; 
i = 1:N_v;

y = mc4'; 
x = 1:numel(y);
std_dev = 1;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

% Change fill color to transparent yellow and only use one shaded region
fill(x2, inBetween, [1, 0.6, 0.6],'facealpha',0.5,'edgecolor','none');



plot(i,mc2,'LineStyle', '-',Color='green');
plot(i,mc3,'LineStyle', '-',Color='#A020F0');
plot(i,mc4,'LineStyle', '-',Color='red');

% ylim([0,ceil(cost_comparison(end,1,1))]);
set(gca,'xtick', 1:N_v);set(gca,'ytick', 1:ceil(cost_comparison(end,1,1)));
toc;


% OPTIMISATION FUNCTIONS


function path=two_Opt_Optimal(nodes,path)
    cost=getCost(nodes,path);
    path=two_Opt(nodes,path);
    cost2=getCost(nodes,path);
    while(~(cost2==cost))
        cost2=getCost(nodes,path);
        path=two_Opt(nodes,path);
        cost=getCost(nodes,path);
    end
    
end

function path_sa=sa(nodes)
    
min_cost_sa=inf;

% Set initial temperature and other parameters
T0=1;
coolingRate = 0.9125; % cooling rate = 91.25%
numIterations = 9; % number of iterations
numAttempts = 6; % number of attempts to find new solution at each temperature





% Loop through iterations
for i = 1:numIterations
    % Cooling schedule
    T = T0*(coolingRate^i);
    
 

    currentPath=randperm(length(nodes))';
    currentPath(currentPath==1)=[];
    currentPath = [1;currentPath;1]; % choose random path
    currentDist = getCost(nodes,currentPath); 
    for j = 1:numAttempts
        % Generate a new solution
        path_sa = two_Opt(nodes,currentPath);
        newDist = getCost(nodes, path_sa);
        deltaE = newDist - currentDist;
        if deltaE < 0
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
    min_path=path_sa;
    end
end




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

function [ indices, cost ] = min_perfect_matching( G )
% Minimum Perfect Matching Tool
%
% Synopsis
%   [indices, cost] = min_perfect_matching( G ) 
%   
% Description
%   Function to solve the Minimum Perfect Matching on non-biparite graphs
%   problem using Integer linear programming.
%
%   Returns vector of matched indices and cost of the match. Requires
%   symmetric adjacent matrix of even rank.
% 
% Author: Vojtech Knyttl, knyttvoj@fel.cvut.cz

    %% dimensions
    if( ~isequal(G, G') )
      e = MException('PerfectMatching:MxNotSymetric', 'Matrix is not symmetric.' );
      throw(e);
    end;

    len = uint16(length( G ));
    if( mod(len,2) == 1 )
      e = MException('PerfectMatching:MxNotEvenRank', 'Matrix must be of even rank.' );
      throw(e);
    end;
    
    llen = (len*len)-len;
 
    %% function to minimize
    f = zeros(1,llen);
    for i=1:len-1
      f((i-1)*len+1:i*len)=G(i,:);
    end;
       
    %% respect matrix
    A = zeros(len,llen);
    for i=1:len
      for j=1:len*len-len
        idiv = idivide(j-1,len,'floor')+1;
        if idiv == i
          if and(mod(j-1,len)>=i,mod(j,len)~=i)
            A(i,j) = 1;
          end;
        elseif idiv < i
          if mod(j-1,len)+1==i
            A(i,j) = 1;  
          end;
        end;
      end;
    end;
    b = ones( len, 1 );
    
    %% removing zero columns
    remove_cols = find(all(A==0));
    f(:,remove_cols)=[];
    A(:,remove_cols)=[];
    
    %% glpk ilinprog settings
    sense=1;                            % minimization
    ctype = repmat( 'S', 1, len );      % equalities
    lb = zeros( llen/2, 1 );            % lower bound
    ub = ones( llen/2, 1 );             % upper bound
    i = 1:llen/2; 
    e=2^-24;
    % vartype = repmat( 'I', 1, llen/2 ); % integral vals
   % param.msglev = 1;  
   % param.itlim = 100;
    
    xmin = IP1(f,[],[],A,b,lb,ub,i,e);

    %xmin = glpk(f,A,b,lb,ub,ctype,vartype,sense,param);
    
    %% adding remove columns
    match = find(xmin==1)';
    
    for i=remove_cols
      match(match>=i)=match(match>=i)+1;
    end;
    
    %% forming the result
    indices = zeros(1,len); cost = 0;
    for i=match
      x = idivide(i-1,len,'floor')+1;
      y = mod(i-1,len)+1;
      indices(x)=y;
      indices(y)=x;
      cost = cost + G(x,y);
    end;
end

% By Sherif A. Tawfik, Faculty of Engineering, Cairo University
% [x,val,status]=IP1(f,A,b,Aeq,beq,lb,ub,M,e)
% this function solves the following mixed-integer linear programming problem
%   min f*x
%  subject to
%        A*x <=b
%        Aeq * x = beq   
%        lb <= x <= ub
%        M is a vector of indeces for the variables that are constrained to be integers
%        e is the integarilty tolerance
% the return variables are :
% x : the solution
% val: value of the objective function at the optimal solution
% status =1 if successful
%        =0 if maximum number of iterations reached in he linprog function
%        =-1 if there is no solution
% Example:
%        maximize 17 x1 + 12 x2 
%        subject to
% 	             10 x1 + 7 x2 <=40
%                   x1 +   x2 <= 5
%                   x1, x2 >=0 and are integers
% f=[-17, -12]; %take the negative for maximization problems
% A=[ 10  7; 1 1];
% B=[40; 5];
% lb=[0 0];
% ub=[inf inf];
% M=[1,2];
% e=2^-24;
% [x v s]= IP(f,A,B,[],[],lb,ub,M,e)
function [x,val,status]=IP1(f,A,b,Aeq,beq,lb,ub,M,e)
options = optimset('display','off');
bound=inf; % the initial bound is set to +ve infinity
[x0,val0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 
[x,val,status,b]=rec(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound);
end % a recursive function that processes the BB tree 
function [xx,val,status,bb]=rec(f,A,b,Aeq,beq,lb,ub,x,v,M,e,bound) 
options = optimset('display','off');
% x is an initial solution and v is the corressponding objective function value
% solve the corresponding LP model with the integarily constraints removed
[x0,val0,status0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 

% if the solution is not feasible or the value of the objective function is
% higher than the current bound return with the input intial solution
if status0<=0 | val0 > bound  
    xx=x; val=v; status=status0; bb=bound;
    return;
end
% if the integer-constraint variables turned to be integers within the
% input tolerance return
ind=find( abs(x0(M)-round(x0(M)))>e ); 
if isempty(ind)
    status=1;        
    if val0 < bound    % this solution is better than the current solution hence replace
        x0(M)=round(x0(M));
        xx=x0;        
        val=val0;
        bb=val0;
    else
        xx=x;  % return the input solution
        val=v;
        bb=bound;
    end
    return
end
% if we come here this means that the solution of the LP relaxation is
% feasible and gives a less value than the current bound but some of the
% integer-constraint variables are not integers. 
% Therefore we pick the first one that is not integer and form two LP problems
% and solve them recursively by calling the same function (branching)
% first LP problem with the added constraint that Xi <= floor(Xi) , i=ind(1)
br_var=M(ind(1));
br_value=x(br_var);
if isempty(A)
    [r c]=size(Aeq);
else
    [r c]=size(A);
end
A1=[A ; zeros(1,c)];
A1(end,br_var)=1;
b1=[b;floor(br_value)];
% second LP problem with the added constraint that Xi >= ceil(Xi) , i=ind(1)
A2=[A ;zeros(1,c)];
A2(end,br_var)=-1;
b2=[b; -ceil(br_value)];
% solve the first LP problem
[x1,val1,status1,bound1]=rec(f,A1,b1,Aeq,beq,lb,ub,x0,val0,M,e,bound);
status=status1;
if status1 >0 & bound1<bound % if the solution was successfull and gives a better bound
   xx=x1;
   val=val1;
   bound=bound1;
   bb=bound1;
else
    xx=x0;
    val=val0;
    bb=bound;
end
    
% solve the second LP problem
% [x2,val2,status2,bound2]=rec(f,A2,b2,Aeq,beq,lb,ub,x0,val0,M,e,bound);
% if status2 >0 & bound2<bound % if the solution was successfull and gives a better bound
%     status=status2;
%     xx=x2;
%     val=val2;
%     bb=bound2;
% end
end