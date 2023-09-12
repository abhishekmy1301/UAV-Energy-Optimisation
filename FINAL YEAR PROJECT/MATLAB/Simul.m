clc
clear
close all
rng(144); % Could be any random natural number
seed=rng;
N=1; %Number of nodes excluding the origin
nodes=rand(N,2); % Coordinates of nodes

%Initalizing coordinates
rng(seed);
xvals=nodes(:,1);
yvals=nodes(:,2);
origin=[0.5 0.5];
nodes=[origin;nodes;origin]; % adding origin coordinates to the front and the back of the array to complete the path
% now length of nodes is N+2





% Original Path
% figure('Name','Path in Given Order of Nodes (Random)','NumberTitle','off');
% plotnodes(nodes,origin);
% for i=1:N+1 % not N-1 or N because we are including the origin as well
%     p1=nodes(i,1:2);
%     p2=nodes(i+1,1:2);
%     dp=p2-p1;
%     quiver(p1(1),p1(2),dp(1),dp(2),0,'black'); %prev 4 lines is for creating arrow pointers
% end
% hold off


%Nearest Neighbour Heuristic
figure('Name','Nearest Neighbour Heuristic','NumberTitle','off');
plotnodes(nodes,origin);
traversed=nodes(1:N+1,:);
current=1; % try starting from different points
minpos=1;
path_nnh=0;
for i=1:N+1 % because N+2 -1 = N+1 path length
    p1=nodes(current,1:2);
    for k=2:N+1
        if(isTraversed(traversed,nodes(k,:))==1 )
            min=dist(p1,nodes(k,1:2));
            break;
        end
    end
    for j=1:N+1
        if(isTraversed(traversed,nodes(j,:))==1 )
            p2=nodes(j,1:2);
            d=dist(p1,p2);
            if(d<=min)
                min=d; %distance of nearest node
                minpos=j; %index of nearest node
                min_p=p2; %nearest node

            end
        end
    end

    current=minpos;

    dp=min_p-p1;
    traversed(minpos,:)=[0 0];
    quiver(p1(1),p1(2),dp(1),dp(2),0,'black');
    path_nnh=path_nnh+dist(p1,min_p);
end


p1=min_p;
p2=origin;
dp=p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'black'); % Adding the last arrow manually to complete the path
path_nnh=path_nnh+dist(p1,p2);
fprintf('Cost of NNH: %f\n',path_nnh);



% Greedy Heuristic Approach
figure('Name','Greedy Heuristic','NumberTitle','off');
nodes(end,:) = [];

plotnodes(nodes,origin);
distmat=zeros(N+1,N+1); % matrix of distances between two edges


for i=1:N+1 % matrix of distances between two nodes
    p1=nodes(i,:);
    for j=1:i
        p2=nodes(j,:);
        distmat(i,j)=dist(p1,p2);
    end
end



traversed=zeros(N+1,1);
path_gh=0;

for j=1:(N+1)^2
    x=min_nz(distmat);
    if(x(3)==0)
        break;
    end

    if(traversed(x(1))<2 && traversed(x(2))<2 && ~areIndirectlyConnected(distmat,x(1),x(2)) )




        p1=nodes(x(1),:);p2=nodes(x(2),:);

        line(p1,p2);
        path_gh=path_gh+dist(p1,p2);
        path_nnh=path_nnh+dist(p1,p2);
        traversed(x(1))=traversed(x(1))+1;
        traversed(x(2))=traversed(x(2))+1;

        distmat(x(1),x(2))=-1;
        distmat(x(2),x(1))=-1;

    else
        distmat(x(1),x(2))=0;
        distmat(x(2),x(1))=0;
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
        distmat(a,i)=-1;
        distmat(i,a)=-1;
    end
end
line(p1,p2);
path_gh=path_gh+dist(p1,p2);
path=getPath(distmat,1);
fprintf('Cost of GH: %f\n',path_gh);
hold off


%LOCAL SEARCH

% Random Swap
% figure('Name','Local Search','NumberTitle','off');
% path=getPath(distmat,1); % chose greedy heuristic
% plotnodes(nodes,origin);
%
% path;
%
% temp=path(5);
% path(5)=path(3);
% path(3)=temp;
%
%
% plotPath(nodes,path);
% cost_rs=getCost(nodes,path);
%
% fprintf('Cost after Random Swaps: %f\n',cost_rs);


% 2-Opt

figure('Name','Greedy Heuristic - 2-Opt Optimal','NumberTitle','off');
path=getPath(distmat,1); %  greedy heuristic
cost_2opt=getCost(nodes,path);
plotnodes(nodes,origin);
plotPath(nodes,path);
for i=1:length(path)-5
    for j=i+2:length(path)-2
        greenLine(nodes(path(i),:),nodes(path(i+1),:));
        greenLine(nodes(path(j),:),nodes(path(j+1),:));
        path=swapEdges(path,i,j);
        c=getCost(nodes,path);
        if(c<cost_2opt)
            cost_2opt=c;
        else

            path=swapEdges(path,i,j);
        end
        
        plotPath(nodes,path);
    end
end
fprintf('Cost After 2-Opt Optimization: %f\n',cost_2opt);
% 3-Opt

% K-Opt

% FUNCTIONS

function path=swapEdges(path,p1,p2)
path=path(1:end-1);


path(p1+1:p2)=flip(path(p1+1:p2));
path=[path;path(1)];


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

function path=getPath(distmat,p1)
start=p1;
path=p1;

for k=1:length(distmat(1,:))
    for i=1:length(distmat(1,:))

        if(distmat(p1,i)==-1)
            if(length(path)>1 && i==path(end-1))

            else
                p1=i;

                path=[path;i];
                if(i==start)
                    return;
                end

            end




        end
    end
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
plotnodes(nodes,[0.5 0.5]);
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



function plotnodes(nodes,origin)
scatter(nodes(:,1),nodes(:,2),10,'blue','filled');
hold on
scatter(origin(1),origin(2),250,'black','x');
text(origin(1),origin(2),'Depot','HorizontalAlignment','left','FontSize',10,'VerticalAlignment','top');
for i=2:(length(nodes(:,1)))
    text(nodes(i,1),nodes(i,2),num2str(i),'HorizontalAlignment','left','FontSize',10,'VerticalAlignment','top');

end
xlim([0 1]); ylim([0,1]);
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