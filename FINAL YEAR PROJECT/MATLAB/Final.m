clc;clear;close all;
tic
rng(144) % Could be any random natural number
seed=rng;
N=19; %Number of nodes excluding the origin
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
cost_rnd2=getCost(nodes,path_rnd);
while(~(cost_rnd2==cost_rnd1))
    cost_rnd2=getCost(nodes,path_rnd);
    path_rnd=two_Opt(nodes,path_rnd);
    cost_rnd1=getCost(nodes,path_rnd);
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
path_nnh=two_Opt(nodes,path_nnh);
cost_nnh2=getCost(nodes,path_nnh);
while(~(cost_nnh2==cost_nnh1))
    cost_nnh2=getCost(nodes,path_nnh);
    path_nnh=two_Opt(nodes,path_nnh);
    cost_nnh1=getCost(nodes,path_nnh);
end
 plotNodes(nodes);

plotPathDirection(nodes,path_nnh);


fprintf('Cost After 2-Opt Optimization: %f\n',getCost(nodes,path_nnh));
plotPathDirection(nodes,path_nnh);





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


 path_gh_2opt=two_Opt(nodes,path_gh_2opt); path_gh_2opt=two_Opt(nodes,path_gh_2opt);



plotNodes(nodes,origin);

plotPathDirection(nodes,path_gh_2opt);
 fprintf('Cost After 2-Opt Optimization: %f\n',getCost(nodes,path_gh_2opt));





toc
% OPTIMISATION FUNCTIONS

function min_path=two_Opt(nodes,path)
cost_2opt=getCost(nodes,path);
min_path=path;

for i=1:length(path)-3 % first edge
    for j=i+2:length(path)-1    % second edge will always be after the first edge
        
%         g=dist(nodes(path(j),:),nodes(path(j+1),:));
        path=swapEdges(path,i,j);
%         h=dist(nodes(path(i),:),nodes(path(i+1),:));
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

