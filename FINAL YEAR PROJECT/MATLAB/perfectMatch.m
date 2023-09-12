function G = perfectMatch(adjMat1,nodes)
adjMatx = convertBipartite(adjMat1,nodes);
% converts the given adjacency matrix into an M x N matrix where M and N are the sizes of the partite
% sets of the input graph. In case of regular bipartite graphs M = N. 
% Check the function "convertBipartite" defined below, first.
% We use the logical array adjMatx to divide the input graph into two sets as below
[~,partiteSet1] = find(adjMatx);
[~,partiteSet2] = find(~adjMatx);
M = length(partiteSet1);
N = length(partiteSet2);
adjMat = zeros(M,N); %this M x N matrix will hold the transformed adjacency matrix
for i = 1:M
    for j = 1:N
        if (adjMat1(partiteSet1(i),partiteSet2(j)) == 1)
            adjMat(i,j) = 1;
        end
    end
end
% adjMat(i,j) = 1 implies ith vertex of partite set 1 is connected to jth vertex in partite set 2
matchR = zeros(1,N);
for u = 1:M
    %for each vertex u in partite set 1, create a boolean vector seen
    %which donates the vertices in partite set 2 as false i.e. not explored from u
    seen = false(1,N); 
    bpm(adjMat,u,seen,N); %this allots a vertex for u in the partite set 2
end
A = matchR; %A(a) = b implies partiteSet1(a) is connected to partiteSet2(b) in the matching
B = zeros(1,M);
C = zeros(1,N);
for x = 1:length(A)
    B(x) = partiteSet1(x); %vertices in partite set 1 which are included in the matching
    C(x) = partiteSet2(A(x)); %the corresponding vertices in the other partite set
end
G = graph(B,C);   
    function x = bpm(adjMat,u1,seen1,N)
        for v = 1:N %for each vertex in partite set 2
            if (adjMat(u1,v) == 1 && seen1(1,v) == false)
                seen1(1,v) = true; %if v is connected to u1 from partite set 1 then mark seen1 as true
            if(matchR(v) == 0 || bpm(adjMat,matchR(v),seen1,N) == true)
                %if vertex v is not yet alloted to any vertex in partite set 1 or
                %previously alloted vertex to v can be assigned some different vertex in set 2,then
                matchR(1,v) = u1; %allot v to u1
                x = true;
                return;
            end
            end
        end
        x = false;
        return;
    end
end


% define this separately, this function returns an array color of size 1 x (no. of nodes) of 0s and 1s
% color(i) = j implies vertex i is in the jth partite set where j = 1 or 0 
function A = convertBipartite(adjMat,nodes)
    % This is a DFS based approach to identify the partite sets of a bipartite set by 
    % coloring them with two colors 0 and 1
    color = zeros(1,nodes) - 1; % initialize by coloring with -1
    % start is vertex 1; 
    pos = 1; 
    colorGraph(adjMat,nodes,pos,1);
    A = color;
    % Below function checks if input graph is two - colorable (bipartite) with two colors 1 and 0 
    function x = colorGraph(adjMat,V,pos,c)
        if(color(pos) ~= -1 && color(pos) ~= c)
             x = false; %if color of vertex number pos is not c or -1 that means 
             % the previous vertex and pos has same color thus graph is not bipartite, return.
             return;
        end
            
            % color this vertex pos as c and all its neighbours as 1-c
            color(pos) = c;
            ans1 = true; %this holds the value true or false according to fact that input graph is
                        % two - colorable or not.
            for i=1:V %loop through all vertices
                if(adjMat(pos,i) == 1) % check if vertex i is a neighbor of pos
                    if(color(i) == -1) % if yes and i is uncolored
                        ans1 = ans1 && colorGraph(adjMat,V,i,1-c);
                        %recurse with the vertex i colored as 1-c
                    end
                    if(color(i) ~=-1 && color(i) ~= 1-c)
                        %if neighbor i of pos is colored but its color is not 1-c, then we know that
                        % input graph is not two - colorable i . e. bipartite
                            x = false;  
                            return;
                    end
                 end
                 if (ans1 == false) %if at any point we encounter a situation like this, return false
                      x = false;
                      return;
                 end
            end
       %else, return true value                       
       x = true;
       return;
    end
end     