function [total_dist,total_dist_v_sqe,total_dist_v_cb,total_dist_v_cos,total_dist_v_cor,t] = Final(n_values)
close all
%Creating Values
tic
rng('default')
rng(1)
n_values=10;
[x,y] = CreateVals(n_values,10);

x_o = 5;
y_o = 5;

%Only here when calculating time complexity
total_dist_v_sqe = 0;
total_dist_v_cb= 0;
total_dist_v_cos = 0;
total_dist_v_cor = 0;

k = 2;
cost = 2;

idx = {};

for q = 1:5
    %Find Centers
    rng shuffle
    c = CreateCenters(x,y,k);
    
    %Clustering - First pass
    [m,idx{q}] = Clustering1(x,y,c);
    
    %Clustering - Second pass
    idx{q} = Clustering2(idx{q},x_o,y_o,x,y);
    
    %Clustering - Third pass
    idx{q} = Clustering3(idx{q},x,y);
    
    %Path optimisation for all clusters
    nodes = {};
    D = {};
    path = {};
    dist = [];
    for i = 1:length(c)
        nodes{i} = [x_o y_o;x(idx{q}==i).' y(idx{q}==i).'];
        D{i} = MatrixDist(nodes{i});
        path{i} = NNH(1,D{i});
        path{i} = twoopt(path{i},D{i});
        dist(i) = distTSP(path{i},D{i});
    end
    
    total_dist(q) = sum(dist);
    total_dist_s(q) = sumsqr(dist);

end

[total_dist,iter] = min(total_dist);
idx = idx{iter};

%Path optimisation for all clusters
nodes_t = {};
D_f = {};
path_f = {};
dist_f = [];
for i = 1:length(c)
    nodes_t{i} = [x_o y_o;x(idx==i).' y(idx==i).'];
    D_f{i} = MatrixDist(nodes_t{i});
    path_f{i} = NNH(1,D_f{i});
    path_f{i} = twoopt(path_f{i},D_f{i});
    dist_f(i) = distTSP(path_f{i},D_f{i});
end

total_dist = sum(dist_f);
total_dist_s = sumsqr(dist_f);

%Plotting
C = {'b','r','g','y','m','c','k'};
C2 = {'b*','r*','g*','y*','m*','c*','k*'};
figure(1)
for i = 1:length(c)
    PrintSol(path{i},nodes{i},C2{i},C{i});
    hold on
end
title(['Clustering (d = ',num2str(total_dist),')'])



%{
MUTED FOR THE TIME COMPLEXITY
%VALIDATION - sqeuclidean
idx_v_sqe = kmeans([x.' y.'],k,'Distance','sqeuclidean');
idx_v_sqe = idx_v_sqe.';

nodes_v_sqe = {};
D_v_sqe = {};
path_v_sqe = {};
dist_v_sqe = [];
for i = 1:length(c)
    nodes_v_sqe{i} = [x_o y_o;x(idx_v_sqe==i).' y(idx_v_sqe==i).'];
    D_v_sqe{i} = MatrixDist(nodes_v_sqe{i});
    path_v_sqe{i} = NNH(1,D_v_sqe{i});
    path_v_sqe{i} = twoopt(path_v_sqe{i},D_v_sqe{i});
    dist_v_sqe(i) = distTSP(path_v_sqe{i},D_v_sqe{i});
end



total_dist_v_sqe = sum(dist_v_sqe)
total_dist_v_s_sqe = sumsqr(dist_v_sqe);

figure(2)
for i = 1:length(c)
    PrintSol(path_v_sqe{i},nodes_v_sqe{i},C2{i},C{i});
    hold on
end
title(['Clustering squared eucledian distance validation (d = ',num2str(total_dist_v_sqe),')'])

%VALIDATION - cityblock
idx_v_cb = kmeans([x.' y.'],k,'Distance','cityblock');
idx_v_cb = idx_v_cb.';

nodes_v_cb = {};
D_v_cb = {};
path_v_cb = {};
dist_v_cb = [];
for i = 1:length(c)
    nodes_v_cb{i} = [x_o y_o;x(idx_v_cb==i).' y(idx_v_cb==i).'];
    D_v_cb{i} = MatrixDist(nodes_v_cb{i});
    path_v_cb{i} = NNH(1,D_v_cb{i});
    path_v_cb{i} = twoopt(path_v_cb{i},D_v_cb{i});
    dist_v_cb(i) = distTSP(path_v_cb{i},D_v_cb{i});
end



total_dist_v_cb = sum(dist_v_cb)
total_dist_v_s_cb = sumsqr(dist_v_cb);

figure(3)
for i = 1:length(c)
    PrintSol(path_v_cb{i},nodes_v_cb{i},C2{i},C{i});
    hold on
end
title(['Clustering cityblock distance validation (d = ',num2str(total_dist_v_cb),')'])

%VALIDATION - cosine
idx_v_cos = kmeans([x.' y.'],k,'Distance','cosine');
idx_v_cos = idx_v_cos.';

nodes_v_cos = {};
D_v_cos = {};
path_v_cos = {};
dist_v_cos = [];
for i = 1:length(c)
    nodes_v_cos{i} = [x_o y_o;x(idx_v_cos==i).' y(idx_v_cos==i).'];
    D_v_cos{i} = MatrixDist(nodes_v_cos{i});
    path_v_cos{i} = NNH(1,D_v_cos{i});
    path_v_cos{i} = twoopt(path_v_cos{i},D_v_cos{i});
    dist_v_cos(i) = distTSP(path_v_cos{i},D_v_cos{i});
end



total_dist_v_cos = sum(dist_v_cos)
total_dist_v_s_cos = sumsqr(dist_v_cos);

figure(4)
for i = 1:length(c)
    PrintSol(path_v_cos{i},nodes_v_cos{i},C2{i},C{i});
    hold on
end
title(['Clustering cosine distance validation (d = ',num2str(total_dist_v_cos),')'])

%VALIDATION - correlation
idx_v_cor = kmeans([x.' y.'],k,'Distance','correlation');
idx_v_cor = idx_v_cor.';

nodes_v_cor = {};
D_v_cor = {};
path_v_cor = {};
dist_v_cor = [];
for i = 1:length(c)
    nodes_v_cor{i} = [x_o y_o;x(idx_v_cor==i).' y(idx_v_cor==i).'];
    D_v_cor{i} = MatrixDist(nodes_v_cor{i});
    path_v_cor{i} = NNH(1,D_v_cor{i});
    path_v_cor{i} = twoopt(path_v_cor{i},D_v_cor{i});
    dist_v_cor(i) = distTSP(path_v_cor{i},D_v_cor{i});
end



total_dist_v_cor = sum(dist_v_cor)
total_dist_v_s_cor = sumsqr(dist_v_cor);

figure(5)
for i = 1:length(c)
    PrintSol(path_v_cor{i},nodes_v_cor{i},C2{i},C{i});
    hold on
end
title(['Clustering correlation distance validation (d = ',num2str(total_dist_v_cor),')'])
%}

t=toc



%PROGRAM FUNCTIONS
%Fuction 1 - Creating Values
function [x,y] = CreateVals(n_values,range)
    x = rand(1,n_values)*range;
    y = rand(1,n_values)*range;
end

%Function 2 - Find cluster centers
function [c] = CreateCenters(x,y,k)
%Select a random point as the first cluster
    init_c = randi(length(x));
    c(:,1) = [x(init_c);y(init_c)];

    %For each point, find it's squared distance to c1 and create a
    %probability distribution with those numbers. Find second center
    %randomly with that given distribution
    s_distance = (x-c(1,1)).^2+(y-c(2,1)).^2;
    prob = s_distance./sum(s_distance);
    index_c2 = find(rand<cumsum(prob),1,'first');
    c(:,2) = [x(index_c2);y(index_c2)];

    %To find the rest of the centers, do the same thing but only use the
    %minimum distance from each point to a cluster center
    for i=3:k
        %Find the distance to each cluster center
        for j = 1:length(x)
            dist_2_center = zeros(1,length(c));
            for m = 1:length(c)
                dist_2_center(m) = (x(j)-c(1,m))^2+(y(j)-c(2,m))^2;
            end
            final_s_dist(j) = min(dist_2_center);
        end
        prob = final_s_dist./sum(final_s_dist);
        index = find(rand<cumsum(prob),1,'first');
        c(:,i) = [x(index);y(index)];
    end
end

%Function 3 - Clustering (FIRST UPDATE)
function[m,idx] = Clustering1(x,y,c)
    idx = [];
    k = length(c);
    cond = 0;
    m = c;

    flag = 0;

    while (cond == 0 && flag < 10000)
        flag = flag + 1;
        m_old = m;
        %For every point, find the distance to each cluster center and
        %assign it to the cluster with minimum distance

        for i = 1:length(x)
            point = [x(i) y(i)];
            for j = 1:k
                d(j) = Distance2Points(point,c(:,j));
            end
            [m,cluster] = min(d);
            idx = [idx cluster];
        end
        
        %Find new cluster centers
        for i = 1:k
            m(1,i) = mean(x(idx==i));
            m(2,i) = mean(y(idx==i));
        end

        %Check if cluster centers keep changing
        if (m == m_old)
            idx = [];
        else
            cond = 1;
        end

    end
end

%Function 4 - Clustering (SECOND UPDATE)
function [idx] = Clustering2(idx,x_o,y_o,x,y)
    kx = {};
    ky = {};
    k = max(idx);
    
        for i = 1:length(x)
            for p = 1:k
                kx{p} = x(idx==p);
                ky{p} = y(idx==p);
            end
    
            for j = 1:k
                nodes = [x_o y_o;kx{j}.' ky{j}.'];
                D = MatrixDist(nodes);
                path = NNH(1,D);
                path = twoopt(path,D);
                dist(j) = distTSP(path,D)^cost;
            end
            old_sum_of_path = sum(dist);
            idx_new = idx;
            for j = 1:k
                if j ~= idx(i)
                    idx_new(j) = j;
                    for n = 1:k
                        kx{n} = x(idx_new==n);
                        ky{n} = y(idx_new==n);
                        nodes = [x_o y_o;kx{n}.' ky{n}.'];
                        D = MatrixDist(nodes);
                        path = NNH(1,D);
                        path = twoopt(path,D);
                        new_dist(n) = distTSP(path,D)^cost;
                    end
                    new_sum_of_path(j) = sum(new_dist);
                else
                    new_sum_of_path(j) = old_sum_of_path;
                end
            end
            [~,cluster] = min(new_sum_of_path);
            idx(i) = cluster;
        end
end

%Clustering (THIRD UPDATE)
function [idx] = Clustering3(idx,x,y)
    k = max(idx);
    %Find cluster centers
    for i = 1:k
            m(1,i) = mean(x(idx==i));
            m(2,i) = mean(y(idx==i));
    end
    cond = 0;
    %Perform first update again
    [m,idx] = Clustering1(x,y,m);
end

%AUXILIARY FUNCITONS
%Funciton 1 - Finding the distance between 2 points
function [d] = Distance2Points(a,b)
    d = sqrt((a(1)-b(1)).^2 + (a(2)-b(2)).^2);
end
end