


function path = eulerian_path(dist_matrix)
    % Create adjacency matrix from distance matrix
    adj_matrix = dist_matrix;
    
    % Find any vertex with odd degree
    odd_vertices = find(mod(sum(adj_matrix), 2));
    start_vertex = odd_vertices(1);
    if isempty(start_vertex)
        start_vertex = 1; % If no odd degree vertices, start anywhere
    end
    
    % Initialize stack and path
    stack = start_vertex;
    path = [];
    
    while ~isempty(stack)
        % Take the last vertex on the stack
        current_vertex = stack(end);
        
        % Find all unvisited edges from the current vertex
        unvisited_edges = find(adj_matrix(current_vertex, :) == 1);
        visited_edges = find(adj_matrix(current_vertex, :) == 2);
        
        if ~isempty(unvisited_edges)
            % Visit the first unvisited edge
            next_vertex = unvisited_edges(1);
            
            % Mark the edge as visited
            adj_matrix(current_vertex, next_vertex) = 2;
            adj_matrix(next_vertex, current_vertex) = 2;
            
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
