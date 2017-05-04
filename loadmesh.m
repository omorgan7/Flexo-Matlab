function [vertices,vertex_indices] = loadmesh(mesh_name)
fileID = fopen(mesh_name);
vertices = zeros(1,3);
vertex_indices = zeros(1,3);

fileline = 1;
while fileline ~= -1
    fileline = fgets(fileID);
    if fileline(1) == 'v'
        vertices = [vertices; cell2mat(textscan(fileline(3:end),'%f %f %f'))];
    end
    if fileline(1) == 'f'
        vertex_indices = [vertex_indices; cell2mat(textscan(fileline(3:end),'%d %d %d'))];
    end
end
vertices = vertices(2:end,:);
vertex_indices = vertex_indices(2:end,:);
fclose(fileID);
end