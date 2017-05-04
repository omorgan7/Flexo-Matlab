function [N_B_H,nbh_i] = findlocalNBH(vertices, vertex_indices, vi, vj)
%[NBH , NBH_indices] = findlocalNBH(vertices, vertex_indices, POI)
%returns a cell array containing the neighbour hood of the point of
%interest (POI), as well as the indices of all the neighbours in the
%neighbourhood.
nbh_1 = mod(find(vertex_indices==vi),length(vertex_indices));
nbh_2 = mod(find(vertex_indices==vj),length(vertex_indices));

count = 1;

if(length(nbh_2)>length(nbh_1))
    iterator = nbh_2;
    smaller = nbh_1;
else
    iterator = nbh_1;
    smaller = nbh_2;
end

for i = 1:length(smaller)
    if(nnz(smaller(i) == iterator))
        foundIndices(count) = iterator(smaller(i) == iterator);
        count = count+1;
    end
end
foundIndices(foundIndices == 0) = length(vertex_indices);
if(length(foundIndices) == 1)
    N_B_H = cell(3,1);
    N_B_H{1} = vertices(vi,:);
    N_B_H{2} = vertices(vj,:);
    for i = 1:3
       if(vertex_indices(foundIndices,i) ~= vi && vertex_indices(foundIndices,i) ~= vj)
          N_B_H{3} = vertices(vertex_indices(foundIndices,i),:);
          nbh_i(1) = vi;
          nbh_i(2) = vj;
          nbh_i(3) = vertex_indices(foundIndices,i);
          break;
       end
    end
    return;
end

index_set = [vertex_indices(foundIndices(1),:);vertex_indices(foundIndices(2),:)];
nbh_i = unique(index_set);
mask_1 = nbh_i == vi;
mask_2 = nbh_i == vj;
vlr = nbh_i(~(mask_1+mask_2));
N_B_H = cell(4,1);
N_B_H{1} = vertices(vi,:);
N_B_H{2} = vertices(vj,:);
N_B_H{4} = vertices(vlr(1),:);
N_B_H{3} = vertices(vlr(2),:);
nbh_i(1) = vi;
nbh_i(2) = vj;
nbh_i(3) = vlr(2);
nbh_i(4) = vlr(1);


end

