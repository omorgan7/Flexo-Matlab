function mesh(mesh_name)
%%%%%%%%%%%%%%%%%%%%%%%
%runs the mesh deformation 
%method implemented in Igarashi and Igarashi's 2009 paper.
%input: a string containing the location of the mesh.

%load mesh
[v,vertexIndices] = loadmesh(mesh_name);

%draw mesh
close all
figure(1)
hold on;
x = v(:,1);
v(:,1) = v(:,2);
v(:,2) = -x;
v = v(:,1:2);
for i = 1:length(vertexIndices)
    points = vertexIndices(i,:);
    DrawBetweenTwoPoints(v(points(1),1:2),v(points(2),1:2));
    DrawBetweenTwoPoints(v(points(2),1:2),v(points(3),1:2));
    DrawBetweenTwoPoints(v(points(3),1:2),v(points(1),1:2));
end

%get input points
[x,y] = ginput(3);
plot(x,y,'r.','MarkerSize',20);
handleIndices = zeros(3,1);
for j = 1:3
    for i = 1:length(v)
        if((norm(v(i,:) - [x(j), y(j)])^2)<0.01)
            handleIndices(j) = i;
            break;
        end
    end
end
newMove = ginput(1);

%work out which verticies these correspond to, simple euclidean norm.
handleIndex = 1;
for i = 2:3
    if((norm(v(handleIndices(i),:)-newMove)^2)<= (norm(v(handleIndices(handleIndex),:)-newMove)^2))
        handleIndex = i;
    end
end
while(1)
    plot(newMove(1),newMove(2),'b.','MarkerSize',20);
    
    %build A matrix, b vector.
    w = 1000;
    A = zeros(numel(vertexIndices)*2 + 6,numel(v));
    b = zeros(numel(vertexIndices)*2 + 6, 1);
    edgeNBH = zeros(numel(vertexIndices),5);
    G_cell = cell(numel(vertexIndices),1);
    for i = 1:3
        if(i == handleIndex)
            b(end-6+2*(i-1)+1) = w*newMove(1);
            b(end-6+2*i) = w*newMove(2);
        else
            b(end-6+2*(i-1)+1) = w*v(handleIndices(i),1);
            b(end-6+2*i) = w*v(handleIndices(i),2);
        end
        A(end-6+2*(i-1)+1,handleIndices(i)*2-1) = w;
        A(end-6+2*i,handleIndices(i)*2) = w;
    end
    %%% Algorithm 1
    for i = 1:length(vertexIndices)
        for j = 1:3
            %for each edge, get the edge neighbourhood.
            [NBH,NBH_indices]=findlocalNBH(v,vertexIndices,vertexIndices(i,j),vertexIndices(i,mod(j,3)+1));
            numIndices = length(NBH_indices);
            edgeNBH(3*(i-1)+j,1:numIndices) = NBH_indices;
            edgeNBH(3*(i-1)+j,5) = numIndices;
            edgeArr = zeros(2,numIndices*2);
            G = zeros(numIndices*2,4);
            for k = 1:2:2*numIndices
                G(k,1:2) = NBH{ceil(k/2)};
                G(k+1,1:2) = fliplr(NBH{ceil(k/2)});
                G(k+1,2) = -G(k+1,2);
                G(k:k+1,3:4) = eye(2);
            end
            
            H = (G/(G'*G))';
            H = H(1:2,:);
            edgeArr(:,1:4) = [-1,0,1,0;0,-1,0,1];
            G_cell{3*(i-1)+j} =H;
            edge = NBH{2} - NBH{1};
            edgeMat = [edge(1), edge(2);edge(2),-edge(1)];
            
            h = edgeArr - edgeMat*H;
            for k = 1:numIndices
                A(6*(i-1)+2*j-1,2*NBH_indices(k)-1:2*NBH_indices(k)) = h(1,2*k-1:2*k);
                A(6*(i-1)+2*j,2*NBH_indices(k)-1:2*NBH_indices(k)) = h(2,2*k-1:2*k);
            end
        end
    end
    newV = ((A'*b)'/(A'*A))';
    newV = reshape(newV',[2,length(v)])';
    %%% end of algorithm 1

    %algorithm 2
    %build the new A, and bx,by
    A = zeros(numel(vertexIndices)+3,length(v));
    bx = zeros(numel(vertexIndices)+3,1);
    by = zeros(numel(vertexIndices)+3,1);
    for i = 1:3
        if(i == handleIndex)
            bx(end-3+i) = w*newMove(1);
            by(end-3+i) = w*newMove(2);
        else
            bx(end-3+i) = w*v(handleIndices(i),1);
            by(end-3+i) = w*v(handleIndices(i),2);
        end
        A(end-3+i,handleIndices(i)) = w;    
    end
    
    for i = 1:numel(vertexIndices)
        numIndices = edgeNBH(i,5);
        vi_nbh = edgeNBH(i,1:numIndices);
        v_nbh = zeros(numIndices,2);
        for j =1:numIndices
            v_nbh(j,:) = newV(vi_nbh(j),:);
        end
        v_nbh = v_nbh';
        cs = G_cell{i}*v_nbh(:);
        T = [cs';fliplr(cs')];
        T(2,1) = -T(2,1);
        T = T/norm(cs);
        edge = v(vi_nbh(2),:)-v(vi_nbh(1),:);
        Te = T*edge';
        bx(i) = Te(1);
        by(i) = Te(2);
        A(i,vi_nbh(1)) = -1;
        A(i,vi_nbh(2)) = 1;
    end
    
    newVx = ((A'*A)^-1)*(A'*bx);
    newVy = ((A'*A)^-1)*(A'*by);
    newV = [newVx,newVy];
    %end of algorithm 2
    
    %draw:
    clf
    hold on
    for i = 1:length(vertexIndices)
        points = vertexIndices(i,:);
        DrawBetweenTwoPoints(newV(points(1),:),newV(points(2),:));
        DrawBetweenTwoPoints(newV(points(2),:),newV(points(3),:));
        DrawBetweenTwoPoints(newV(points(3),:),newV(points(1),:));
    end
    
    x(handleIndex) = newMove(1);
    y(handleIndex) = newMove(2);
    plot(x,y,'r.','MarkerSize',20);
    newMove = ginput(1);
end
