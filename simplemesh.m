function simplemesh(mesh_name)
[v,vi] = loadmesh(mesh_name);
close all
figure(1)
hold on;
x = v(:,1);
v(:,1) = v(:,2);
v(:,2) = -x;
v = v(:,1:2);

for i = 1:length(vi)
    points = vi(i,:);
    DrawBetweenTwoPoints(v(points(1),1:2),v(points(2),1:2));
    DrawBetweenTwoPoints(v(points(2),1:2),v(points(3),1:2));
    DrawBetweenTwoPoints(v(points(3),1:2),v(points(1),1:2));
end

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
handleIndex = 1;
for i = 2:3
    if((norm(v(handleIndices(i),:)-newMove)^2)<= (norm(v(handleIndices(handleIndex),:)-newMove)^2))
        handleIndex = i;
    end
end

while(1)
    plot(newMove(1),newMove(2),'b.','MarkerSize',20);
    tic;
    w = 1000;

    bx = zeros(3+length(vi)*3, 1);
    by = bx;

    % A = sparse(3+length(vi)*3,length(v));
    A = zeros(3+length(vi)*3,length(v));

    for i = 1:length(vi)
        bx(3*(i-1)+1) = v(vi(i,2),1) - v(vi(i,1),1);
        bx(3*(i-1)+2) = v(vi(i,3),1) - v(vi(i,2),1);
        bx(3*(i-1)+3) = v(vi(i,1),1) - v(vi(i,3),1);
        by(3*(i-1)+1) = -v(vi(i,1),2) + v(vi(i,2),2);
        by(3*(i-1)+2) = -v(vi(i,2),2) + v(vi(i,3),2);
        by(3*(i-1)+3) = -v(vi(i,3),2) + v(vi(i,1),2);
        A(3*(i-1)+1,vi(i,1)) = -1;
        A(3*(i-1)+2,vi(i,2)) = -1;
        A(3*(i-1)+3,vi(i,3)) = -1;
        A(3*(i-1)+1,vi(i,2)) = 1;
        A(3*(i-1)+2,vi(i,3)) = 1;
        A(3*(i-1)+3,vi(i,1)) = 1;
    end

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

    newVx = ((A'*A)^-1)*(A'*bx);
    newVy = ((A'*A)^-1)*(A'*by);
    newV = [newVx,newVy];
    clf
    toc
    hold on
    for i = 1:length(vi)
        points = vi(i,:);
        DrawBetweenTwoPoints(newV(points(1),:),newV(points(2),:));
        DrawBetweenTwoPoints(newV(points(2),:),newV(points(3),:));
        DrawBetweenTwoPoints(newV(points(3),:),newV(points(1),:));
    end

    x(handleIndex) = newMove(1);
    y(handleIndex) = newMove(2);
    plot(x,y,'r.','MarkerSize',20);
    newMove = ginput(1);
end
end