function rigidinterp(mesh1, mesh2)
[v1,vi] = loadmesh(mesh1);
[v2,~] = loadmesh(mesh2);

close all

x = v1(:,1);
v1(:,1) = v1(:,2);
v1(:,2) = -x;
v1 = v1(:,1:2);

x = v2(:,1);
v2(:,1) = v2(:,2);
v2(:,2) = -x;
v2 = v2(:,1:2);
figure(1);
axis([-2, 2, -4, 4]);
t = 1/(10);
for i = 0:t:1
    newV = (1-i)*v1 + i*v2;
    clf;
    axis([-3, 3, -4, 4]);
    hold on;
    for j = 1:length(vi)
        points = vi(j,:);
        DrawBetweenTwoPoints(newV(points(1),:),newV(points(2),:));
        DrawBetweenTwoPoints(newV(points(2),:),newV(points(3),:));
        DrawBetweenTwoPoints(newV(points(3),:),newV(points(1),:));
    end
    drawnow;
end
end
