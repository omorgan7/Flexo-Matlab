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

A_gamma = cell(length(vi));
P_inv_saved = cell(length(vi));
A = zeros(4*length(vi),2*(length(v1)-1));
for i = 1:length(vi)
    p = [v1(vi(i,1),:);v1(vi(i,2),:);v1(vi(i,3),:)];
    q = [v2(vi(i,1),:);v2(vi(i,2),:);v2(vi(i,3),:)];
    [P_inv, a_out] = trerp(p,q);
    A_gamma{i} = a_out;
    P_inv_saved{i} = P_inv;
    for k = 1:4
        for j = 1:3
            if(vi(i,j)>1)
                A(4*i-4+k,2*vi(i,j)-3:2*vi(i,j)-2) = P_inv(k,2*j-1:2*j);
            end
        end
    end
end
A_inv = (A'*A)^(-1)*A';
newV = zeros(2*length(v1),1);
figure(1);
hold on;
for i = 0:t:1
    newV(1) = (1-i)*v1(1,1) + i*v2(1,1);
    newV(2) = (1-i)*v1(1,2) + i*v2(1,2);
    b = zeros(4*length(vi),1);
    for j = 1:length(vi)
        A_out = A_gamma_calc(A_gamma{j},i)';
        b(4*j-3:4*j) = A_out(:);
        if(nnz(vi(j,:) == 1))
            idx = find(vi(j,:) == 1);
            P_inv = P_inv_saved{j};
            b(4*j-3:4*j) = b(4*j-3:4*j) - P_inv(:,2*idx-1)*newV(1) - P_inv(:,2*idx)*newV(2);
        end
    end
    newV(3:end) = A_inv * b;
    newV = reshape(newV,[2,length(v1)]);
    newV = newV';
    
    clf;
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
