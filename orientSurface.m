function F = orientSurface(V,F);
%   orient surface
%
%   we will orient a surface by checking if points are inside or outside
%


% keyboard

draw = 0;
if nargin == 0
    x = -2:2;
    [X,Y,Z] = meshgrid(x);
    I = double( (X.^2 + Y.^2 + Z.^2) <= 0.5^2 );
    [F,V] = isosurface(I,0.5);
    for i = 1 : size(F,1)
        if rand > 0.5
            F(i,[1 2 3]) = F(i,[2 1 3]);
        end
    end
    draw = 1;
end


if draw
patch('faces',F,'vertices',V,'facecolor','c','edgecolor','k','facealpha',1,'edgealpha',1)
axis image
hold on;
end

% how far to step
epsilon = 1e-9;

% consider a ray in direction u (any direction should be okay)
u = randn(1,3);
u = u/norm(u);


% calculate all the face centers and normals
c = ( V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:) )/3;
n = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:))/2;
n = bsxfun(@rdivide,n,sqrt( sum(n.^2,2) ));
if draw
    scatter3(c(:,1),c(:,2),c(:,3))
    quiver3(c(:,1),c(:,2),c(:,3),n(:,1),n(:,2),n(:,3))
end

% at each point c + epsilon n, we check if the point is inside or outside
p = c + epsilon*n;

% we have to count how many times a ray intersects the surface
count = zeros(size(F,1),1);


% we do this by looping over each triangle
for i = 1 : size(F,1)
    % define points in our triangle
    q = V(F(i,1),:);
    v = V(F(i,2),:) - V(F(i,1),:);
    w = V(F(i,3),:) - V(F(i,1),:);
    
    
    % we have to solve the equation
    % p + au = q + bv + cw
    % p-q = -au + bv + cw
    pmq = bsxfun(@minus,p,q);
    nUVW = [-u;v;w];
    % so we have abc*nUVW = pmq
    abc = pmq/(nUVW);
    
    % now we put it in barycentric corrdinates
    triCoord = [1-(abc(:,2)+abc(:,3)),abc(:,2),abc(:,3)];
    intersect = ( (sum( abs(triCoord - 0.5) <= 0.5 , 2) == 3)  & (abc(:,1)>=0));
    if draw
        a = abc(:,1);
        if sum(intersect)
            patch('vertices',V,'faces',F(i,:),'facecolor','r','facealpha',1)
        else
            patch('vertices',V,'faces',F(i,:),'facecolor','b','facealpha',1)
        end
        scatter3(p(intersect,1) + a(intersect)*u(1),p(intersect,2) + a(intersect)*u(2),p(intersect,3) + a(intersect)*u(3),'m.');
        drawnow
    end
    count = count + intersect;
end

% keyboard
% now if the count is odd, we are inside, and we should flip the face
% for i = 1 : size(F,1)
%     if mod(count(i),2) % if odd
%         F(i,[1 2 3]) = F(i,[2 1 3]);
%         disp(['Flipping face ' num2str(i)]);
%     end
% end

ind = mod(count,2) == 1; % odd (i.e. p is inside)
F(ind,[1 2 3]) = F(ind,[2 1 3]);


% calculate all the face centers and normals
if draw
    n = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:))/2;
    n = bsxfun(@rdivide,n,sqrt( sum(n.^2,2) ));
    quiver3(c(:,1),c(:,2),c(:,3),n(:,1),n(:,2),n(:,3),'r:')
end

