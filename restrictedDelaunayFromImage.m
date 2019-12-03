function [F,V] = restrictedDelaunayFromImage(I,OPT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Given a surface, I want to downsample it, and compute its restricted
%   delaunay triangulation
%
%   When using an image, all we need is a function to check if a line
%   segment intersects the boundary
%
%   Note I originally used d = 0.2 for most structures, and d = 0.3
%   neocortex
%   
%   I will switch to d = 0.1 for most structures, and d = 0.2 for neocortex
%   AND cerebellum
%
%   Note, I should consider using the neighbours function in the Delaunay
%   Triangulation class.  This may speed things up
%
%   In restrictedDelaunayFromImage I make this into a nice function with an
%   image and options
%   OPT.dx      - pixel size (default 1,1,1) (in the order X Y Z, not row col slice)
%   OPT.isoval  - isosurface value (default max(I(:))/2)
%   OPT.dilate  - number of times to dilate (default 0)
%   OPT.nIdeal  - target number of vertices, usually we get around 0.7*nIdeal (default 1400)
%   OPT.minD    - minimum distance regardless of nIdeal (default 0)
%   OPT.display - draw pictures (default 0)
%   OPT.useLight- use a light in display, this is very slow on some machines (default 0) 
%   OPT.average - average from neighbouring points (default 1)
%   OPT.useDjikstra-compute shortest paths on the surface, rather than Euclidean distance.  Good for thin structures.  In general quite slow.
%   OPT.vStart- points to include to start, they won't be removed, defaults to empty
%   OPT.verbose - display progress
if nargin == 0
    error('must have at least one input argument I');
end
I = double(I);% be sure, or we'll convert every time interp3 is called
if nargin == 1
    OPT = struct;
end

% default options
dx = [1,1,1];
isoval = max(I(:))/2;
dilate = 0;
nIdeal = 1400;
minD = 0;
display = 0;
useLight = 0;
average = 1;
useDijkstra = 0;
vStart = [];
verbose = 0;

if isfield(OPT,'dx')
    dx = OPT.dx;
end
if isfield(OPT,'isoval')
    isoval = OPT.isoval;
end
if isfield(OPT,'dilate')
    dilate = OPT.dilate;
end
if isfield(OPT,'nIdeal')
    nIdeal = OPT.nIdeal;
end
if isfield(OPT,'minD')
    minD = OPT.minD;
end
if isfield(OPT,'display')
    display = OPT.display;
end
if isfield(OPT,'useLight')
    useLight = OPT.useLight;
end
if isfield(OPT,'average')
    average = OPT.average;
end
if isfield(OPT,'useDijkstra')
    useDijkstra = OPT.useDijkstra;
end
if useDijkstra
    addpath /cis/home/dtward/Functions/dijkstra
end
if isfield(OPT,'vStart')
    vStart = OPT.vStart;
end
if isfield(OPT,'verbose')
    verbose = OPT.verbose;
end

    
% dilate the image if necessary
nX = [size(I,2),size(I,1),size(I,3)];
x_ = (0:nX(1)-1)*dx(1);
y_ = (0:nX(2)-1)*dx(2);
z_ = (0:nX(3)-1)*dx(3);
if dilate
    for dilateLoop = 1 : dilate
        I = convn(I,ones(3,3,3),'same');
        I = double(I > 0);
    end
end

% create the isosurface using marching cubes
[tris,triloc] = isosurface(x_,y_,z_,I,isoval);
triloc = double(triloc);% seems to come out as single!
nvertex = size(triloc,1);


% find the distance between points by comparing surface area to desired
% number of points
SA = sum( sqrt( sum( (cross(triloc(tris(:,2),:) - triloc(tris(:,1),:), triloc(tris(:,3),:) - triloc(tris(:,1),:))/2).^2 , 2) ) );
% we never get up to nIdeal
% looks like we get to almost 0.7*nIdeal
% roughly dmin^2*nIdeal = SA
dmin = sqrt(SA/nIdeal);
% we should have a real minimum of 0.1, dno't make the triangles smaller
% than the voxels
dmin = max(dmin,minD);


% display it
if display
patch('faces',tris,'vertices',triloc,'facecolor','none','edgecolor','c','edgealpha',1,'linewidth',1)
% view(0,0)
axis image
xlabel x;
ylabel y;
zlabel z;
if useLight
light
end
end



% select a subset of vertices
% start with the first vertex
% then choose the closest vertex that is more than dmin away
chooseOldWay = 1;
if chooseOldWay
nmax = nIdeal*2;% allocate arrays
x = zeros(nmax,3);
x(1,:) = triloc(1,:);
count = 1;
ind = 1;


% this is just for one case
% when done, uncomment this dist, and delete everything below
try
  dist = nan(nmax,nvertex);
%   dist = nan(nmax*1000,nvertex*1000); % this will definitely give out of memory
%   disp('about to test error')
%   error('test error')
catch
% % dist = nan(5000,nvertex); % will reallocate...
% nvertex
% nmax

% feb 2, 2017
% I don't think the dist varaible needs to be more than 1 row
dist = nan(1,nvertex);
end
% dist = nan(round(nmax/4),nvertex);
while count < nmax
    if verbose
        if ~mod(count,10)
            disp(['Adding vertex ' num2str(count) ' of at most ' num2str(nmax) ]);
        end
    end
    % find the closest vertex more than dmin away
    if ~useDijkstra
      if size(dist,1)>1
        dist(count,:) = pdist2(x(count,:),triloc);
      else
        dist(1,:) = pdist2(x(count,:),triloc);
      end
    else
      tmp = dijkstra(triloc,tris,ind);
      if size(dist,1)>1
        dist(count,:) = tmp;
      else
        dist(1,:) = tmp;
      end
    end
    %%%%%
%     mindist = min(dist(1:count,:),[],1);
    % the above is silly, I should be able to save the previous
    if count == 1
      if size(dist,1)>1
        mindist = min(dist(1:count,:),[],1);
      else
        mindist = min(dist(1,:),[],1);
      end
    else
      if size(dist,1) > 1
        mindist = min([mindistOld;dist(count,:)],[],1);
      else
        mindist = min([mindistOld;dist(1,:)],[],1);
      end
    end
    mindistOld = mindist;
    % this speeds it up SOO SOO much.  Stops it from being quadratic in
    % time.
    %%%%%
    
    
    themin = min(mindist(mindist>dmin));
    if isempty(themin)
        break;
    end
    ind = find(  mindist ==  themin , 1, 'first');    
    if isempty(ind)
        break;
    end
%     ind = ind(1); % switched to 1, 'first'

    count = count + 1;
    x(count,:) = triloc(ind,:);
end
x = x(1:count,:);
n = size(x,1);
% daniel's note, september 25, 2015
% This part is rediculously slow
% I think the problem is that I keep checking ALL of triloc
% in fact, if a vertex is too close, it will always be too close
% I really should trim them as I go somehow
% but this won't work with my rectangular array for dist%
%
%
% the new way is NOT faster, and doesn't give the outputs I need
%
%
else
nmax = nIdeal*2;% allocate arrays
x = zeros(nmax,3);
x(1,:) = triloc(1,:);
count = 1;
ind = 1;
triloc_ = triloc; % we will take vertices away from here
dist = nan(1,nvertex);
while count < nmax
    if verbose
        if ~mod(count,10)
            disp(['Adding vertex ' num2str(count) ' of at most ' num2str(nmax) ]);
            disp([num2str(size(triloc_,1)) ' vertices left' ]);
        end
    end

    % compute the distance from everyone in x to everone in triloc_
    % start by computing the distance from the last added point
    distNew = pdist2(x(count,:),triloc_);
    % check for anyone who is too close
    notTooClose = distNew >= dmin;
    if sum(notTooClose) == 0 % nothing left
        break;
    end
    % anyone who is too close to this point is ALWAYS too close
    % this changing size makes things too slow
    tic;
    dist = dist(:,notTooClose);
    t = toc;
    if ~mod(count,10);disp(['Removing columns from dist, ' num2str(t) ' s']);end
    tic;
    triloc_ = triloc_(notTooClose,:);
    t = toc;
    if ~mod(count,10);disp(['Removing rows from triloc_, ' num2str(t) ' s']);end
    distNew = distNew(notTooClose);
    % now we get all the distances that are not too close
    if count > size(dist,1)
        tic;
        dist = [dist;nan(size(dist))];
        t = toc;
        if ~mod(count,10);disp(['Expanding dist, ' num2str(t) ' s']);end
    end
    dist(count,:) = distNew;
    % and we look for the closest
    tic;
    mindist = min(dist(1:count,:),[],1);
    themin = min(mindist);
    ind = find(themin == mindist,1,'first');
    t = toc;
    if ~mod(count,10);disp(['Finding closest vertex, ' num2str(t) ' s']);end
    % add
    count = count + 1;
    x(count,:) = triloc_(ind,:);
end
x = x(1:count,:);
n = size(x,1);

end


if display
hold on;
scatter3(x(1:n,1),x(1:n,2),x(1:n,3))
axis image;
end

% I should smooth them, take the average of all the points closest to you
if average && size(dist,1)>1 % we can't do this if we don't have a distance matrix
dist = dist(1:n,:);
closest = zeros(nvertex,1);
for i = 1 : nvertex
    inds = find(dist(:,i) == min(dist(:,i)));
    closest(i) = inds(1);
end
xSave = x;
for i = 1 : n
    x(i,:) = mean( triloc(closest == i,:), 1);
end
if display
scatter3(x(:,1),x(:,2),x(:,3),'b.')
axis image;
end
end


% add some points so that nothing is on the convex hull
hull = [1 1 1;1 1 -1;1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1]*1e6;
x = [x;hull];

% we take a 3D delaunay triangulation
% we only keep a face if its dual edge intersects the original surface
if ~isempty(vStart)
% we want to remove from x the ones which are too close to vStart
pairwiseD = pdist2(vStart,x); % this has one column for each x, and one row for each vStart
tooClose = sum(pairwiseD<dmin/2,1);% a bit smaller than dmin is okay
DT = delaunayTriangulation([vStart;x(~tooClose,:)]);
x = DT.Points; % DT will include only the unique points
n = size(x,1)-8;
else
DT = delaunayTriangulation(x);
end
[V,R] = voronoiDiagram(DT);
nmax = 0;
for i = 1 : size(R,1)
    nmax = max(nmax,length(R{i}));
end

% figure;
% patch('faces',R_,'vertices',V,'edgecolor','b','facecolor','r')



% here is how to find the voronoi edges for a given face
% remember that the points V are the circumcenters of each tetrahedra
% they have the same size as the tetrahedra and the same indexing
% every face defines two neigboring tetrahedra, and a dual edge connecting
% their circumcenters

% so, loop through each tetrahedra
% loop through each face
% find a tetrahedra with the same face
% find the edge connecting circumcenters
% check if it intersects the original surface
% if so, add it!
% figure;axis image;
T = DT.ConnectivityList;
Tsort = sort(T,2);
F = nan(size(T,1)*4,3);
FList = nan(size(T,1)*4,3);% maximum size
count = 0;
listCount = 0;
if verbose
    disp('Beginning to check tetrahedrons')
end
for tLoop = 1 : size( T , 1)
    if verbose
        if ~mod(tLoop,100);
            disp(['Checking tetrahedron ' num2str(tLoop) ' of ' num2str( size( DT.ConnectivityList , 1))]);
        end
    end
    
    tet = T(tLoop,:);
    c = V(tLoop+1,:);
    faces = tet([1 2 3;2 1 4;1 3 4;3 2 4]);
    for fLoop = 1 : 4
        % find the other tetrahedra with this face
        myFace = sort( faces(fLoop,:) ,2) ;
        
        % check if we've already looked at this face
        if ~isempty(FList) && ( sum( myFace(1) == FList(1:listCount,1) & myFace(2) == FList(1:listCount,2) & myFace(3) == FList(1:listCount,3) ) )
            continue
        end
        listCount = listCount + 1;
        FList(listCount,:) = myFace;

        % if this face is using the "hull"
        if sum( myFace > n )
            TF = 0;
            continue;
        end
        
        test = (Tsort(:,1) == myFace(1) & Tsort(:,2) == myFace(2) & Tsort(:,3) == myFace(3)) | ...
            (Tsort(:,1) == myFace(1) & Tsort(:,2) == myFace(2) & Tsort(:,4) == myFace(3)) | ...
            (Tsort(:,1) == myFace(1) & Tsort(:,3) == myFace(2) & Tsort(:,4) == myFace(3)) | ...
            (Tsort(:,2) == myFace(1) & Tsort(:,3) == myFace(2) & Tsort(:,4) == myFace(3));


        ind = find(test);
        ind = ind(ind ~= tLoop);
        if isempty(ind)
            continue;
        end
        c2 = V(ind+1,:);
        
        % now we check if the line segment c c2 intersects the original
        % surface
%         TF = lineSurfaceIntersect([c;c2],tris,triloc);
%         TF = 1;
%         if TF
%         patch('faces',[1 2],'vertices',[c;c2],'edgecolor','r','linewidth',2)
%         else
%             patch('faces',[1 2],'vertices',[c;c2],'edgecolor','k')
%         end

        % the slowest part of the code by far is interp3
        % the slowest part of that is permuting the image [2 1 3]! 
        % why should that permutation be done every time
%         val = interp3(x_,y_,z_,I,[c(1),c2(1)],[c(2),c2(2)],[c(3),c2(3)],'nearest',0);
        % I will use the same technique as the interp3 function
        if ~exist('Interpolator','var')
            [X_,Y_,Z_] = meshgrid(x_,y_,z_);
            Interpolator = griddedInterpolant(permute(X_,[2 1 3]), permute(Y_,[2 1 3]), permute(Z_,[2 1 3]), permute(I,[2 1 3]), 'linear','nearest');% interp, extrap
        end
        
        % test if there is an intersection of the edge with the boundary
        val = Interpolator([c(1),c2(1)],[c(2),c2(2)],[c(3),c2(3)]);
        val1 = val(1);
        val2 = val(2);
        if (val1>isoval && val2<isoval) || (val1<isoval && val2>isoval)
            TF = 1;
        else
            TF = 0;
        end
        
        % daniel want's to do more tests to be extra thorough (Nov 1, 2014)
%         nTests = 5;
%         val = Interpolator(linspace(c(1),c2(1),nTests),linspace(c(2),c2(2),nTests),linspace(c(3),c2(3),nTests));
%         tests = val>isoval;
%         % if we change sign at all, we are good.  This will get thin regions 
%         if sum(tests) == nTests || sum(tests)==0
%             TF = 0;
%         else
%             TF = 1;
%         end
%         % if we change sign an odd number of times we're good?
%         if ~mod(abs(diff(tests)),2)
%             TF = 0;
%         else
%             TF = 1;
%         end
        

        
        if TF
            count = count + 1;
            F(count,:) = faces(fLoop,:);% don't use myface because it is in the wrong order (well it's in the wrong order anyway...)            
            if display
            if isempty( findobj('tag','mysurface'))
                h = patch('faces',F(1:count,:),'vertices',x,'facecolor','b','edgecolor','k','facealpha',1,'edgealpha',1,'tag','mysurface');
            else
                set(h,'faces',F(1:count,:));
            end
            % plot the edge why not!
%             patch('faces',[1 2],'vertices',[c;c2],'edgecolor','r')
%             axis image
            drawnow;
            end
        end
        
        
    end
end
F = F(1:count,:);
if verbose
    disp('Finished checking tetrahedrons')
end

% are any x's not being used? I should get rid of them here
if verbose
    disp('Beginning to check for unused vertices')
end
if ~(length( unique(F) ) == size(x,1)-8)
    warning('Some vertices not used, removing')
    % need to get rid of them!    
    check = 1;
    while 1
        all = unique(F);
        if ~sum( all == check ) % if vertex i is not here
            x = x([1:check-1,check+1:end],:);
            F(F>=check) = F(F>=check)-1;
        else
            check = check+1;
        end
        if check > length(x)
            break;
        end
    end
    n = size(x,1);
else
    disp('All vertices used')
    x = x(1:n,:);
end
if verbose
    disp('Finished checking for unused vertices')
end

% if display
% figure
% h = patch('faces',F,'vertices',x,'facecolor','b','edgecolor','k','facealpha',1,'edgealpha',0.25,'backfacelighting','lit');
% axis image
% if useLight
% light
% end
% end

% try to orient the faces
if verbose
    disp('Beginning to orient faces')
end
[F] = orientSurface(x,F);
vol = sum( sum( x(F(:,1),:) .* cross(x(F(:,2),:), x(F(:,3),:) ) ) )/6;
if vol < 0
    F(:,[2 3]) = F(:,[3 2]);
end
if verbose
    disp('Finished orienting faces')
end
% the output vertex array
V = x;
if verbose
    disp('Finished restrictedDelaunayFromImage, exiting normally.')
end


