function coordinates = compute_intersections(normal,point_in_plane,P1,P2)
%% doing transpose of all inputs
normal=transpose(normal);
point_in_plane=transpose(point_in_plane);
P1=transpose(P1);
%% making copies of points
normal=repmat(normal,1,length(P2));
point_in_plane=repmat(point_in_plane,1,length(P2));
P1=repmat(P1,1,length(P2));
%% computing dot products and all the intersection points
 u=P1-P2;
w=P2-point_in_plane;
D=dot(normal,u);
N=-dot(normal,w);

sI = N./D;
sI=repmat(sI,3,1);
coordinates = P2+ sI.*u;



end

