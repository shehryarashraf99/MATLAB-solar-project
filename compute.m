function output = compute(X,Y,Z)
%% check if Tracker half Z is shading Tracker half X. with Shadow Coordinates of Tracker half Z in Y

% get plane equation of tracker half x, use three points.
p1=X(:,1);
p2=X(:,2);
p3=X(:,length(X));

%% create multiple copies of normal vector
n=cross(p1-p2,p1-p3);
n=repmat(n,1,length(Y));
p1=repmat(p1,1,length(Y));
%% calculate intersections for all the shading lines of tracker half Z on plane denoted by tracker half X corners

u=Z-Y;
w=Y-p1;

D=dot(n,u);
N=-dot(n,w);



sI = N./D;
sI=repmat(sI,3,1);
I = Y+ sI.*u;

% get range of values in all 3 dimensions to check if intersection lies in
% tracker half X
x_max=max(X(1,:));
x_min=min(X(1,:));
y_max=max(X(2,:));
y_min=min(X(2,:));
z_max=max(X(3,:));
z_min=min(X(3,:));

shaded=zeros(1,length(X));


%% check if coordinates of intersection lie in tracker half X

 for i=1:length(I)
      if sI(i) < 0 || sI(i) > 1
          intersect=0;
      else
          intersect=1;
      end
if I(1,i)>=x_min && I(1,i)<=x_max && I(2,i)<=y_max && I(2,i)>=y_min && I(3,i)>=z_min && I(3,i)<=z_max   && intersect==1
    shaded(i)=1;
else
    shaded(i)=0;
end
 end
 
output=any(shaded); % if any panel in tracker half z shades any part of tracker half x, tracker half z shaded tracker half x.




end