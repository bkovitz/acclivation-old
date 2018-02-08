function result = seefit(filename, step)
% See a 3D mesh plot of a fitness function

D = importdata(filename);

% Y = D(:,2);
% lb = Y(1);
% step = round((Y(2) - Y(1)) * 1000000) / 1000000;
% ub = Y(end);
% len=sqrt(length(Y));
% [X,Y] = meshgrid(lb:step:ub, lb:step:ub);
% Z = reshape(D(:,3),[len,len]);
% %size(X),size(Y),size(Z)
% result = mesh(X,Y,Z);

X = D(:,1);
Y = D(:,2);
Z = D(:,3);
tri = delaunay(X,Y);
trimesh(tri,X,Y,Z)
axis vis3d

% Relevant Matlab functions:
%   view                   Gets/sets current azimuth and elevation
%   print -depsc filename  Saves current plot to filename
