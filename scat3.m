function result = scat3(filename)
% See a 3D scatterplot of function
%
% Filename must contain three columns: X, Y, Z

D = importdata(filename);
X = D(:,1);
Y = D(:,2);
Z = D(:,3);
scatter3(X,Y,Z,1,Z);
axis vis3d
