function result = scat3(filename)
% Saves a 3D scatterplot of function
%
% Filename must contain three columns: X, Y, Z

D = importdata(filename);
output_filename = strcat(filename, '.pdf');

X = D(:,1);
Y = D(:,2);
Z = D(:,3);
%scatter3(X,Y,Z,1,Z);
scatter3(X,Y,Z,10);
axis vis3d

%fig = gcf;
%set(fig, 'PaperPositionMode', 'auto')
%print('-dpdf', '-r0', output_filename);

fig = gcf;
set(fig, 'PaperPositionMode', 'auto')
print('-dpdf', '-r0', output_filename);
