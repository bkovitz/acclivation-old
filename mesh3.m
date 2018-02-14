function result = mesh3(filename)
% Saves a 3D mesh plot of a fitness function

D = importdata(filename);
output_filename = strcat(filename, '.pdf');

X = D(:,1);
Y = D(:,2);
Z = D(:,3);
tri = delaunay(X,Y);
trimesh(tri,X,Y,Z);
axis vis3d;

fig = gcf;
set(fig, 'PaperPositionMode', 'auto')
print('-dpdf', '-r0', output_filename);

% Relevant Matlab functions:
%   view                   Gets/sets current azimuth and elevation
%   print -depsc filename  Saves current plot to filename
