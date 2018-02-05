D = importdata('w.data')
[X,Y] = meshgrid(-1:0.01:1, -1:0.01:1)
Z = reshape(D(:,3),[201,201])
surf(X,Y,Z)
