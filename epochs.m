function result = epochs(filename)
% Plots phenotype fitness and virtual fitness against epoch number.

fid = fopen(filename);
D = textscan(fid, '%f%f%f%s%s%s%d%s%s%f');
fclose(fid);

epoch = D{1};
ph_accl = D{2};
v_accl = D{3};
edges = D{7};
fitness = D{10};

plot(epoch,ph_accl,epoch,v_accl)
legend('phenotype acclivity', 'virtual acclivity', 'Location', 'nw');
fig = gcf;
set(fig,'PaperUnits','inches');
set(fig,'PaperPosition',[0.5 0.5 8 2]);
print('-dpdf', strcat(filename, '.pdf'));
