function result = epochs(filename)
% Plots phenotype fitness and virtual fitness against epoch number.
%
% Expected data format corresponds to the (genotype-data) function in
% Clojure.

output_filename = strcat(filename, '.pdf');

fid = fopen(filename);
D = textscan(fid, '%d%d%f%f%s%s%s%d%s%s%f');
fclose(fid);

epoch = D{1};
generation = D{2};
ph_accl = D{3};
v_accl = D{4};
edges = D{8};
fitness = D{11};

plot(epoch,ph_accl,epoch,v_accl);
legend('phenotype acclivity', 'virtual acclivity', 'Location', 'nw');
fig = gcf;

%set(fig,'PaperUnits','inches');
%set(fig,'PaperPosition',[0.5 0.5 8 2]);
%print('-dpdf', strcat(filename, '.pdf'));

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print('-dpdf', '-r0', output_filename);
