data = load('../data/data.csv');
data(end,:) = [];
data = reshape(data, [40, 100, 100, 100]);

for t = 1:size(data,1)
    m = data(t,:);
    m = reshape(m, [100, 100, 100]);
    colormap("default")
    s = slice(m, 50, 50, 50);
    %s.FaceColor = "interp";
    colorbar
    clim([-10 10])
    drawnow
end