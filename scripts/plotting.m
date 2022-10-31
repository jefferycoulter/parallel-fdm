data = load('../data/data.csv');
data(end,:) = [];
data = reshape(data, [40, 100, 100, 100]);

m = data(20,:);
m = reshape(m, [100, 100, 100]);
slice(m, [50], [50], [50])