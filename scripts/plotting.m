data = load('../data/data_uv.csv');
data(end,:) = [];
data = reshape(data, [40, 100, 100, 100]);

video = VideoWriter("../graphics/animations/uv3d", 'MPEG-4');
video.FrameRate = 10;
open(video);

for t = 1:size(data,1)
    m = data(t,:);
    m = reshape(m, [100, 100, 100]);
    colormap("default")
    s = slice(m, 50, 50, 50);
    %s.FaceColor = "interp";
    colorbar
    %clim([0 10^-6])
    pause(0.005)
    frame = getframe(gcf);

    writeVideo(video, frame);
    
    drawnow
end

close(video)