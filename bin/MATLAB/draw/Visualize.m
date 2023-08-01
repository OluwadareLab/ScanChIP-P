function [ output_args ] = Visualize( Data )
% This function allows for the visualization of data in Heat Map
colormap('sky');
imagesc(Data);
colorbar;
xlabel('Genomic bin (resolution: 40kb)');
ylabel('Genomic bin (resolution: 40kb)');
title('HDBSCAN - Cross Features (len(matrix)-4)');
end

