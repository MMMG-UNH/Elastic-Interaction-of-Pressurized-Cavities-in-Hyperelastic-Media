clc;clear;close all;

base_dir = 'D:\ElasticInteraction\files';
P = 5;
d_min = 0.22;
d_max = 0.22;
d_increment = 0.02;
R = 0.1;
d = d_min:d_increment:d_max;
dR = d/R;

% Define the number of files to read
numFiles = length(dR);
% Define size of the file to store area values within the loop
DeformedAreas = zeros(1,numFiles);
W = zeros(1,numFiles);
% Preallocate a cell array to store the data
data = cell(1, numFiles);
% Iterate over the files
for fileIndex = 1:numFiles
    % Construct the file name
    fileName = sprintf('deformed_coordinates-%d.csv', fileIndex);
    % Read the data from the file
    filePath = fullfile(base_dir, fileName);
    data{fileIndex} = readtable(filePath);
    % Coordinates of the discrete points (x, y)
    raw = data{fileIndex};
    Xcoord = raw(:, 2);
    points(:,1) = Xcoord.XCoordinate;
    Ycoord = raw(:, 3);
    points(:,2) = Ycoord.YCoordinate;
    % Triangulate the points
    Tri = delaunayTriangulation(points(:,1), points(:,2));
    tri = Tri.ConnectivityList;
    % Plot the shape and its comprising triangles
    figure;
    axis equal;
    hold on
    triplot(Tri);
    hold off
    % Calculate triangle areas
    v1 = [points(tri(:, 1), :) ones(size(tri, 1), 1)];
    v2 = [points(tri(:, 2), :) ones(size(tri, 1), 1)];
    v3 = [points(tri(:, 3), :) ones(size(tri, 1), 1)];
    triangle_areas = 0.5 * abs(cross(v2 - v1, v3 - v1));
    % Calculate the total surface area
    surface_area = sum(triangle_areas);
    disp(['Surface Area = ' num2str(surface_area(3))]);
    DeformedAreas(1,fileIndex) = surface_area(3);
end