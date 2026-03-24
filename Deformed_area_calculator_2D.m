% =========================================================================
% Deformed_area_calculator_2D.m
%
% Computes the deformed cross-sectional area of the left pressurized cavity
% from finite element node coordinates extracted via ABAQUS post-processing.
%
% Purpose:
%   The total potential energy of the system is (Eq. 9 of paper):
%       PE = SE - P * ΔA
%   where SE is the total strain energy and ΔA is the change in cavity area
%   due to pressurization. This script computes the deformed area A_deformed
%   for each separation distance, from which ΔA = A_deformed - pi*R^2.
%
% Method:
%   For each separation distance, the deformed cavity boundary is defined by
%   a set of nodes whose coordinates were extracted from ABAQUS. These nodes
%   are connected in no particular order, so Delaunay triangulation is used
%   to triangulate the enclosed region, and the total area is computed by
%   summing triangle areas using the cross-product formula.
%
% Inputs (from ABAQUS post-processing):
%   deformed_coordinates-{index}.csv — CSV files containing node IDs and
%   deformed (X, Y) coordinates of the left cavity boundary nodes, produced
%   by Section 3 of Two_pressurized_holes_ABAQUS_script_2D.py.
%
%   File naming: deformed_coordinates-1.csv corresponds to the smallest
%   separation (d = d_min), and the index increments with separation distance.
%
% Parameters (must match those used in the ABAQUS simulation):
%   base_dir    : directory containing the deformed_coordinates CSV files
%   P           : applied internal pressure
%   d_min       : minimum center-to-center distance
%   d_max       : maximum center-to-center distance
%   d_increment : step in center-to-center distance
%   R           : undeformed hole radius
%
% Output:
%   DeformedAreas (1 x numFiles): deformed cavity area at each separation
%
% Note: The Delaunay triangulation assumes the node coordinates form a
% simply connected region (a disk). If the cavity deforms significantly,
% check that the triangulation correctly fills the cavity interior.
%
% Reference:
%   Saeedi & Kothari (2025). J. Appl. Mech., 92(5), 051008.
%   DOI: 10.1115/1.4067855
% =========================================================================

clc; clear; close all;

% -------------------------------------------------------------------------
% Parameters (must match the ABAQUS simulation settings)
% -------------------------------------------------------------------------
base_dir    = 'YOUR PATH\files';
P           = 5;
d_min       = 0.22;
d_max       = 2.0;
d_increment = 0.02;
R           = 0.1;

% Derive the array of center-to-center distances and their normalized values
d  = d_min:d_increment:d_max;
dR = d / R;  % non-dimensional center-to-center distance (eta/R)

% -------------------------------------------------------------------------
% Pre-allocate output arrays
% -------------------------------------------------------------------------
numFiles     = length(dR);
DeformedAreas = zeros(1, numFiles);  % deformed area of the left cavity
W             = zeros(1, numFiles);  % (placeholder; not currently used)
data          = cell(1, numFiles);   % raw CSV data for each file

% =========================================================================
% Main loop: read each CSV file and compute the deformed cavity area
% =========================================================================
for fileIndex = 1:numFiles

    % Construct the filename for this separation distance.
    % Files are indexed sequentially starting from 1 (d = d_min).
    fileName = sprintf('deformed_coordinates-%d.csv', fileIndex);
    filePath = fullfile(base_dir, fileName);

    % Read the deformed node coordinates from CSV
    data{fileIndex} = readtable(filePath);
    raw = data{fileIndex};

    % Extract X and Y coordinates into a 2-column points array
    points(:, 1) = raw.XCoordinate;
    points(:, 2) = raw.YCoordinate;

    % --- Delaunay triangulation of the cavity interior ---
    % Since the node coordinates are unordered, direct polygon area formulas
    % are not applicable. Instead, we triangulate the point cloud and sum
    % triangle areas.
    %
    % Note: Delaunay triangulation fills the convex hull of the points. For a
    % nearly circular cavity, this is identical to the cavity interior. For
    % highly deformed cavities, verify visually that no spurious triangles
    % are generated outside the true cavity boundary.
    Tri = delaunayTriangulation(points(:,1), points(:,2));
    tri = Tri.ConnectivityList;  % N_tri x 3 array of vertex indices

    % Visualize the triangulation for this separation distance
    % (useful for verifying mesh quality and cavity shape)
    figure;
    axis equal;
    hold on
    triplot(Tri);
    title(sprintf('Triangulation for dR = %.2f', dR(fileIndex)));
    hold off

    % --- Compute area of each triangle using the cross-product formula ---
    % For triangle with vertices v1, v2, v3:
    %   Area = 0.5 * ||(v2-v1) × (v3-v1)||
    % In 2D, we augment with a z=1 column to use the 3D cross product;
    % the magnitude of the z-component of the cross product gives the area.
    v1 = [points(tri(:,1), :) ones(size(tri,1), 1)];
    v2 = [points(tri(:,2), :) ones(size(tri,1), 1)];
    v3 = [points(tri(:,3), :) ones(size(tri,1), 1)];

    % cross() returns a 3-component vector per row; index 3 is the z-component
    triangle_areas = 0.5 * abs(cross(v2 - v1, v3 - v1));

    % Sum all triangle areas to get the total enclosed area
    % triangle_areas(:,3) is the z-component (the 2D area contribution)
    surface_area = sum(triangle_areas);
    disp(['Surface Area = ' num2str(surface_area(3))]);

    % Store the computed area for this separation distance
    DeformedAreas(1, fileIndex) = surface_area(3);

end  % end loop over separation distances
