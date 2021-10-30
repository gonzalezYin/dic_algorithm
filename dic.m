% Read the images
% image1 = imread('testdata2/0019.tif');
% image2 = imread('testdata2/0020.tif');
% Get image 
image_folder = "/home/dulab/Desktop/noisy_simulation_data/hexagon01/noisy01/";
imds = imageDatastore(image_folder, "FileExtensions", ".tif");
image1 = readimage(imds, 1);
image2 = readimage(imds, 2);

% Find the SURF features
points1 = detectSURFFeatures(image1);
points2 = detectSURFFeatures(image2);

% Extract the features: Returns extracted feature vectors/descriptors, 
% and their corresponding locations. 
[feature1, valid_points1] = extractFeatures(image1, points1);
[feature2, valid_points2] = extractFeatures(image2, points2);

% Retrieve the locations of matched points
indexPairs = matchFeatures(feature1, feature2);
matchedPoints1 = valid_points1(indexPairs(:, 1));
matchedPoints2 = valid_points2(indexPairs(:, 2));
featurePoints1 = feature1(indexPairs(:, 1),:);
featurePoints2 = feature2(indexPairs(:, 2),:);

% Display the matching points
figure;
showMatchedFeatures(image1, image2, matchedPoints1, matchedPoints2);
legend('matched points 1, matched points 2');

% Get the Eculidean distance 
numOfPairs = matchedPoints1.length();
dist = zeros(numOfPairs, 2);
for index = 1: numOfPairs
    featureVector = (featurePoints1(index,:) - featurePoints2(index,:)).^2;
    sumOfFeatureVector = sum(featureVector);
    dist(index, :) = [sqrt(sumOfFeatureVector), index];
    
end
sortedFeaturePointDist = sortrows(dist);
indexOfPoints = sortedFeaturePointDist(:, 2);

% Get best-matched feature vector queue: (di, xi, yi, ui, vi, өi)
featurePointQueue = zeros(numOfPairs, 6);
for index = 1: numOfPairs
    tempPoint = matchedPoints1(indexOfPoints(index));
    featurePointQueue(index,:) = [sortedFeaturePointDist(index, 1), tempPoint.Location(1), tempPoint.Location(2), 0, 190, 0];
end

% Get the best-matched feature point (the first loop begins) 
best_matched_point = featurePointQueue(1, :);
d1 = best_matched_point(1);
x1 = best_matched_point(2);
y1 = best_matched_point(3);
u1 = best_matched_point(4);
v1 = best_matched_point(5);
rotation1 = best_matched_point(6);


% Get the calculation point distribution: per 5*5=25 grids has a calculation point
numOfRow = floor(size(image1, 1) / 5);
numOfCol = floor(size(image1, 2) / 5);
n = numOfRow * numOfCol;
deltaRow = -2;
deltaCol = -2;
calPoint = zeros(n, 2);
indexOfCalPoint = 1; 
for row = 1: numOfRow
    deltaRow = deltaRow + 5;
    deltaCol = -2;
    for col = 1: numOfCol
        deltaCol = deltaCol + 5;
        calPoint(indexOfCalPoint, :) = [deltaRow, deltaCol];
        indexOfCalPoint = indexOfCalPoint + 1;
    end
end

% Get the seed point
min_dist = 1000000;
seed_point = [1, 1];
for index = 1: n
    dist_vector = (calPoint(index, :) - [x1, y1]).^2;
    temp_dist = sqrt(sum(dist_vector));
    if temp_dist < min_dist
        min_dist = temp_dist;
        seed_point = calPoint(index, :);
    end
end
xseed = seed_point(1);
yseed = seed_point(2);

% Initial guess
u0seed = u1 + (cos(rotation1)-1) * (xseed-x1) + sin(rotation1) * (yseed-y1);
u0xseed = cos(rotation1)-1;
u0yseed = sin(rotation1);
v0seed = v1 - sin(rotation1) * (xseed-x1) + (cos(rotation1)-1) * (yseed-y1);
v0xseed = -1 * sin(rotation1);
v0yseed = cos(rotation1)-1;
p0seed = [u0seed, u0xseed, u0yseed, v0seed, v0xseed, v0yseed];

% p0seed -> pseed (optimized by IC-GN: advanced sub_pixel registration algorithm)
pseed = p0seed;

c1 = c1';   c2 = c2';
% Concatenate images
[h1,w1,~] = size(image1);
[h2,~,~] = size(image2);
if h1 > h2
    image2 = padarray(image2, [h1-h2 0], 'post');
else
    image1 = padarray(image1, [h2-h1 0], 'post');
end
imc = cat(2, image1, image2);

% figure;
imshow(imc);
hold on; axis image off;
c2(1,:) = c2(1,:) + w1; % adjust the x coordinate by width
plot(c1(1,:), c1(2,:), 'r.');
plot(c2(1,:), c2(2,:), 'r.');
for i = 1:size(c1,2)    
    line([c1(1,i) c2(1,i)], ([c1(2,i) c2(2,i)]),'color','g','LineWidth',1);
end




%imshow(image1);
%hold on;
%plot(points1);
%figure;
%imshow(image1);
%hold on;
%plot(valid_points2);

