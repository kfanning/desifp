%% read file

filename = 'K:\Google Drive\DESI\tools_manuals\Andor iKon\magic_n_0000.fits';
filename = 'K:\Google Drive\DESI\tools_manuals\Andor iKon\filter2.fits';
info = fitsinfo(filename);
fitsdisp(filename);
data = fitsread(filename,'primary');
dataRenorm = data/max(data(:));
displayRange = [min(data(:)), max(data(:))];

% % display
% colormap gray;
% imagesc(data);
% % axis(data);
% iptsetpref('ImshowAxesVisible','on')
% imshow(data, displayRange, 'InitialMagnification', 30)

%% ID

% [optimizer, metric]  = imregconfig('monomodal');
% 0.1 for single source, 0.25? for mesh
dataBW = im2bw(dataRenorm,0.30);
ellipseStats = regionprops(dataBW, 'Centroid', 'MajorAxisLength', ...
	'MinorAxisLength', 'Eccentricity', 'Area','PixelIdxList');

%% processing

% extract centroid for plotting
centroid = cat(1, ellipseStats.Centroid); 
% separate x and y positions
% ellipseStats.CentroidX = ellipseStats.Centroid(1);
% ellipseStats.CentroidY = ellipseStats.Centroid(2);


% mean of background intensity in mask 0-valued region including
% everywhere outside the elliptic ROIs
IntensityBgd = mean2(data(dataBW<1));
% clear background from data
dataCleaned = data - IntensityBgd;
dataCleanedMasked = dataCleaned.*dataBW;

% statistics of elliptic ROIs
% fill ellipseStats "array of struct" - pita
for i=1:size(centroid,1)
	valuesROI = dataCleanedMasked(ellipseStats(i).PixelIdxList);
	% max intensity in ROI;
	ellipseStats(i).IntensityMax = max(valuesROI);
	% mean intensity in ROI
	ellipseStats(i).IntensityMean = mean2(valuesROI);
	% total flux through ROI
	ellipseStats(i).IntensityTotal = sum(valuesROI);
end

%% show images for comparision

% cleaned data above background
imagesc(dataCleanedMasked); colormap(gray)
% BW mask
imagesc(dataBW); colormap(gray)

% data original
iptsetpref('ImshowAxesVisible','on')
imshow(data, displayRange, 'InitialMagnification', 30)
hold on
plot(centroid(:,1), centroid(:,2), 'b*')
hold off

% data BW
figure
iptsetpref('ImshowAxesVisible','on')
imshow(dataBW, 'InitialMagnification', 30)
hold on
plot(centroid(:,1), centroid(:,2), 'b*')
hold off

%% plots

% surface plots of data cleaned
[x, y] = meshgrid(1:1:2048);
plotSurf = surf(x,y,dataCleaned);
set(plotSurf,'LineStyle','none')
colormap jet

% contour plot
figure;
% meanThreshold = max(cat(1, ellipseStats.IntensityMean));
% plotContour = contour(dataCleaned,[meanThreshold,meanThreshold]);
plotContour = contour(dataCleaned,3);

% pseudocolor (checkerboard) plot
figure;
plotColor = pcolor(x,y,dataCleanedMasked);
set(plotColor,'LineStyle','none')
colorbar

%% print results
disp(sprintf('\n'));
disp(ellipseStats)
disp(sprintf('\n'));
