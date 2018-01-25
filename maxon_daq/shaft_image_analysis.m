inputdir = 'K:\Google Drive\DESI\shaft_images\20150730 Small Motions\3500\crop';
filelist = dir(inputdir);
angle = NaN(length(filelist)-4, 1);
for i=3 : length(filelist)-2
	
	before = imread(fullfile(inputdir,filelist(i).name));
	after = imread(fullfile(inputdir,filelist(i+1).name));

	original = rgb2gray(before);
	distorted = rgb2gray(after);
	%% Find Matching Features Between Images
	% Detect features in both images.
	
	ptsOriginal  = detectSURFFeatures(original);
	ptsDistorted = detectSURFFeatures(distorted);

	%%
	% Extract feature descriptors.
	[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
	[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);

	%%
	% Match features by using their descriptors.
	indexPairs = matchFeatures(featuresOriginal, featuresDistorted);

	%%
	% Retrieve locations of corresponding points for each image.
	matchedOriginal  = validPtsOriginal(indexPairs(:,1));
	matchedDistorted = validPtsDistorted(indexPairs(:,2));

	%%
% 	Show putative point matches.
% 	figure;
% 	showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
% 	title('Putatively matched points (including outliers)');

	%% Estimate Transformation
	% Find a transformation corresponding to the matching point pairs using the
	% statistically robust M-estimator SAmple Consensus (MSAC) algorithm, which
	% is a variant of the RANSAC algorithm. It removes outliers while computing
	% the transformation matrix. You may see varying results of the
	% transformation computation because of the random sampling employed by the
	% MSAC algorithm.
	[tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(...
		matchedDistorted, matchedOriginal, 'similarity');

	%%
	% Display matching point pairs used in the computation of the
	% transformation.
% 	figure; 
% 	showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
% 	title('Matching points (inliers only)');
% 	legend('ptsOriginal','ptsDistorted');

	%% Solve for Scale and Angle
	% Use the geometric transform, tform, to recover the scale and angle.
	% Since we computed the transformation from the distorted to the original
	% image, we need to compute its inverse to recover the distortion.
	%
	%  Let sc = s*cos(theta)
	%  Let ss = s*sin(theta)
	%
	%  Then, Tinv = [sc -ss  0;
	%                ss  sc  0;
	%                tx  ty  1]
	%
	%  where tx and ty are x and y translations, respectively.
	%

	%%
	% Compute the inverse transformation matrix.
	Tinv  = tform.invert.T;

	ss = Tinv(2,1);
	sc = Tinv(1,1);
	scaleRecovered = sqrt(ss*ss + sc*sc)
	thetaRecovered = atan2(ss,sc)*180/pi

	%%
	% The recovered values should match your scale and angle values selected in
	% *Step 2: Resize and Rotate the Image*.

	%% Recover the Original Image
	% Recover the original image by transforming the distorted image.
	%%
	% Compare |recovered| to |original| by looking at them side-by-side in a
	% montage.
% 	outputView = imref2d(size(original));
% 	recovered  = imwarp(distorted,tform,'OutputView',outputView);
% 	figure, imshowpair(original,recovered,'montage')

	%%
	% The |recovered| (right) image quality does not match the |original|
	% (left) image because of the distortion and recovery process. In
	% particular, the image shrinking causes loss of information. The artifacts
	% around the edges are due to the limited accuracy of the transformation.
	% If you were to detect more points in *Step 4: Find Matching Features
	% Between Images*, the transformation would be more accurate. For example,
	% we could have used a corner detector, detectFASTFeatures, to complement
	% the SURF feature detector which finds blobs. Image content and image size
	% also impact the number of detected features.

% 	displayEndOfDemoMessage(mfilename)
	angle(i-2) = thetaRecovered;
end