inputdir = 'K:\Google Drive\DESI\shaft_images\Small Motions\3500\full';
outputdir = 'K:\Google Drive\DESI\shaft_images\Small Motions\3500\crop';
filelist = dir(fullfile(inputdir));
for i=3:length(filelist)
	[~,name] = fileparts(filelist(i).name);
	image = imread(fullfile(inputdir,filelist(i).name));
	image_crop = image(70:3080,730:3730,:);
%  	imshow(image_crop)
	imwrite(image_crop, fullfile(outputdir,[name '_crop.jpg']))
end