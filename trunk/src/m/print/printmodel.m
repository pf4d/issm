function printmodel(filename,format,varargin)
%PRINTMODEL - save an image of current figure
%
%   filename: output name of image file (no extension)
%   format: image format (ex: 'tiff','jpg','pdf') 
%
%   List of options to printfmodel: 
%
%   figure: number of figure to print (default: current figure)
%   resolution: use higher resolution to anti-alias (default 2)
%   margin: add margin around final image  
%   marginsize: size of margin around final image (default 5)
%   frame: add frame around final image
%   framesize: size of frame around final image (default 5)
%   framecolor: color of frame around final image (default 'black')
%   trim: trim empty space around image (default 'off')
%   hardcopy: 'off' to impose MATLAB to use the same colors (default 'off')
%   
%   Usage:
%      printmodel(filename,format,varargin);
%
%   Examples:
%      printmodel('image','tiff')
%      printmodel('image','eps','margin','on','frame','on','hardcopy','on')

%get options: 
options=pairoptions(varargin{:});

%set defaults
options=addfielddefault(options,'figure','gcf');
options=addfielddefault(options,'format','tiff');
options=addfielddefault(options,'resolution',1);
options=addfielddefault(options,'margin','on');
options=addfielddefault(options,'marginsize',25);
options=addfielddefault(options,'frame','on');
options=addfielddefault(options,'framesize',3);
options=addfielddefault(options,'framecolor','black');
options=addfielddefault(options,'trim','on');
options=addfielddefault(options,'hardcopy','off');

%get fig: 
fig=getfieldvalue(options,'figure');
if ischar(fig),
	fig=gcf;
else
	figure(fig);
	fig=gcf;
end

%In auto mode, MATLAB prints the figure the same size as it appears on the computer screen, centered on the page
set(fig, 'PaperPositionMode', 'auto');

%InvertHardcopy off imposes MATLAB to use the same colors
set(fig, 'InvertHardcopy', getfieldvalue(options,'hardcopy'));

%we could have several formats, as a cell array of strings.
formats=format;
if ~iscell(formats),
	formats={formats};
end

%loop on formats:
for i=1:length(formats),
	format=formats{i};

	%Use higher resolution to anti-alias and use zbuffer to have smooth colors
	print(fig, '-zbuffer','-dtiff',['-r' num2str(get(0,'ScreenPixelsPerInch')*getfieldvalue(options,'resolution'))],filename);

	%some trimming involved? 
	if ~strcmpi(format,'pdf'),
		if strcmpi(getfieldvalue(options,'trim'),'on'),
			system(['convert -trim ' filename '.tif ' filename '.tif']);
		end
	end

	%margin?
	if ~strcmpi(format,'pdf'),
		if strcmpi(getfieldvalue(options,'margin'),'on'),
			marginsize=getfieldvalue(options,'marginsize');
			system(['convert -border ' num2str(marginsize) 'x' num2str(marginsize) ' -bordercolor "white" ' filename '.tif ' filename '.tif']);
		end
	end

	%frame?
	if ~strcmpi(format,'pdf'),
		if strcmpi(getfieldvalue(options,'frame'),'on'),
			framesize=getfieldvalue(options,'framesize');
			framecolor=getfieldvalue(options,'framecolor');
			system(['convert -border ' num2str(framesize) 'x' num2str(framesize) ' -bordercolor "' framecolor '" ' filename '.tif ' filename '.tif']);
		end
	end

	%convert image to correct format
	if ~strcmpi(format,'tiff') & ~strcmpi(format,'tif'),
		system(['convert ' filename '.tif ' filename '.' format]);
		delete([ filename '.tif']);
	end
end
