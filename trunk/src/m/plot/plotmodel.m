function plotmodel(md,varargin)
%At command prompt, type plotdoc for help on documentation.

%First process options
options=plotoptions(varargin{:});

%get number of subplots
subplotwidth=ceil(sqrt(options.numberofplots));

%Get figure number and number of plots
figurenumber=options.figurenumber;
numberofplots=options.numberofplots;

%if nlines and ncols specified, then bypass.
if exist(options.list{1},'nlines'),
	nlines=getfieldvalue(options.list{1},'nlines');
else 
	nlines=ceil(numberofplots/subplotwidth);
end

if exist(options.list{1},'ncols'),
	ncols=getfieldvalue(options.list{1},'ncols');
else 
	ncols=subplotwidth;
end

%check that nlines and ncols were given at the same time!
if ((exist(options.list{1},'ncols') & ~exist(options.list{1},'ncols')) | (~exist(options.list{1},'ncols') & exist(options.list{1},'ncols')))
	error('plotmodel error message: nlines and ncols  need to be specified together, or not at all');
end

%go through subplots
if numberofplots,

	%Create figure 
	f=figure(figurenumber);clf;
	if strcmpi(getfieldvalue(options.list{1},'visible','on'),'off'),
		set(f,'Visible','Off');
	end

	if exist(options.list{1},'figposition'), % {{{
		figposition=getfieldvalue(options.list{1},'figposition');
		if ischar(figposition),
			if strcmpi(figposition,'larour'),
				set(gcf,'Position',[1604 4 1594 1177]);
			elseif strcmpi(figposition,'larour2'),
				set(gcf,'Position',[756    62   827   504]);
			elseif strcmpi(figposition,'mathieu'),
				set(gcf,'Position',[300 1 1580 1150]);
			elseif strcmpi(figposition,'fullscreen'),
				set(gcf,'Position',get(0,'ScreenSize'));
			elseif strcmpi(figposition,'halfright'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				set(gcf,'Position',fix([left+widt/2 bott widt/2 heig]));
			elseif strcmpi(figposition,'halfleft'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				set(gcf,'Position',fix([left bott widt/2 heig]));
			elseif strcmpi(figposition,'square'),
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=min(screen(3)-25,screen(4)-25);
				set(gcf,'Position',fix([left+(screen(3)-widt) bott widt widt]));
			elseif strcmpi(figposition,'portrait'),
				%reformat with letter paper size (8.5" x 11")
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				portrait=fix([left+widt-(heig*8.5/11) bott heig*8.5/11 heig]);
				set(gcf,'Position',portrait)
			elseif strcmpi(figposition,'landscape'),
				%reformat with letter paper size (8.5" x 11")
				screen=get(0,'ScreenSize');
				left=screen(1); bott=screen(2); widt=screen(3); heig=screen(4)-25;
				landscape=fix([left+widt-(heig*11/8.5) bott heig*11/8.5 heig]);
				set(gcf,'Position',landscape)
			else
				disp('''figposition'' string not supported yet');
			end
		else
			set(gcf,'Position',figposition);
		end
	end % }}}

	%Use zbuffer renderer (snoother colors) and white background
	set(f,'Renderer','zbuffer','color',getfieldvalue(options.list{1},'figurebackgroundcolor','w'));

	%Go through all data plottable and close window if an error occurs
	try,
		for i=1:numberofplots,
			plot_manager(getfieldvalue(options.list{i},'model',md),options.list{i},subplotwidth,nlines,ncols,i);
			%List all unused options
			displayunused(options.list{i})
		end
	catch me,
		%figure(figurenumber),close;
		rethrow(me);
	end
else
	error('plotmodel error message: no output data found. ');
end
