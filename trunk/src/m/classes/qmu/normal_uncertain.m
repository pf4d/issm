%
%  definition for the normal_uncertain class.
%
%  [nuv]=normal_uncertain(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
%    mean          (double, mean, NaN)
%    stddev        (double, standard deviation, NaN)
%  and the optional varargin and defaults are:
%    lower         (double, lower bound, -Inf)
%    upper         (double, upper bound,  Inf)
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and three or more
%  arguments constructs a new instance from the arguments.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
classdef normal_uncertain
    properties
        descriptor='';
        mean      = NaN;
        stddev    = NaN;
        lower     =-Inf;
        upper     = Inf;
    end

    methods
        function [nuv]=normal_uncertain(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'normal_uncertain')
                        nuv=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'normal_uncertain');
                    end

%  not enough arguments

                case 2
                    error('Construction of ''%s'' class object requires at least %d inputs.',...
                        'normal_uncertain',3)

%  create the object from the input

                otherwise
                    asizec=num2cell(array_size(varargin{1:min(nargin,5)}));
                    nuv(asizec{:})=normal_uncertain;
                    clear asizec

                    if ischar(varargin{1})
                        varargin{1}=cellstr(varargin{1});
                    end
                    for i=1:numel(nuv)
                        if (numel(varargin{1}) > 1)
                            nuv(i).descriptor=varargin{1}{i};
                        else
                            if numel(nuv)==1,
								nuv(i).descriptor=char(varargin{1});
							else
								nuv(i).descriptor=[char(varargin{1}) num2str(i)];
							end
                        end
                        if (numel(varargin{2}) > 1)
                            nuv(i).mean      =varargin{2}(i);
                        else
                            nuv(i).mean      =varargin{2};
                        end
                        if (numel(varargin{3}) > 1)
                            nuv(i).stddev    =varargin{3}(i);
                        else
                            nuv(i).stddev    =varargin{3};
                        end
                    end

                    if (nargin >= 4)
                        for i=1:numel(nuv)
                            if (numel(varargin{4}) > 1)
                                nuv(i).lower     =varargin{4}(i);
                            else
                                nuv(i).lower     =varargin{4};
                            end
                        end
                        if (nargin >= 5)
                            for i=1:numel(nuv)
                                if (numel(varargin{5}) > 1)
                                    nuv(i).upper     =varargin{5}(i);
                                else
                                    nuv(i).upper     =varargin{5};
                                end
                            end
                            if (nargin > 5)
                                warning('normal_uncertain:extra_arg',...
                                    'Extra arguments for object of class ''%s''.',...
                                    class(nuv));
                            end
                        end
                    end
            end

        end

        function []=disp(nuv)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(nuv)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(nuv),inputname(1),string_dim(nuv,i)));
                disp(sprintf('    descriptor: ''%s'''  ,nuv(i).descriptor));
                disp(sprintf('          mean: %g'      ,nuv(i).mean));
                disp(sprintf('        stddev: %g'      ,nuv(i).stddev));
                disp(sprintf('         lower: %g'      ,nuv(i).lower));
                disp(sprintf('         upper: %g\n'    ,nuv(i).upper));
            end

        end

        function [desc]  =prop_desc(nuv,dstr)
            desc=cell(1,numel(nuv));
            for i=1:numel(nuv)
                if ~isempty(nuv(i).descriptor)
                    desc(i)=cellstr(nuv(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(nuv,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(nuv,i,'vector')]);
                else
                    desc(i)=cellstr(['nuv'        string_dim(nuv,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [initpt]=prop_initpt(nuv)
            initpt=[];
        end
        function [lower] =prop_lower(nuv)
            lower=zeros(1,numel(nuv));
            for i=1:numel(nuv)
                lower(i)=nuv(i).lower;
            end
            lower=allequal(lower,-Inf);
        end
        function [upper] =prop_upper(nuv)
            upper=zeros(1,numel(nuv));
            for i=1:numel(nuv)
                upper(i)=nuv(i).upper;
            end
            upper=allequal(upper, Inf);
        end
        function [mean]  =prop_mean(nuv)
            mean=zeros(1,numel(nuv));
            for i=1:numel(nuv)
                mean(i)=nuv(i).mean;
            end
        end
        function [stddev]=prop_stddev(nuv)
            stddev=zeros(1,numel(nuv));
            for i=1:numel(nuv)
                stddev(i)=nuv(i).stddev;
            end
        end
        function [initst]=prop_initst(nuv)
            initst=[];
        end
        function [stype] =prop_stype(nuv)
            stype={};
        end
        function [scale] =prop_scale(nuv)
            scale=[];
        end
    end

    methods (Static)
        function []=dakota_write(fidi,dvar)

%  collect only the variables of the appropriate class

            nuv=struc_class(dvar,'normal_uncertain');

%  write variables

            vlist_write(fidi,'normal_uncertain','nuv',nuv);
        end
    end
end
