%
%  definition for the uniform_uncertain class.
%
%  [uuv]=uniform_uncertain(varargin)
%
%  where the required varargin are:
%    descriptor    (char, description, '')
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
classdef uniform_uncertain
    properties
        descriptor='';
        lower     =-Inf;
        upper     = Inf;
    end

    methods
        function [uuv]=uniform_uncertain(varargin)

            switch nargin

%  create a default object

                case 0

%  copy the object

                case 1
                    if isa(varargin{1},'uniform_uncertain')
                        uuv=varargin{1};
                    else
                        error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
                            inputname(1),class(varargin{1}),'uniform_uncertain');
                    end

%  not enough arguments

                case 2
                    error('Construction of ''%s'' class object requires at least %d inputs.',...
                        'uniform_uncertain',3)

%  create the object from the input

                otherwise
                    asizec=num2cell(array_size(varargin{1:min(nargin,3)}));
                    uuv(asizec{:})=uniform_uncertain;
                    clear asizec

                    if ischar(varargin{1})
                        varargin{1}=cellstr(varargin{1});
                    end
                    for i=1:numel(uuv)
                        if (numel(varargin{1}) > 1)
                            uuv(i).descriptor=varargin{1}{i};
                        else
                            uuv(i).descriptor=[char(varargin{1}) string_dim(uuv,i,'vector')];
                        end
                        if (numel(varargin{2}) > 1)
                            uuv(i).lower     =varargin{2}(i);
                        else
                            uuv(i).lower     =varargin{2};
                        end
                        if (numel(varargin{3}) > 1)
                            uuv(i).upper     =varargin{3}(i);
                        else
                            uuv(i).upper     =varargin{3};
                        end
                    end
            end

        end

        function []=disp(uuv)

%  display the object

            disp(sprintf('\n'));
            for i=1:numel(uuv)
                disp(sprintf('class ''%s'' object ''%s%s'' = \n',...
                    class(uuv),inputname(1),string_dim(uuv,i)));
                disp(sprintf('    descriptor: ''%s'''  ,uuv(i).descriptor));
                disp(sprintf('         lower: %g'      ,uuv(i).lower));
                disp(sprintf('         upper: %g\n'    ,uuv(i).upper));
            end

        end

        function [desc]  =prop_desc(uuv,dstr)
            desc=cell(1,numel(uuv));
            for i=1:numel(uuv)
                if ~isempty(uuv(i).descriptor)
                    desc(i)=cellstr(uuv(i).descriptor);
                elseif ~isempty(inputname(1))
                    desc(i)=cellstr([inputname(1) string_dim(uuv,i,'vector')]);
                elseif exist('dstr','var')
                    desc(i)=cellstr([dstr         string_dim(uuv,i,'vector')]);
                else
                    desc(i)=cellstr(['uuv'        string_dim(uuv,i,'vector')]);
                end
            end
            desc=allempty(desc);
        end
        function [initpt]=prop_initpt(uuv)
            initpt=[];
        end
        function [lower] =prop_lower(uuv)
            lower=zeros(1,numel(uuv));
            for i=1:numel(uuv)
                lower(i)=uuv(i).lower;
            end
            lower=allequal(lower,-Inf);
        end
        function [upper] =prop_upper(uuv)
            upper=zeros(1,numel(uuv));
            for i=1:numel(uuv)
                upper(i)=uuv(i).upper;
            end
            upper=allequal(upper, Inf);
        end
        function [mean]  =prop_mean(uuv)
            mean=[];
        end
        function [stddev]=prop_stddev(uuv)
            stddev=[];
        end
        function [initst]=prop_initst(uuv)
            initst=[];
        end
        function [stype] =prop_stype(uuv)
            stype={};
        end
        function [scale] =prop_scale(uuv)
            scale=[];
        end
    end

    methods (Static)
        function []=dakota_write(fidi,dvar)

%  collect only the variables of the appropriate class

            uuv=struc_class(dvar,'uniform_uncertain');

%  write variables

            vlist_write(fidi,'uniform_uncertain','uuv',uuv);
        end
    end
end
