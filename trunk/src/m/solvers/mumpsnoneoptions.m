function mumps=mumpsoptions(varargin)
%MUMPSOPTIONS - return MUMPS direct solver  petsc options
%
%   Usage:
%      options=mumpsoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
mumps=struct();

%default mumps options
PETSC_VERSION=IssmConfig('_PETSC_MAJOR_');
if PETSC_VERSION==2.,
	mumps.toolkit='petsc';
	mumps.mat_type=getfieldvalue(options,'mat_type','aijmumps');
	mumps.ksp_type=getfieldvalue(options,'ksp_type','preonly');
	mumps.pc_type=getfieldvalue(options,'pc_type','lu');
	mumps.mat_mumps_icntl_14=getfieldvalue(options,'mat_mumps_icntl_14',120);
	mumps.pc_factor_shift_positive_definite=getfieldvalue(options,'pc_factor_shift_positive_definite','true');
end

if PETSC_VERSION==3.,
	mumps.toolkit='petsc';
	mumps.mat_type=getfieldvalue(options,'mat_type','mpiaij');
	mumps.ksp_type=getfieldvalue(options,'ksp_type','preonly');
	mumps.pc_type=getfieldvalue(options,'pc_type','lu');
	mumps.pc_factor_mat_solver_package=getfieldvalue(options,'pc_factor_mat_solver_package','mumps');
	mumps.mat_mumps_icntl_14=getfieldvalue(options,'mat_mumps_icntl_14',120);
	mumps.pc_factor_shift_positive_definite=getfieldvalue(options,'pc_factor_shift_positive_definite','true');
	mumps.mat_mumps_icntl_28=2; %1=serial, 2=parallel
	mumps.mat_mumps_icntl_29=2; %parallel ordering 1 = ptscotch, 2 = parmetis
end
