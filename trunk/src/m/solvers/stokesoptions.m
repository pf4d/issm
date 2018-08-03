function stokes=stokesoptions(varargin)
%STOKESOPTIONS - return STOKES multi-physics solver petsc options
%
%   Usage:
%      options=stokesoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
stokes=struct();

stokes.toolkit='petsc';
stokes.mat_type=getfieldvalue(options,'mat_type','mpiaij');
stokes.ksp_max_it=getfieldvalue(options,'ksp_max_it',1000);
stokes.ksp_type=getfieldvalue(options,'ksp_type','gmres');
stokes.pc_type=getfieldvalue(options,'pc_type','fieldsplit');
stokes.pc_field_split_type=getfieldvalue(options,'pc_field_split_type','schur');
stokes.fieldsplit_0_pc_type=getfieldvalue(options,'fieldsplit_0_pc_type','hypre');
stokes.fieldsplit_0_ksp_type=getfieldvalue(options,'fieldsplit_0_ksp_type','gmres');
stokes.fieldsplit_0_pc_hypre_type=getfieldvalue(options,'fieldsplit_0_pc_hypre_type','boomerang');
stokes.fieldsplit_1_pc_type=getfieldvalue(options,'fieldsplit_1_pc_type','jacobi');
stokes.fieldsplit_1_ksp_type=getfieldvalue(options,'fieldsplit_1_ksp_type','preonly');
stokes.issm_option_solver=getfieldvalue(options,'issm_option_solver','stokes');
