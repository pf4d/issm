function [element_partitioning, node_partitioning] = MeshPartition(md.mesh,numpartitions);
%MESHPARTITION - Partition mesh according to the number of areas, using Metis library.
%
%	   Usage:
%			[element_partitioning,node_partitioning]=MeshPartition(md.mesh,numpartitions)");
%
%	   element_partitioning: Vector of partitioning area numbers, for every element.
%	   node_partitioning: Vector of partitioning area numbers, for every node.

% Check usage
if nargin~=2
	help MeshPartition
	error('Wrong usage (see above)');
end

% Call mex module
[element_partitioning, node_partitioning] = MeshPartition_matlab(md.mesh,numpartitions);
