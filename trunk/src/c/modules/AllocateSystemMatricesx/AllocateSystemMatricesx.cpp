/*!\file AllocateSystemMatricesx
 * \brief retrieve vector from inputs in elements
 */

#include "./AllocateSystemMatricesx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void AllocateSystemMatricesx(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,Vector<IssmDouble>** pdf,Vector<IssmDouble>** ppf,FemModel* femmodel){

	/*Intermediary*/
	int  fsize,ssize,flocalsize,slocalsize;
	int  connectivity, numberofdofspernode;
	int  configuration_type;
	int  m,n,M,N;
	int *d_nnz = NULL;
	int *o_nnz = NULL;

	/*output*/
	Matrix<IssmDouble> *Kff  = NULL;
	Matrix<IssmDouble> *Kfs  = NULL;
	Vector<IssmDouble> *pf   = NULL;
	Vector<IssmDouble> *df   = NULL;

	bool oldalloc=false;
	char* toolkittype=NULL;

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&connectivity,MeshAverageVertexConnectivityEnum);

	/*retrieve node info*/
	fsize      = femmodel->nodes->NumberOfDofs(configuration_type,FsetEnum);
	ssize      = femmodel->nodes->NumberOfDofs(configuration_type,SsetEnum);
	flocalsize = femmodel->nodes->NumberOfDofsLocal(configuration_type,FsetEnum);
	slocalsize = femmodel->nodes->NumberOfDofsLocal(configuration_type,SsetEnum);

	numberofdofspernode=femmodel->nodes->MaxNumDofs(configuration_type,GsetEnum);

	/*if our matrices are coming from issm, we don't do dynamic allocation like Petsc 
	 * does, and this routine is essentially useless. Force standard alloc in this case: */
	toolkittype=ToolkitOptions::GetToolkitType();

	if(oldalloc){
		if(pKff) Kff=new Matrix<IssmDouble>(fsize,fsize,connectivity,numberofdofspernode);
		if(pKfs) Kfs=new Matrix<IssmDouble>(fsize,ssize,connectivity,numberofdofspernode);
		if(pdf)  df =new Vector<IssmDouble>(fsize);
		if(ppf)  pf =new Vector<IssmDouble>(fsize);
	}
	else{
		if(pKff){
			m=flocalsize; n=flocalsize; /*local  sizes*/
			M=fsize;      N=fsize;      /*global sizes*/
			if(strcmp(toolkittype,"issm")==0){
				Kff=new Matrix<IssmDouble>(m,n,M,N,NULL,NULL);
			}
			else{
				MatrixNonzeros(&d_nnz,&o_nnz,femmodel,FsetEnum,FsetEnum);
				Kff=new Matrix<IssmDouble>(m,n,M,N,d_nnz,o_nnz);
				xDelete<int>(d_nnz);
				xDelete<int>(o_nnz);
			}
		}
		if(pKfs){
			m=flocalsize; n=slocalsize; /*local  sizes*/
			M=fsize;      N=ssize;      /*global sizes*/
			if(strcmp(toolkittype,"issm")==0){
				Kfs=new Matrix<IssmDouble>(m,n,M,N,NULL,NULL);
			}
			else{
				MatrixNonzeros(&d_nnz,&o_nnz,femmodel,FsetEnum,SsetEnum);
				Kfs=new Matrix<IssmDouble>(m,n,M,N,d_nnz,o_nnz);
				xDelete<int>(d_nnz);
				xDelete<int>(o_nnz);
			}
		}
		if(pdf) df =new Vector<IssmDouble>(flocalsize,fsize);
		if(ppf) pf =new Vector<IssmDouble>(flocalsize,fsize);
	}
	
	/*Free ressources: */
	xDelete<char>(toolkittype);

	/*Allocate output pointers*/
	if(pKff) *pKff = Kff;
	if(pKfs) *pKfs = Kfs;
	if(pdf)  *pdf  = df;
	if(ppf)  *ppf  = pf;
}

void MatrixNonzeros(int** pd_nnz,int** po_nnz,FemModel* femmodel,int set1enum,int set2enum){

	/*Intermediary*/
	int      i,j,k,index,offset,count;
	int      configuration_type;
	int      d_nz,o_nz;
	Element *element            = NULL;
	Load    *load               = NULL;
	int     *head_e             = NULL;
	int     *next_e             = NULL;
	int     *count2offset_e     = NULL;
	int     *head_l             = NULL;
	int     *next_l             = NULL;
	int     *count2offset_l     = NULL;
	int     *lidlist            = NULL;

	/*output*/
	int *d_nnz = NULL;
	int *o_nnz = NULL;

	/*retrive parameters: */
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	/*Get vector size and number of nodes*/
	int numnodes            = femmodel->nodes->NumberOfNodes(configuration_type);
	int localnumnodes       = femmodel->nodes->Size();
	int numberofdofspernode = femmodel->nodes->MaxNumDofs(configuration_type,GsetEnum);
	int M                   = femmodel->nodes->NumberOfDofs(configuration_type,set1enum);
	int N                   = femmodel->nodes->NumberOfDofs(configuration_type,set2enum);
	int m                   = femmodel->nodes->NumberOfDofsLocal(configuration_type,set1enum);
	int n                   = femmodel->nodes->NumberOfDofsLocal(configuration_type,set2enum);
	int numnodesperelement  = femmodel->elements->MaxNumNodes();
	int numnodesperload     = femmodel->loads->MaxNumNodes(configuration_type);

	/*First, we are building chaining vectors so that we know what nodes are
	 * connected to what elements. These vectors are such that:
	 *   for(int i=head[id];i!=-1;i=next[i])
	 * will loop over all the elements that are connected to the node number
	 * id*/
	head_e         = xNew<int>(localnumnodes); for(i=0;i<localnumnodes;i++) head_e[i]=-1;
	next_e         = xNew<int>(femmodel->elements->Size()*numnodesperelement);
	count2offset_e = xNew<int>(femmodel->elements->Size()*numnodesperelement);

	k=0;
	for(i=0;i<femmodel->elements->Size();i++){
		element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		lidlist = xNew<int>(element->GetNumberOfNodes());
		element->GetNodesLidList(lidlist);

		for(j=0;j<element->GetNumberOfNodes();j++){
			index = lidlist[j];
			_assert_(index>=0 && index<numnodes);

			count2offset_e[k]=i;
			next_e[k]=head_e[index];
			head_e[index]=k++;
		}
		for(j=0;j<numnodesperelement-element->GetNumberOfNodes();j++) k++;

		xDelete<int>(lidlist);
	}

	/*Chain for loads*/
	head_l         = xNew<int>(localnumnodes); for(i=0;i<localnumnodes;i++) head_l[i]=-1;
	next_l         = xNew<int>(femmodel->loads->Size(configuration_type)*numnodesperload);
	count2offset_l = xNew<int>(femmodel->loads->Size(configuration_type)*numnodesperload);
	k=0;
	for(i=0;i<femmodel->loads->Size();i++){
		load = xDynamicCast<Load*>(femmodel->loads->GetObjectByOffset(i));
		if(!load->InAnalysis(configuration_type)) continue;
		lidlist = xNew<int>(load->GetNumberOfNodes());
		load->GetNodesLidList(lidlist);

		for(j=0;j<load->GetNumberOfNodes();j++){
			index = lidlist[j];
			_assert_(index>=0 && index<numnodes);

			count2offset_l[k]=i;
			next_l[k]=head_l[index];
			head_l[index]=k++;
		}
		for(j=0;j<numnodesperload-load->GetNumberOfNodes();j++) k++;

		xDelete<int>(lidlist);
	}

	/*OK now count number of dofs and flag each nodes for each node i*/
	bool *flags                  = xNew<bool>(localnumnodes);
	int  *flagsindices           = xNew<int>(localnumnodes);
	int  *d_connectivity         = xNewZeroInit<int>(numnodes);
	int  *o_connectivity         = xNewZeroInit<int>(numnodes);
	int  *connectivity_clone     = xNewZeroInit<int>(numnodes);
	int  *all_connectivity_clone = xNewZeroInit<int>(numnodes);

	/*Resetting flags to false at eahc iteration takes a lot of time, so we keep track of the flags
	 * to reset in flagsindices, initialized with -1*/
	for(i = 0;i<localnumnodes;i++) flags[i]        = false;
	for(i = 0;i<localnumnodes;i++) flagsindices[i] = -1;

	/*Create connectivity vector*/
	for(i=0;i<femmodel->nodes->Size();i++){
		Node* node=xDynamicCast<Node*>(femmodel->nodes->GetObjectByOffset(i));
		if(node->InAnalysis(configuration_type)){

			/*Reinitialize flags to false*/
			j=0;
			while(true){
				if(flagsindices[j]>=0){
					flags[flagsindices[j]] = false;
					flagsindices[j]        = -1;
					j++;
				}
				else{
					break;
				}
			}

			//for(j=0;j<localnumnodes;j++) flags[j]=false;

			/*Loop over elements that hold node number i*/
			//if(head_e[node->Lid()]==-1 && head_l[node->Lid()]==-1){
			//	printf("[%i] vertex %i\n",IssmComm::GetRank(),node->Lid()+1);
			//}
			for(j=head_e[node->Lid()];j!=-1;j=next_e[j]){
				offset=count2offset_e[j];
				element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(offset));
				element->SetwiseNodeConnectivity(&d_nz,&o_nz,node,flags,flagsindices,set1enum,set2enum);
				if(node->IsClone()){
					connectivity_clone[node->Sid()]+=d_nz+o_nz;
				}
				else{
					d_connectivity[node->Sid()]+=d_nz;
					o_connectivity[node->Sid()]+=o_nz;
				}
			}
			for(j=head_l[node->Lid()];j!=-1;j=next_l[j]){
				offset=count2offset_l[j];
				load=xDynamicCast<Load*>(femmodel->loads->GetObjectByOffset(offset));
				load->SetwiseNodeConnectivity(&d_nz,&o_nz,node,flags,flagsindices,set1enum,set2enum);
				if(node->IsClone()){
					connectivity_clone[node->Sid()]+=d_nz+o_nz;
				}
				else{
					d_connectivity[node->Sid()]+=d_nz;
					o_connectivity[node->Sid()]+=o_nz;
				}
			}
		}
	}
	xDelete<bool>(flags);
	xDelete<int>(flagsindices);
	xDelete<int>(count2offset_e);
	xDelete<int>(head_e);
	xDelete<int>(next_e);
	xDelete<int>(count2offset_l);
	xDelete<int>(head_l);
	xDelete<int>(next_l);

	/*sum over all cpus*/
	ISSM_MPI_Allreduce((void*)connectivity_clone,(void*)all_connectivity_clone,numnodes,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	xDelete<int>(connectivity_clone);

	if(set1enum==FsetEnum){
		count=0;
		d_nnz=xNew<int>(m);
		o_nnz=xNew<int>(m);
		for(i=0;i<femmodel->nodes->Size();i++){
			Node* node=xDynamicCast<Node*>(femmodel->nodes->GetObjectByOffset(i));
			if(node->InAnalysis(configuration_type) && !node->IsClone()){
				for(j=0;j<node->indexing.fsize;j++){
					_assert_(count<m);
					d_nnz[count]=numberofdofspernode*(d_connectivity[node->Sid()] + all_connectivity_clone[node->Sid()]);
					o_nnz[count]=numberofdofspernode*(o_connectivity[node->Sid()] + all_connectivity_clone[node->Sid()]);
					if(d_nnz[count]>n)   d_nnz[count]=n;
					if(o_nnz[count]>N-n) o_nnz[count]=N-n;
					count++;
				}
			}
		}
		_assert_(m==count);
	}
	else{
		_error_("STOP not implemented");
	}
	xDelete<int>(d_connectivity);
	xDelete<int>(o_connectivity);
	xDelete<int>(all_connectivity_clone);

	/*Allocate ouptput pointer*/
	*pd_nnz=d_nnz;
	*po_nnz=o_nnz;
}
