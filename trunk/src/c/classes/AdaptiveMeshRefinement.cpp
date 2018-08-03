/*!\file AdaptiveMeshrefinement.cpp
 * \brief: implementation of the adaptive mesh refinement tool based on NeoPZ library: github.com/labmec/neopz
 */

#ifdef HAVE_CONFIG_H
    #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./AdaptiveMeshRefinement.h"

/*Constructor, copy, clean up and destructor*/
AdaptiveMeshRefinement::AdaptiveMeshRefinement(){/*{{{*/
	this->Initialize();
}
/*}}}*/
AdaptiveMeshRefinement::AdaptiveMeshRefinement(const AdaptiveMeshRefinement &cp){/*{{{*/
   this->Initialize(); 
	this->operator =(cp);
}
/*}}}*/
AdaptiveMeshRefinement & AdaptiveMeshRefinement::operator =(const AdaptiveMeshRefinement &cp){/*{{{*/

	/*Clean all attributes*/
	this->CleanUp();

	/*Copy all data*/
	this->fathermesh     = new TPZGeoMesh(*cp.fathermesh);
	this->previousmesh   = new TPZGeoMesh(*cp.previousmesh); 
	this->levelmax       = cp.levelmax;
	this->elementswidth  = cp.elementswidth;
	this->regionlevel1   = cp.regionlevel1;
	this->regionlevelmax = cp.regionlevelmax;
	return *this;

}
/*}}}*/
AdaptiveMeshRefinement::~AdaptiveMeshRefinement(){/*{{{*/
	
	bool ismismip = true;
	if(ismismip){//itapopo
		TPZFileStream fstr;
		std::stringstream ss;
	    
		ss << this->levelmax;
		std::string AMRfile	= "/home/santos/Misomip2/L" + ss.str() + "_tsai/amr.txt"; 
	
		fstr.OpenWrite(AMRfile.c_str());
		int withclassid = 1;
		this->Write(fstr,withclassid);
	}
	this->CleanUp();
	gRefDBase.clear();
}
/*}}}*/
void AdaptiveMeshRefinement::CleanUp(){/*{{{*/

    /*Verify and delete all data*/
	if(this->fathermesh)    delete this->fathermesh;
   if(this->previousmesh)  delete this->previousmesh;
	this->levelmax			= -1;
	this->elementswidth  = -1;
	this->regionlevel1	= -1;
	this->regionlevelmax = -1;
}
/*}}}*/
void AdaptiveMeshRefinement::Initialize(){/*{{{*/

	/*Set pointers to NULL*/
	this->fathermesh		= NULL;
	this->previousmesh	= NULL;
	this->levelmax			= -1;
	this->elementswidth	= -1;
	this->regionlevel1	= -1;
	this->regionlevelmax = -1;
}
/*}}}*/
int AdaptiveMeshRefinement::ClassId() const{/*{{{*/
    return 13829430; //Antartic area with ice shelves (km^2)
}
/*}}}*/
void AdaptiveMeshRefinement::Read(TPZStream &buf, void *context){/*{{{*/

    try
    {
        /* Read the id context*/
        TPZSaveable::Read(buf,context);

        /* Read class id*/
        int classid;
        buf.Read(&classid,1);
        
        /* Verify the class id*/
        if (classid != this->ClassId() )
        {
            std::cout << "Error in restoring AdaptiveMeshRefinement!\n";
            std::cout.flush();
            DebugStop();
        }
        
        /* Read simple attributes */
        buf.Read(&this->levelmax,1);
        buf.Read(&this->elementswidth,1);
        buf.Read(&this->regionlevel1,1);
        buf.Read(&this->regionlevelmax,1);
        
		/* Read geometric mesh*/
        TPZSaveable *sv1 = TPZSaveable::Restore(buf,0);
        this->fathermesh = dynamic_cast<TPZGeoMesh*>(sv1);
        
        TPZSaveable *sv2 = TPZSaveable::Restore(buf,0);
        this->previousmesh = dynamic_cast<TPZGeoMesh*>(sv2);
    }
    catch(const std::exception& e)
    {
        std::cout << "Exception catched! " << e.what() << std::endl;
        std::cout.flush();
        DebugStop();
    }
}
/*}}}*/
template class TPZRestoreClass<AdaptiveMeshRefinement,13829430>;/*{{{*/
/*}}}*/
void AdaptiveMeshRefinement::Write(TPZStream &buf, int withclassid){/*{{{*/
    
    try
    {
        /* Write context (this class) class ID*/
        TPZSaveable::Write(buf,withclassid);

        /* Write this class id*/
        int classid = ClassId();
        buf.Write(&classid,1);

        /* Write simple attributes */
        buf.Write(&this->levelmax,1);
        buf.Write(&this->elementswidth,1);
        buf.Write(&this->regionlevel1,1);
        buf.Write(&this->regionlevelmax,1);
			
        /* Write the geometric mesh*/
        this->fathermesh->Write(buf, this->ClassId());
        this->previousmesh->Write(buf, this->ClassId());
    }
    catch(const std::exception& e)
    {
        std::cout << "Exception catched! " << e.what() << std::endl;
        std::cout.flush();
        DebugStop();
    }
}
/*}}}*/

/*Mesh refinement methods*/
#include "TPZVTKGeoMesh.h" //itapopo
#include "../shared/shared.h" //itapopo
void AdaptiveMeshRefinement::ExecuteRefinement(int &type_process,double *vx, double *vy, double *masklevelset, int &nvertices, int &nelements, int &nsegments, double** px, double** py, double** pz, int** pelements, int** psegments){/*{{{*/

	/*IMPORTANT! pelements (and psegments) are in Matlab indexing*/
	/*NEOPZ works only in C indexing*/

    _assert_(this->fathermesh);
    _assert_(this->previousmesh);
    
    /*Calculate the position of the grounding line using previous mesh*/
    std::vector<TPZVec<REAL> > GLvec;
    this->CalcGroundingLinePosition(masklevelset, GLvec);
    
   // std::ofstream file1("/home/santos/mesh0.vtk");
   // TPZVTKGeoMesh::PrintGMeshVTK(this->fathermesh,file1 );
    
    /*run refinement or unrefinement process*/
    TPZGeoMesh *newmesh;
    switch (type_process) {
        case 0: newmesh = this->previousmesh; break;                    // refine previous mesh
        case 1: newmesh = new TPZGeoMesh(*this->fathermesh); break;     // refine mesh 0 (unrefine process)
        default: DebugStop(); break;//itapopo verificar se irá usar _assert_
    }
    
    this->RefinementProcess(newmesh,GLvec);
	
    //std::ofstream file2("/home/santos/mesh1.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(this->previousmesh,file2 );
    
    /*Set new mesh pointer. Previous mesh just have uniform elements*/
    if(type_process==1){
        if(this->previousmesh) delete this->previousmesh;
        this->previousmesh = newmesh;
    }
    
    /*Refine elements to avoid hanging nodes*/
	//TPZGeoMesh *nohangingnodesmesh = new TPZGeoMesh(*newmesh);//itapopo testando, este era o original
   TPZGeoMesh *nohangingnodesmesh = this->CreateRefPatternMesh(newmesh);//itapopo testando, este eh novo metodo
    
    //std::ofstream file3("/home/santos/mesh2.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(this->previousmesh,file3);
    
    this->RefineMeshToAvoidHangingNodes(nohangingnodesmesh);
    
	 //std::ofstream file4("/home/santos/mesh3.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(nohangingnodesmesh,file4);
    
    /*Get new geometric mesh in ISSM data structure*/
    this->GetMesh(nohangingnodesmesh,nvertices,nelements,nsegments,px,py,pz,pelements,psegments);
	 
    /*Verify the new geometry*/
    this->CheckMesh(nvertices,nelements,nsegments,this->elementswidth,px,py,pz,pelements,psegments);

	 _printf_("\trefinement process done!\n\n");

    delete nohangingnodesmesh;
}
/*}}}*/
void AdaptiveMeshRefinement::RefinementProcess(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec){/*{{{*/
    
    /*Refine mesh levelmax times*/
   _printf_("\n\trefinement process (level max = " << this->levelmax << ")\n");
	_printf_("\tprogress:  ");
	for(int hlevel=1;hlevel<=this->levelmax;hlevel++){
        
        /*Set elements to be refined using some criteria*/
        std::vector<int> ElemVec; //elements without children
        this->SetElementsToRefine(gmesh,GLvec,hlevel,ElemVec);
        
        /*Refine the mesh*/
        this->RefineMesh(gmesh, ElemVec);
		  
		  _printf_("*  ");
    }
    _printf_("\n");
}
/*}}}*/
void AdaptiveMeshRefinement::RefineMesh(TPZGeoMesh *gmesh, std::vector<int> &ElemVec){/*{{{*/

	/*Refine elements in ElemVec: uniform pattern refinement*/
	for(int i = 0; i < ElemVec.size(); i++){
		
        /*Get geometric element and verify if it has already been refined*/
        int index = ElemVec[i];
        TPZGeoEl * geoel = gmesh->Element(index);
        if(geoel->HasSubElement()) DebugStop();                              //itapopo _assert_(!geoel->HasSubElement());
        if(geoel->MaterialId() != this->GetElemMaterialID()) DebugStop();   //itapopo verificar se usará _assert_
        
        /*Divide geoel*/
        TPZVec<TPZGeoEl *> Sons;
		  geoel->Divide(Sons);
        
        /*If a 1D segment is neighbor, it must be divided too*/
        if(this->elementswidth != 3) DebugStop(); //itapopo verificar o segment para malha 3D
        
        std::vector<int> sides(3);
        sides[0] = 3; sides[1] = 4; sides[2] = 5;
        for(int j = 0; j < sides.size(); j++ ){
            
            TPZGeoElSide Neighbour = geoel->Neighbour(sides[j]);
            
            if( Neighbour.Element()->MaterialId() == this->GetBoundaryMaterialID() && !Neighbour.Element()->HasSubElement() ){
                TPZVec<TPZGeoEl *> pv2;
                Neighbour.Element()->Divide(pv2);
            }
        }
	}
    
    gmesh->BuildConnectivity();

}
/*}}}*/
void AdaptiveMeshRefinement::RefineMeshToAvoidHangingNodes(TPZGeoMesh *gmesh){/*{{{*/
   
	 _printf_("\trefine to avoid hanging nodes...\n");
    /*Refine elements to avoid hanging nodes: non-uniform refinement*/
	 const int NElem = gmesh->NElements();
    for(int i = 0; i < NElem; i++){
        
        /*Get geometric element and verify if it has already been refined. Geoel may not have been previously refined*/
        TPZGeoEl * geoel = gmesh->Element(i);
        if(!geoel) continue;
        if(geoel->HasSubElement()) continue;
        if(geoel->MaterialId() != this->GetElemMaterialID()) continue;
        
        /*Get the refinement pattern for this element and refine it*/
        TPZAutoPointer<TPZRefPattern> refp = TPZRefPatternTools::PerfectMatchRefPattern(geoel);
        if(refp){
            TPZVec<TPZGeoEl *> Sons;
            geoel->SetRefPattern(refp);
            geoel->Divide(Sons);
        }
        
    }
    
    gmesh->BuildConnectivity();
    
}
/*}}}*/
void AdaptiveMeshRefinement::GetMesh(TPZGeoMesh *gmesh,int &nvertices,int &nelements,int &nsegments,double** px,double** py,double** pz, int** pelements, int** psegments){/*{{{*/

	/*IMPORTANT! pelements (and psegments) are in Matlab indexing*/
	/*NEOPZ works only in C indexing*/

	/* vertices */
    int ntotalvertices = gmesh->NNodes();//total
    
    /* mesh coords */
	 double* newmeshX = xNew<IssmDouble>(ntotalvertices);
    double* newmeshY = xNew<IssmDouble>(ntotalvertices);
    double* newmeshZ = xNew<IssmDouble>(ntotalvertices);
   	
   /* getting mesh coords */
    for(int i = 0; i < ntotalvertices; i++ ){
        TPZVec<REAL> coords(3,0.);
        gmesh->NodeVec()[i].GetCoordinates(coords);
        newmeshX[i] = coords[0];
        newmeshY[i] = coords[1];
        newmeshZ[i] = coords[2];
    }
    
	/* elements */
    std::vector<TPZGeoEl*> GeoVec; GeoVec.clear();
    for(int i = 0; i < gmesh->NElements(); i++){ 
        if( gmesh->ElementVec()[i]->HasSubElement() ) continue;
        if( gmesh->ElementVec()[i]->MaterialId() != this->GetElemMaterialID() ) continue;
        GeoVec.push_back( gmesh->ElementVec()[i]);
    }
    
    int ntotalelements = (int)GeoVec.size();
    int* newelements = xNew<int>(ntotalelements*this->elementswidth);

    if ( !(this->elementswidth == 3) && !(this->elementswidth == 4) && !(this->elementswidth == 6) ) DebugStop();

    for(int i=0;i<GeoVec.size();i++){
        for(int j=0;j<this->elementswidth;j++) newelements[i*this->elementswidth+j]=(int)GeoVec[i]->NodeIndex(j)+1;//C to Matlab indexing
	 }
    
    /* segments */
    std::vector<TPZGeoEl*> SegVec; SegVec.clear();
    for(int i = 0; i < gmesh->NElements(); i++){
        if( gmesh->ElementVec()[i]->HasSubElement() ) continue;
        if( gmesh->ElementVec()[i]->MaterialId() != this->GetBoundaryMaterialID() ) continue;
        SegVec.push_back( gmesh->ElementVec()[i]);
    }
    
    int ntotalsegments = (int)SegVec.size();
    int *newsegments=NULL;
	 if(ntotalsegments>0) newsegments=xNew<int>(ntotalsegments*3);
    
    for(int i=0;i<SegVec.size();i++){
        
        for(int j=0;j<2;j++) newsegments[i*3+j]=(int)SegVec[i]->NodeIndex(j)+1;//C to Matlab indexing
        
        int neighborindex = SegVec[i]->Neighbour(2).Element()->Index();
        int neighbourid = -1;
        
        for(int j = 0; j < GeoVec.size(); j++){
            if( GeoVec[j]->Index() == neighborindex || GeoVec[j]->FatherIndex() == neighborindex){
                neighbourid = j;
                break;
            }
        }
        
        if(neighbourid==-1) DebugStop(); //itapopo talvez passar para _assert_
        newsegments[i*3+2] = neighbourid+1;//C to Matlab indexing
    }
    
    //setting outputs
    nvertices  = ntotalvertices;
    nelements  = ntotalelements;
    nsegments  = ntotalsegments;
    *px		   = newmeshX;
    *py		   = newmeshY;
    *pz		   = newmeshZ;
    *pelements = newelements;
    *psegments = newsegments;
    
}
/*}}}*/
void AdaptiveMeshRefinement::CalcGroundingLinePosition(double *masklevelset,std::vector<TPZVec<REAL> > &GLvec){/*{{{*/
    
    /* Find grounding line using elments center point */
    GLvec.clear();
    for(int i=0;i<this->previousmesh->NElements();i++){
        
        if(this->previousmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
        if(this->previousmesh->Element(i)->HasSubElement()) continue;
        
        //itapopo apenas malha 2D triangular!
        int vertex0 = this->previousmesh->Element(i)->NodeIndex(0);
        int vertex1 = this->previousmesh->Element(i)->NodeIndex(1);
        int vertex2 = this->previousmesh->Element(i)->NodeIndex(2);
        
		  //itapopo inserir uma verificação para não acessar fora da memória
        double mls0 = masklevelset[vertex0];
        double mls1 = masklevelset[vertex1];
        double mls2 = masklevelset[vertex2];
        
        if( mls0*mls1 < 0. || mls1*mls2 < 0. ){
            const int side = 6;
            TPZVec<double> qsi(2,0.);
            TPZVec<double> X(3,0.);
            this->previousmesh->Element(i)->CenterPoint(side, qsi);
            this->previousmesh->Element(i)->X(qsi, X);
            GLvec.push_back(X);
        }
    }
    
//    itapopo apenas para debugar
//    std::ofstream fileGL("/Users/santos/Desktop/gl.nb");
//    fileGL << "ListPlot[{";
//    for(int i = 0; i < GLvec.size(); i++){
//        fileGL << "{" << GLvec[i][0] << "," << GLvec[i][1] << /*"," << 0. << */"}";
//        if(i != GLvec.size()-1) fileGL << ",";
//    }
//    fileGL << "}]";
//    fileGL.flush();
//    fileGL.close();
    
}
/*}}}*/
void AdaptiveMeshRefinement::SetElementsToRefine(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec,int &hlevel,std::vector<int> &ElemVec){/*{{{*/

    if(!gmesh) DebugStop(); //itapopo verificar se usará _assert_

    ElemVec.clear();
    
    // itapopo inserir modo de encontrar criterio TESTING!!!! Come to false!
	 if(false) this->TagAllElements(gmesh,ElemVec); //uniform, refine all elements!

    /* Adaptive refinement. This refines some elements following some criteria*/
    this->TagElementsNearGroundingLine(gmesh, GLvec, hlevel, ElemVec);

}
/*}}}*/
void AdaptiveMeshRefinement::TagAllElements(TPZGeoMesh *gmesh,std::vector<int> &ElemVec){/*{{{*/
    
    /* Uniform refinement. This refines the entire mesh */
    int nelements = gmesh->NElements();
    for(int i=0;i<nelements;i++){
        if(gmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
        if(gmesh->Element(i)->HasSubElement()) continue;
        ElemVec.push_back(i);
    }
}
/*}}}*/
void AdaptiveMeshRefinement::TagElementsNearGroundingLine(TPZGeoMesh *gmesh,std::vector<TPZVec<REAL> > &GLvec,int &hlevel,std::vector<int> &ElemVec){/*{{{*/
    
    /* Tag elements near grounding line */ 
	 double D1		= this->regionlevel1;
	 double Dhmax  = this->regionlevelmax;
	 int hmax		= this->levelmax;
    double alpha	= (hmax==1) ? 0. : log(D1/Dhmax)/(hmax-1.);
	 double Di		= D1/exp(alpha*(hlevel-1));
    
    for(int i=0;i<gmesh->NElements();i++){
        
        if(gmesh->Element(i)->MaterialId()!=this->GetElemMaterialID()) continue;
        if(gmesh->Element(i)->HasSubElement()) continue;
        if(gmesh->Element(i)->Level()>=hlevel) continue;
        
        const int side2D = 6;
        TPZVec<REAL> qsi(2,0.);
        TPZVec<REAL> centerPoint(3,0.);
        gmesh->Element(i)->CenterPoint(side2D, qsi);
        gmesh->Element(i)->X(qsi, centerPoint);
        
        REAL distance = Di;
        
        for (int j = 0; j < GLvec.size(); j++) {
            
            REAL value = ( GLvec[j][0] - centerPoint[0] ) * ( GLvec[j][0] - centerPoint[0] ); // (x2-x1)^2
            value += ( GLvec[j][1] - centerPoint[1] ) * ( GLvec[j][1] - centerPoint[1] );// (y2-y1)^2
            value = std::sqrt(value); //Radius
            
            //finding the min distance to the grounding line
            if(value < distance) distance = value;
            
        }
        
        if(distance < Di) ElemVec.push_back(i);
    }
    
}
/*}}}*/
void AdaptiveMeshRefinement::CreateInitialMesh(int &nvertices,int &nelements,int &nsegments,int &width,double* x,double* y,double* z,int* elements,int* segments){/*{{{*/

	/*IMPORTANT! elements come in Matlab indexing*/
	/*NEOPZ works only in C indexing*/
	
	_assert_(nvertices>0);
   _assert_(nelements>0);
	this->SetElementWidth(width);

    /*Verify and creating initial mesh*/
   if(this->fathermesh) _error_("Initial mesh already exists!");
    
   this->fathermesh = new TPZGeoMesh();
	this->fathermesh->NodeVec().Resize( nvertices );

    /*Set the vertices (geometric nodes in NeoPZ context)*/
	for(int i=0;i<nvertices;i++){  
      /*x,y,z coords*/
		TPZManVector<REAL,3> coord(3,0.);
      coord[0]= x[i];
      coord[1]= y[i];
      coord[2]= z[i];
      /*Insert in the mesh*/
      this->fathermesh->NodeVec()[i].SetCoord(coord);
		this->fathermesh->NodeVec()[i].SetNodeId(i);
	}
	
	/*Generate the elements*/
	long index;
   const int mat = this->GetElemMaterialID();
   TPZManVector<long> elem(this->elementswidth,0);
    
	for(int iel=0;iel<nelements;iel++){

		for(int jel=0;jel<this->elementswidth;jel++) elem[jel]=elements[iel*this->elementswidth+jel]-1;//Convert Matlab to C indexing

      /*reftype = 0: uniform, fast / reftype = 1: uniform and non-uniform (avoid hanging nodes), it is not too fast */
      const int reftype = 0;
      switch(this->elementswidth){
			case 3: this->fathermesh->CreateGeoElement(ETriangle, elem, mat, index, reftype);	break;
         case 4: this->fathermesh->CreateGeoElement(ETetraedro, elem, mat, index, reftype); DebugStop(); break;
			case 6: this->fathermesh->CreateGeoElement(EPrisma, elem, mat, index, reftype); DebugStop(); break;
         default:	DebugStop();//itapopo _error_("mesh not supported yet");
		}
        
      /*Define the element ID*/        
      this->fathermesh->ElementVec()[index]->SetId(iel); 
	}
    
   /*Generate the 1D segments elements (boundary)*/
   const int matboundary = this->GetBoundaryMaterialID();
   TPZManVector<long> boundary(2,0.);
    
   for(int iel=nelements;iel<nelements+nsegments;iel++){     
		boundary[0] = segments[(iel-nelements)*2+0]-1;//Convert Matlab to C indexing
      boundary[1] = segments[(iel-nelements)*2+1]-1;//Convert Matlab to C indexing
      /*reftype = 0: uniform, fast / reftype = 1: uniform and non-uniform (avoid hanging nodes), it is not too fast */
      const int reftype = 0;
      this->fathermesh->CreateGeoElement(EOned, boundary, matboundary, index, reftype);//cria elemento unidimensional
      this->fathermesh->ElementVec()[index]->SetId(iel);
	}
    
   /*Build element and node connectivities*/
   this->fathermesh->BuildConnectivity();
    
	/*Create previous mesh as a copy of father mesh*/
   this->previousmesh = new TPZGeoMesh(*this->fathermesh);

}
/*}}}*/
#include "pzgeotriangle.h" //itapopo
#include "pzreftriangle.h" //itapopo
using namespace pzgeom;
TPZGeoMesh* AdaptiveMeshRefinement::CreateRefPatternMesh(TPZGeoMesh* gmesh){/*{{{*/
	
	TPZGeoMesh *newgmesh = new TPZGeoMesh();
   newgmesh->CleanUp();
    
   int nnodes  = gmesh->NNodes();
	int nelem   = gmesh->NElements();
   int mat     = this->GetElemMaterialID();;
   int reftype = 1;
   long index; 
   
	//nodes
	newgmesh->NodeVec().Resize(nnodes);
   for(int i=0;i<nnodes;i++) newgmesh->NodeVec()[i] = gmesh->NodeVec()[i];
    
   //elements
   for(int i=0;i<nelem;i++){
   	TPZGeoEl * geoel = gmesh->Element(i);
      TPZManVector<long> elem(3,0);
      for(int j=0;j<3;j++) elem[j] = geoel->NodeIndex(j);
     
      newgmesh->CreateGeoElement(ETriangle,elem,mat,index,reftype);
      newgmesh->ElementVec()[index]->SetId(geoel->Id());
        
      TPZGeoElRefPattern<TPZGeoTriangle>* newgeoel = dynamic_cast<TPZGeoElRefPattern<TPZGeoTriangle>*>(newgmesh->ElementVec()[index]);
        
      //old neighbourhood
      const int nsides = TPZGeoTriangle::NSides;
      TPZVec< std::vector<TPZGeoElSide> > neighbourhood(nsides);
      TPZVec<long> NodesSequence(0);
      for(int s = 0; s < nsides; s++){
      	neighbourhood[s].resize(0);
      	TPZGeoElSide mySide(geoel,s);
      	TPZGeoElSide neighS = mySide.Neighbour();
         if(mySide.Dimension() == 0){
         	long oldSz = NodesSequence.NElements();
            NodesSequence.resize(oldSz+1);
            NodesSequence[oldSz] = geoel->NodeIndex(s);
         }
      	while(mySide != neighS){
         	neighbourhood[s].push_back(neighS);
            neighS = neighS.Neighbour();
         }
      }
        
      //inserting in new element
      for(int s = 0; s < nsides; s++){
      	TPZGeoEl * tempEl = newgeoel;
         TPZGeoElSide tempSide(newgeoel,s);
         int byside = s;
         for(unsigned long n = 0; n < neighbourhood[s].size(); n++){
         	TPZGeoElSide neighS = neighbourhood[s][n];
            tempEl->SetNeighbour(byside, neighS);
            tempEl = neighS.Element();
            byside = neighS.Side();
         }
         tempEl->SetNeighbour(byside, tempSide);
      }
        
      long fatherindex = geoel->FatherIndex();
      if(fatherindex>-1) newgeoel->SetFather(fatherindex);
        
      if(!geoel->HasSubElement()) continue;
        
      int nsons = geoel->NSubElements();

      TPZAutoPointer<TPZRefPattern> ref = gRefDBase.GetUniformRefPattern(ETriangle);
      newgeoel->SetRefPattern(ref);
        
      for(int j=0;j<nsons;j++){
      	TPZGeoEl* son = geoel->SubElement(j);
         if(!son){
             DebugStop();
         }
         newgeoel->SetSubElement(j,son);
      }
   }
	newgmesh->BuildConnectivity();
    
	return newgmesh;
}
/*}}}*/
void AdaptiveMeshRefinement::SetLevelMax(int &h){/*{{{*/
    this->levelmax = h;
}
/*}}}*/
void AdaptiveMeshRefinement::SetRegions(double &D1,double Dhmax){/*{{{*/
    this->regionlevel1	 = D1;
    this->regionlevelmax = Dhmax;
}
/*}}}*/
void AdaptiveMeshRefinement::SetElementWidth(int &width){/*{{{*/
    this->elementswidth = width;
}
/*}}}*/
void AdaptiveMeshRefinement::CheckMesh(int &nvertices,int &nelements,int &nsegments,int &width,double** px,double** py,double** pz,int** pelements, int** psegments){/*{{{*/

    /*Basic verification*/
    if( !(nvertices > 0) || !(nelements > 0) ) DebugStop(); //itapopo verificar se irá usar o _assert_
    
    if ( !(width == 3) && !(width == 4) && !(width == 6) ) DebugStop(); // itapopo verifcar se irá usar o _assert_
    
    if( !px || !py || !pz || !pelements ) DebugStop(); // itapopo verifcar se irá usar o _assert_
    
    /*Verify if there are orphan nodes*/
    std::set<int> elemvertices;
    elemvertices.clear(); 
    for(int i = 0; i < nelements; i++){
        for(int j = 0; j < width; j++) {
            elemvertices.insert((*pelements)[i*width+j]);
		  }
	 }
    
    if( elemvertices.size() != nvertices ) DebugStop();//itapopo verificar se irá usar o _assert_
	
    //Verify if there are inf or NaN in coords
    for(int i = 0; i < nvertices; i++) if(isnan((*px)[i]) || isinf((*px)[i])) DebugStop();
    for(int i = 0; i < nvertices; i++) if(isnan((*py)[i]) || isinf((*py)[i])) DebugStop();
    for(int i = 0; i < nvertices; i++) if(isnan((*pz)[i]) || isinf((*pz)[i])) DebugStop();
   
	 for(int i = 0; i < nelements; i++){
        for(int j = 0; j < width; j++){
            if( isnan((*pelements)[i*width+j]) || isinf((*pelements)[i*width+j]) ) DebugStop();
        }
    }
    
}
/*}}}*/
