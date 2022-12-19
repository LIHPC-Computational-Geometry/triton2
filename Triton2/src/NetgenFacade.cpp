/*----------------------------------------------------------------------------*/
/** \file    NetgenFacade.t.h
 *  \author  F. LEDOUX
 *  \date    16/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <Triton2/NetgenInterface/NetgenFacade.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace triton;
/*----------------------------------------------------------------------------*/
NetgenFacade::NetgenFacade():ng_mesh_(0)
{}
/*----------------------------------------------------------------------------*/

void NetgenFacade::initializeNetgen(NetgenParams& AParams){
	// netgen init
	NgAddOn_Init();

	// surface parameters
	ng_mesh_param_surf_.maxh =AParams.maxh;
	ng_mesh_param_surf_.fineness = AParams.fineness;
	ng_mesh_param_surf_.secondorder=AParams.secondorder;
	ng_mesh_param_surf_.uselocalh=AParams.uselocalh;
	ng_mesh_param_surf_.grading=AParams.grading;
	ng_mesh_param_surf_.closeedgeenable = AParams.closeedgeenable;
	ng_mesh_param_surf_.elementsperedge=AParams.elementsperedge;
	ng_mesh_param_surf_.elementspercurve=AParams.elementspercurve;

	// volume parameters
	ng_mesh_param_vol_.maxh =AParams.maxh;
	ng_mesh_param_vol_.fineness = AParams.fineness;
	ng_mesh_param_vol_.secondorder= AParams.secondorder;


	// creation of the netgen mesh and the stl geometry to mesh
	ng_mesh_ = Ng_NewMesh();
//	stl_geometry_ = Ng_STL_NewGeometry();

}
/*----------------------------------------------------------------------------*/

void NetgenFacade::initializeNetgen(const double AQuality,
		  gmds::math::Point& APMin,
		  gmds::math::Point& APMax,
		  const double ASize){

	// TODO legoff : look into meaning of quality

	if(AQuality<1.0 ||AQuality>2.0)
		throw gmds::GMDSException("The quality parameter must be included in [1.0,2.0]");


	// TODO legoff : look into relationship between size max for surface triangles and size
	// max for volume tetrahedron.
	// Also look into max tet size and desired tet size in the box.

	// surface parameters
	ng_mesh_param_surf_.maxh =ASize;
//	ng_mesh_param_surf_.maxh =1.;
	ng_mesh_param_surf_.fineness = 0.5;
	ng_mesh_param_surf_.secondorder=0;

//	ng_mesh_param_surf_.uselocalh=1;
//	ng_mesh_param_surf_.grading=1.;

	// volume parameters
	ng_mesh_param_vol_.maxh =ASize*100.;
//	ng_mesh_param_vol_.maxh =2.;
	ng_mesh_param_vol_.fineness = 0.5;
	ng_mesh_param_vol_.secondorder=0;

//	ng_mesh_param_vol_.uselocalh=1;


	// netgen init
	Ng_Init();

	// creation of the netgen mesh and the stl geometry to mesh
	ng_mesh_ = Ng_NewMesh();
	//stl_geometry_ = Ng_STL_NewGeometry();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::initializeNetgen(const double AQuality,
		  const double ASize){

	// TODO legoff : look into meaning of quality

	if(AQuality<1.0 ||AQuality>2.0)
		throw gmds::GMDSException("The quality parameter must be included in [1.0,2.0]");


	// TODO legoff : look into relationship between size max for surface triangles and size
	// max for volume tetrahedron.
	// Also look into max tet size and desired tet size in the box.

	// surface parameters
	ng_mesh_param_surf_.maxh =ASize;
//	ng_mesh_param_surf_.maxh =1.;
	ng_mesh_param_surf_.fineness = 0.5;
	ng_mesh_param_surf_.secondorder=0;

//	ng_mesh_param_surf_.uselocalh=1;
//	ng_mesh_param_surf_.grading=1.;

	// volume parameters
	ng_mesh_param_vol_.maxh =ASize*100.;
//	ng_mesh_param_vol_.maxh =2.;
	ng_mesh_param_vol_.fineness = 0.5;
	ng_mesh_param_vol_.secondorder=0;

//	ng_mesh_param_vol_.uselocalh=1;
//	ng_mesh_param_vol_.grading=1.;

	// netgen init
	nglib::Ng_Init();

	// creation of the netgen mesh and the stl geometry to mesh
	ng_mesh_ = nglib::Ng_NewMesh();
	//stl_geometry_ = nglib::Ng_STL_NewGeometry();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::finalizeNetgen(){
	// we delete the netgen mesh
	Ng_DeleteMesh(ng_mesh_);
	// nothing is given in Netgen to delete STL objects, but we delete it too
//	delete stl_geometry_;

	Ng_Exit();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMesh(
		gmds::geom::GeomVolume& AVol,
		gmds::geom::GeomServiceAbstractFactory& AMan,
		gmds::IGMesh& AVolMesh,
		const double AQuality)
{
	initializeNetgen(AQuality);
	// first of all we build a surfacic mesh from ABoundaryMesh
	// we create a STL geometry...
	buildSTLSurfaceFrom(AVol,AMan);

	createGMDSTriMesh(AVolMesh);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMesh(
		gmds::geom::GeomVolume& AVol,
		gmds::geom::GeomServiceAbstractFactory& AMan,
		gmds::IGMesh& AVolMesh,
		const double AQuality,
		const double ASize)
{
	initializeNetgen(AQuality,ASize);
	// first of all we build a surfacic mesh from ABoundaryMesh
	// we create a STL geometry...
	buildSTLSurfaceFrom(AVol,AMan);

	createGMDSTriMesh(AVolMesh);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMesh(gmds::geom::GeomSurface& ASurf,
		gmds::geom::GeomServiceAbstractFactory& AMan,
		gmds::IGMesh& AVolMesh,
		const double AQuality)
{
	initializeNetgen(AQuality);
	// first of all we build a surfacic mesh from ABoundaryMesh
	// we create a STL geometry...
	buildSTLSurfaceFrom(ASurf,AMan);
	createGMDSTriMesh(AVolMesh);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTetMesh(gmds::geom::GeomVolume& AVol,
									  gmds::geom::GeomServiceAbstractFactory& AMan,
									  gmds::IGMesh& AVolMesh,
									  const double AQuality)
{
	initializeNetgen(AQuality);
	// first of all we build a surface mesh from ABoundaryMesh
	// we create a STL geometry...
	buildSTLSurfaceFrom(AVol,AMan);

	createGMDSMesh(AVolMesh);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTetMesh(gmds::IGMesh& ABoundaryMesh,
									 gmds::IGMesh& AVolMesh,
									 const double AQuality){

	//<DIM3|F|N|F2N> -> model of boundary mesh
	if (ABoundaryMesh.getNbQuadrilaterals()!=0 )
		//||ABoundaryMesh.getNbPolygons()!=0 )
		throw gmds::GMDSException("The boundary mesh cannot contain quads or polygons");

	initializeNetgen(AQuality);

	buildSTLSurfaceFrom(ABoundaryMesh);

	createGMDSMesh(AVolMesh);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTetMesh(gmds::geom::GeomVolume& AVol,
									  gmds::geom::GeomServiceAbstractFactory& AMan,
									  gmds::IGMesh& AVolMesh,
									  const double AQuality,
									  gmds::math::Point& APMin,
									  gmds::math::Point& APMax,
									  const double ASize)
{
	initializeNetgen(AQuality,APMin,APMax,ASize);
	// first of all we build a surfacic mesh from ABoundaryMesh
	// we create a STL geometry...
	buildSTLSurfaceFrom(AVol,AMan);

	createGMDSMesh(AVolMesh,APMin,APMax,ASize);

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/
void NetgenFacade::generateTetMesh(
		std::vector<gmds::Node>& ANodes,
		std::vector<gmds::Face>& AFaces,
		gmds::IGMesh& AMesh)
{

	initializeNetgen(1.8);

	int nbInitNodes = ANodes.size();

	std::vector<gmds::TCellID> netgen2gmds_nodes;
	std::map<gmds::TCellID,int> gmds2netgen_nodes;
	//first vector item is not used and thus initiliased
	netgen2gmds_nodes.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i< nbInitNodes; i++){
		gmds::Node n = ANodes[i];
		double c[3];
		c[0] = n.X();
		c[1] = n.Y();
		c[2] = n.Z();
		Ng_AddPoint(ng_mesh_, c);
		netgen2gmds_nodes.push_back(n.getID());
		gmds2netgen_nodes.insert(std::pair<gmds::TCellID,int>(n.getID(),i+1/*+1*/));
	}


	std::vector<gmds::TCellID> netgen2gmds_faces;
	//first vector item is not used and thus initiliased
	netgen2gmds_faces.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i<AFaces.size();i++){
		gmds::Face f = AFaces[i];
		std::vector<gmds::TCellID> f_node_ids;
		f.getIDs<Node>(f_node_ids);
		int ids[3];
		ids[0] = gmds2netgen_nodes[f_node_ids[0]];
		ids[1] = gmds2netgen_nodes[f_node_ids[1]];
		ids[2] = gmds2netgen_nodes[f_node_ids[2]];
		Ng_AddSurfaceElement(ng_mesh_, NG_TRIG, ids);
		netgen2gmds_faces.push_back(f.getID());
	}


	Ng_SurfaceMeshOrientation (ng_mesh_);
	// VOLUME MESHING
	NgAddOn_GenerateVolumeMesh(ng_mesh_,2);
//	NgAddOn_OptimizeVolumeMesh(ng_mesh_,2);

//	// GMDS CELL CREATION
	int nbNodes = Ng_GetNP(ng_mesh_);
	int nbTri   = Ng_GetNSE(ng_mesh_);
	int nbTet   = Ng_GetNE(ng_mesh_);

//	std::cout<<"Result"<<"\n";
//	std::cout<<"\t number of mesh points    : "<<nbNodes<<"\n";
//	std::cout<<"\t number of boundary faces : "<<nbTri<<"\n";
//	std::cout<<"\t number of elements (tets): "<<nbTet<<"\n";

	std::vector<Node> netgen2GMDSNode;
	netgen2GMDSNode.reserve(nbNodes);
	Node n;
	for(int i = 1; i < nbNodes+1; i++){
		/* First, we try to find if a node is already located in this
		 * location. It is so a very unstable (unpredicatable) function
		 * depending on the mesh precision
		 */
		double coord[3];
		Ng_GetPoint(ng_mesh_,i,coord);
//		std::cout<<i<<" - ("<<coord[0]<<", "<<coord[1]<<", "<<coord[2]<<")"<<"\n";
		gmds::math::Point p_netgen(coord[0],coord[1],coord[2]);
		bool found=false;
		gmds::TCellID found_id;
		unsigned int k=0;
		for(;k<ANodes.size() && !found;k++){
			gmds::math::Point p_init = ANodes[k].getPoint();
			if(p_init==p_netgen){
				found =true;
				found_id=ANodes[k].getID();
			}
		}

		if(!found){
			//new node. we add it
			n= AMesh.newNode(coord[0],coord[1],coord[2]);
			netgen2GMDSNode.push_back(n);
		}
		else{
			//existing node
			n = AMesh.get<Node>(found_id);
			netgen2GMDSNode.push_back(n);
		}
	}

//	std::vector<Node*> Netgen2GMDSNode;
//	Netgen2GMDSNode.resize(nbNodes+1);
//
//	for(int i = 1; i< nbInitNodes+1; i++){
//		std::cout<<i<<" -> "<<netgen2gmds_nodes[i];
//		Netgen2GMDSNode[i]= AMesh.get<Node>(netgen2gmds_nodes[i]);
//		Node* n = AMesh.get<Node>(netgen2gmds_nodes[i]);;
//		std::cout<<"\t"<<n.getX()<<" "<<n.getY()<<" "<<n.getZ()<<"\n";
//	}
//
//	for(int i = nbInitNodes+1; i < nbNodes+1; i++){
//		double coord[3];
//		Ng_GetPoint(ng_mesh_,i,coord);
//		Netgen2GMDSNode[i]= AMesh.newNode(coord[0],coord[1],coord[2]);
//
//	}

	for(int i = 1; i <nbTet+1; i++){
		int index[4];
		Ng_GetVolumeElement(ng_mesh_,i,index);

		AMesh.newTet(netgen2GMDSNode[index[0]-1],netgen2GMDSNode[index[1]-1],
				netgen2GMDSNode[index[2]-1],netgen2GMDSNode[index[3]-1]);

//		std::cout<<"tet N"<<i<<" - "<<index[0]<<", "<<index[1]<<", "<<index[2]<<", "<<index[3]<<"\n";
//		std::cout<<"tet G "<<i<<" - "<<netgen2GMDSNode[index[0]-1].getID()<<", "<<
//									netgen2GMDSNode[index[1]-1].getID()<<", "<<
//									netgen2GMDSNode[index[2]-1].getID()<<", "<<
//									netgen2GMDSNode[index[3]-1].getID()<<"\n";
}


	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/


void NetgenFacade::generateTetMesh(
		std::vector<gmds::Node>& ANodes,
		std::vector<gmds::Face>& AFaces,
		std::vector<gmds::Node>& ACreatedNodes,
		std::vector<gmds::Region>& ACreatedRegions,
		gmds::IGMesh& AMesh,
		NetgenParams& AParams)
{

	initializeNetgen(AParams);

	int nbInitNodes = ANodes.size();

	std::vector<gmds::TCellID> netgen2gmds_nodes;
	std::map<gmds::TCellID,int> gmds2netgen_nodes;
	//first vector item is not used and thus initiliased
	netgen2gmds_nodes.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i< nbInitNodes; i++){
		gmds::Node n = ANodes[i];
		double c[3];
		c[0] = n.X();
		c[1] = n.Y();
		c[2] = n.Z();
		Ng_AddPoint(ng_mesh_, c);
		netgen2gmds_nodes.push_back(n.getID());
		gmds2netgen_nodes.insert(std::pair<gmds::TCellID,int>(n.getID(),i+1/*+1*/));
	}


	std::vector<gmds::TCellID> netgen2gmds_faces;
	//first vector item is not used and thus initiliased
	netgen2gmds_faces.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i<AFaces.size();i++){
		gmds::Face f = AFaces[i];
		std::vector<gmds::TCellID> f_node_ids;
		f.getIDs<Node>(f_node_ids);
		int ids[3];
		ids[0] = gmds2netgen_nodes[f_node_ids[0]];
		ids[1] = gmds2netgen_nodes[f_node_ids[1]];
		ids[2] = gmds2netgen_nodes[f_node_ids[2]];
		Ng_AddSurfaceElement(ng_mesh_, NG_TRIG, ids);
		netgen2gmds_faces.push_back(f.getID());
	}


	Ng_SurfaceMeshOrientation (ng_mesh_);
	// VOLUME MESHING
	NgAddOn_GenerateVolumeMesh(ng_mesh_,2);
	NgAddOn_OptimizeVolumeMesh(ng_mesh_,2);


    // déclaration d'une géométrie OCC-netgen
//	netgen::OCCGeometry occgeo;


//    occgeo.occdeflection = 0.01;
//
//
//    TopoDS_Shape shape_loc = AFace;
//    occgeo.shape = shape_loc;
//    occgeo.changed = 1;
//
//    occgeo.BuildFMap();
//
//    // nécessaire car les maillages par défaut OCC y sont générés
//    occgeo.BuildVisualizationMesh();
//
//    // Déclaration d'un maillage Netgen
//    Mesh *maill = new Mesh();
//
//    maill->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));

    ng_mesh_param_surf_.maxh =2;
    ng_mesh_param_surf_.fineness = 0.5;
    ng_mesh_param_surf_.secondorder=0;
    //ng_mesh_param_surf_.

//    netgen::mparam.quad = 0;
//    netgen::mparam.grading = m_grading;
//    netgen::mparam.curvaturesafety = m_curve;
//    netgen::mparam.segmentsperedge = m_edge;
//    netgen::mparam.uselocalh = 1;
//    netgen::mparam.checkoverlap = 0;

//    OCCGenerateMesh( occgeo, maill,
//                             MESHCONST_ANALYSE,
//                             MESHCONST_OPTSURFACE, "maillage");

//	// GMDS CELL CREATION
	int nbNodes = Ng_GetNP(ng_mesh_);
	int nbTri   = Ng_GetNSE(ng_mesh_);
	int nbTet   = Ng_GetNE(ng_mesh_);

	std::vector<Node> netgen2GMDSNode;
	netgen2GMDSNode.reserve(nbNodes);
	Node n;
	for(int i = 1; i < nbNodes+1; i++){
		/* First, we try to find if a node is already located in this
		 * location. It is so a very unstable (unpredicatable) function
		 * depending on the mesh precision
		 */
		double coord[3];
		Ng_GetPoint(ng_mesh_,i,coord);
		gmds::math::Point p_netgen(coord[0],coord[1],coord[2]);
		bool found=false;
		gmds::TCellID found_id;
		unsigned int k=0;
		for(;k<ANodes.size() && !found;k++){
			gmds::math::Point p_init = ANodes[k].getPoint();
			if(p_init==p_netgen){
				found =true;
				found_id=ANodes[k].getID();
			}
		}

		if(!found){
			//new node. we add it
			n= AMesh.newNode(coord[0],coord[1],coord[2]);
			netgen2GMDSNode.push_back(n);
			ACreatedNodes.push_back(n);
		}
		else{
			//existing node
			n = AMesh.get<Node>(found_id);
			netgen2GMDSNode.push_back(n);
		}
	}

	for(int i = 1; i <nbTet+1; i++){
		int index[4];
		Ng_GetVolumeElement(ng_mesh_,i,index);

		Region r =AMesh.newTet(netgen2GMDSNode[index[0]-1],netgen2GMDSNode[index[1]-1],
				netgen2GMDSNode[index[2]-1],netgen2GMDSNode[index[3]-1]);

		ACreatedRegions.push_back(r);
}


	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTetMesh(gmds::IGMesh& AMesh, const double ASize)
{

	initializeNetgen(1.2,ASize);

	// copy nodes and faces
	// TODO legoff : useless, needs to be replaced later
	std::vector<gmds::Node> ANodes;
	std::vector<gmds::Face> AFaces;
	{
		gmds::IGMesh::node_iterator it  = AMesh.nodes_begin();

		for(;!it.isDone();it.next())
		{
			ANodes.push_back(it.value());
		}
	}
	{
		gmds::IGMesh::face_iterator it  = AMesh.faces_begin();

		for(;!it.isDone();it.next())
		{
			AFaces.push_back(it.value());
		}
	}

	// FILL THE NETGEN DS AND KEEP A LINK FROM GMDS TO NETGEN ENTITIES
	int nbInitNodes = ANodes.size();
	int nbInitFaces = AFaces.size();
	std::vector<gmds::TCellID> netgen2gmds_nodes;
	std::map<gmds::TCellID,int> gmds2netgen_nodes;
	//first vector item is not used and thus initiliased
	netgen2gmds_nodes.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i<ANodes.size();i++){
		gmds::Node n = ANodes[i];
		double c[3];
		c[0] = n.X();
		c[1] = n.Y();
		c[2] = n.Z();
		nglib::Ng_AddPoint(ng_mesh_, c);
		netgen2gmds_nodes.push_back(n.getID());
//		std::cout<<"Node "<<i<<" : "<<c[0]<<" "<<c[1]<<" "<<c[2]<<"\n";
		gmds2netgen_nodes.insert(std::pair<gmds::TCellID,int>(n.getID(),i+1));
	}

	std::vector<gmds::TCellID> netgen2gmds_faces;
	//first vector item is not used and thus initiliased
	netgen2gmds_faces.push_back(gmds::NullID);
	// NETGEN ENTITES ARE ORDERED FROM 1 to N
	for(unsigned int i=0; i<AFaces.size();i++){
		gmds::Face f = AFaces[i];
		std::vector<gmds::TCellID> f_node_ids;
		f.getIDs<Node>(f_node_ids);
		int ids[3];
		ids[0] = gmds2netgen_nodes[f_node_ids[0]];
		ids[1] = gmds2netgen_nodes[f_node_ids[1]];
		ids[2] = gmds2netgen_nodes[f_node_ids[2]];
//		std::cout<<"Triangle "<<ids[0]<<" "<<ids[1]<<" "<<ids[2]<<"\n";
		nglib::Ng_AddSurfaceElement(ng_mesh_, nglib::NG_TRIG, ids);
		netgen2gmds_faces.push_back(f.getID());
	}


	// VOLUME MESHING
	Ng_SurfaceMeshOrientation (ng_mesh_);

	 nglib::Ng_Meshing_Parameters mp;
	  mp.maxh = 1;
	  mp.fineness = 1;
	  mp.secondorder = 0;
	nglib::Ng_GenerateVolumeMesh(ng_mesh_,&mp);

	// GMDS CELL CREATION
	int nbNodes = nglib::Ng_GetNP(ng_mesh_);
	int nbTri   = nglib::Ng_GetNSE(ng_mesh_);
	int nbTet   = nglib::Ng_GetNE(ng_mesh_);

//	std::cout<<"Result"<<"\n";
//	std::cout<<"\t number of points: "<<nbNodes<<"\n";
//	std::cout<<"\t number of faces : "<<nbTri<<"\n";
//	std::cout<<"\t number of tet   : "<<nbTet<<"\n";

	std::vector<Node> Netgen2GMDSNode;
	Netgen2GMDSNode.resize(nbNodes+1);

	for(int i = 1; i< nbInitNodes+1; i++){
//		std::cout<<i<<" -> "<<netgen2gmds_nodes[i]<<"\n";
		Netgen2GMDSNode[i]= AMesh.get<Node>(netgen2gmds_nodes[i]);
		Node n = AMesh.get<Node>(netgen2gmds_nodes[i]);;
//		std::cout<<"\t"<<n.getX()<<" "<<n.getY()<<" "<<n.getZ()<<"\n";
	}

	for(int i = nbInitNodes+1; i < nbNodes+1; i++){
		double coord[3];
		nglib::Ng_GetPoint(ng_mesh_,i,coord);
		Netgen2GMDSNode[i]= AMesh.newNode(coord[0],coord[1],coord[2]);

	}

	for(int i = 1; i <nbTet+1; i++){
		int index[4];
		nglib::Ng_GetVolumeElement(ng_mesh_,i,index);
		AMesh.newTet(Netgen2GMDSNode[index[0]],Netgen2GMDSNode[index[1]],
				Netgen2GMDSNode[index[2]],Netgen2GMDSNode[index[3]]);
	}

	finalizeNetgen();
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::createGMDSMesh(gmds::IGMesh& AMesh)
{
	//then we generate the surface mesh
	//Ng_STL_GenerateSurfaceMesh(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);
	// and the volume mesh
	Ng_GenerateVolumeMesh(ng_mesh_,&ng_mesh_param_vol_);
	//finally we build the corresponding GMDS 3D mesh
	buildGMDSOutput(AMesh);
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::createGMDSMesh(gmds::IGMesh& AMesh,
		  gmds::math::Point& APMin,
		  gmds::math::Point& APMax,
		  const double ASize)
{
	//then we generate the surface mesh
//	Ng_STL_GenerateSurfaceMesh(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);

	// and the volume mesh
	double pmin[3];
	pmin[0] = APMin.X();
	pmin[1] = APMin.Y();
	pmin[2] = APMin.Z();
	double pmax[3];
	pmax[0] = APMax.X();
	pmax[1] = APMax.Y();
	pmax[2] = APMax.Z();

	Ng_RestrictMeshSizeBox(ng_mesh_,pmin,pmax,ASize);
	Ng_GenerateVolumeMesh(ng_mesh_,&ng_mesh_param_vol_);
	//finally we build the corresponding GMDS 3D mesh
	buildGMDSOutput(AMesh);
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::createGMDSTriMesh(gmds::IGMesh& AMesh)
{
	//then we generate the surface mesh
//	Ng_STL_GenerateSurfaceMesh(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);
	// and the volume mesh
	//Ng_GenerateVolumeMesh(ng_mesh_,&ng_mesh_param_vol_);
	//finally we build the corresponding GMDS 3D mesh
	buildGMDSOutput(AMesh);
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::buildSTLSurfaceFrom(IGMesh& AMesh)
{
	IGMesh::surfaces_iterator its   = AMesh.surfaces_begin();
	IGMesh::surfaces_iterator its_e = AMesh.surfaces_end();


	while(its != its_e){
		IGMesh::surface& surf = *its;

		std::vector<Face> faces = surf.cells();
		std::vector<Face>::iterator it=faces.begin(), ite=faces.end();
		int face_index = 0;
		// we add a STL triangle for each face
		while(it!=ite){
			Face f_GMDS = *it;
			std::vector<Node> nodes;
			f_GMDS.get<Node>(nodes);
			double N[3][3];

			for(int node_i =0; node_i<nodes.size();node_i++)
			{
				Node ni = nodes[node_i];
				N[node_i][0]= ni.X();
				N[node_i][1]= ni.Y();
				N[node_i][2]= ni.Z();
			}
			it++;
			/* we add the corresponding triangle into the STL geometry without
			 * indicating its normal direction */
	//		Ng_STL_AddTriangle(stl_geometry_,N[0],N[1],N[2],0);
		}
		its++;
	}


//	Ng_STL_InitSTLGeometry(stl_geometry_);
//	Ng_STL_MakeEdges(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);

}
/*----------------------------------------------------------------------------*/
 void NetgenFacade::
buildSTLSurfaceFrom(gmds::geom::GeomVolume& AVol,
					gmds::geom::GeomServiceAbstractFactory& AMan)
{
	std::vector<gmds::geom::GeomSurface* > surfs;
	AVol.get(surfs);

	// we get the triangles from the geometric volume AVol
	gmds::geom::GeomTriangulationService *s = AMan.newGeomTriangulationService();
	std::vector<gmds::Face> triangles;
	s->perform(AVol,triangles);

	for(unsigned int i=0;i<triangles.size();i++)
	{
		double N[3][3];
		std::vector<Node> nodes = triangles[i].get<Node>();
		gmds::math::Point p = nodes[0].getPoint();
		N[0][0] = p.X();
		N[0][1] = p.Y();
		N[0][2] = p.Z();

		p = nodes[1].getPoint();
		N[1][0] = p.X();
		N[1][1] = p.Y();
		N[1][2] = p.Z();

		p = nodes[2].getPoint();
		N[2][0] = p.X();
		N[2][1] = p.Y();
		N[2][2] = p.Z();

		/* we add the corresponding triangle into the STL geometry without
		 * indicating its normal direction */
	//	Ng_STL_AddTriangle(stl_geometry_,N[0],N[1],N[2],0);
	}

//
//	Ng_STL_InitSTLGeometry(stl_geometry_);
//	Ng_STL_MakeEdges(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);

}
/*----------------------------------------------------------------------------*/
 void NetgenFacade::
buildSTLSurfaceFrom(gmds::geom::GeomSurface& ASurf,
					gmds::geom::GeomServiceAbstractFactory& AMan)
{
	// we get the triangles from the geometric volume AVol
	gmds::geom::GeomTriangulationService *s =
			AMan.newGeomTriangulationService();
	std::vector<gmds::Face > triangles;
	s->perform(ASurf,triangles);

	for(unsigned int i=0;i<triangles.size();i++)
	{
		double N[3][3];
		std::vector<Node> nodes = triangles[i].get<Node>();
		gmds::math::Point p = nodes[0].getPoint();
		N[0][0] = p.X();
		N[0][1] = p.Y();
		N[0][2] = p.Z();

		p = nodes[1].getPoint();
		N[1][0] = p.X();
		N[1][1] = p.Y();
		N[1][2] = p.Z();

		p = nodes[2].getPoint();
		N[2][0] = p.X();
		N[2][1] = p.Y();
		N[2][2] = p.Z();

		/* we add the corresponding triangle into the STL geometry without
		 * indicating its normal direction */
	//	Ng_STL_AddTriangle(stl_geometry_,N[0],N[1],N[2],0);
	}
//	std::cout<<"NB TRIANGLES -----> "<<triangles.size()<<"\n";

//	Ng_STL_InitSTLGeometry(stl_geometry_);

	// no edge to build we do not have a volume !!
//	Ng_STL_MakeEdges(stl_geometry_,ng_mesh_,&ng_mesh_param_surf_);

}
/*----------------------------------------------------------------------------*/

void NetgenFacade::buildGMDSOutput(gmds::IGMesh& AMesh)
{

	int nbNodes = Ng_GetNP(ng_mesh_);
	int nbTri   = Ng_GetNSE(ng_mesh_);
	int nbTet   = Ng_GetNE(ng_mesh_);

	std::cout<<"Result"<<"\n";
	std::cout<<"\t number of points: "<<nbNodes<<"\n";
	std::cout<<"\t number of faces : "<<nbTri<<"\n";
	std::cout<<"\t number of tet   : "<<nbTet<<"\n";

	std::vector<Node> Netgen2GMDSNode;
	Netgen2GMDSNode.resize(nbNodes+1);

	for(int i = 1; i < nbNodes+1; i++){
		double coord[3];
		Ng_GetPoint(ng_mesh_,i,coord);
		Netgen2GMDSNode[i]= AMesh.newNode(coord[0],coord[1],coord[2]);

	}

	for(int i = 1; i <nbTet+1; i++){
		int index[4];
		Ng_GetVolumeElement(ng_mesh_,i,index);
		AMesh.newTet(Netgen2GMDSNode[index[0]],Netgen2GMDSNode[index[1]],
				Netgen2GMDSNode[index[2]],Netgen2GMDSNode[index[3]]);
	}


	{
		for(int i = 1; i < nbTri+1; i++){
			int index[3];
					Ng_GetSurfaceElement(ng_mesh_,i,index);
			Face f = AMesh.newTriangle(Netgen2GMDSNode[index[0]],
										Netgen2GMDSNode[index[1]],
										Netgen2GMDSNode[index[2]]);
		}
	}
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::buildGMDSOutput(gmds::IGMesh& AMesh,
			std::set<TCellID>& initNodes,
		 std::map<TCellID,TCellID>& ANodeMap)
{

	int nbNodes = Ng_GetNP(ng_mesh_);
	int nbTri   = Ng_GetNSE(ng_mesh_);
	int nbTet   = Ng_GetNE(ng_mesh_);

	std::cout<<"Result"<<"\n";
	std::cout<<"\t number of points: "<<nbNodes<<"\n";
	std::cout<<"\t number of faces : "<<nbTri<<"\n";
	std::cout<<"\t number of tet   : "<<nbTet<<"\n";

	std::vector<Node> Netgen2GMDSNode;
	Netgen2GMDSNode.resize(nbNodes+1);

	for(int i = 1; i < nbNodes+1; i++){
		double coord[3];
		Ng_GetPoint(ng_mesh_,i,coord);

		Netgen2GMDSNode[i]= AMesh.newNode(coord[0],coord[1],coord[2]);

		//map building
		std::set<TCellID>::iterator it_n;
		bool not_found = true;
		for(it_n = initNodes.begin(); it_n!=initNodes.end() && not_found; it_n++){
			Node n = AMesh.get<Node>(*it_n);
			double sum = 0;
			sum += ((n.X()-Netgen2GMDSNode[i].X()) * (n.X()-Netgen2GMDSNode[i].X()));

			sum += ((n.Y()-Netgen2GMDSNode[i].Y()) * (n.Y()-Netgen2GMDSNode[i].Y()));

			sum += ((n.Z()-Netgen2GMDSNode[i].Z()) * (n.Z()-Netgen2GMDSNode[i].Z()));

			if(sum < 1e-6){
				not_found=false;
				ANodeMap.insert(std::pair<TCellID,TCellID>(Netgen2GMDSNode[i].getID(),n.getID()));
			}
		}

	}

	for(int i = 1; i <nbTet+1; i++){
		int index[4];
		Ng_GetVolumeElement(ng_mesh_,i,index);
		AMesh.newTet(Netgen2GMDSNode[index[0]],Netgen2GMDSNode[index[1]],
				Netgen2GMDSNode[index[2]],Netgen2GMDSNode[index[3]]);
	}


	{
		for(int i = 1; i < nbTri+1; i++){
			int index[3];
					Ng_GetSurfaceElement(ng_mesh_,i,index);
			Face f = AMesh.newTriangle(Netgen2GMDSNode[index[0]],
										Netgen2GMDSNode[index[1]],
										Netgen2GMDSNode[index[2]]);
		}
	}
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMeshFromSTEP(const std::string& AName,
	 	 	 gmds::IGMesh& ASurfMesh,
		 const double AQuality)

{
	initializeNetgen(AQuality);;

	Ng_OCC_Geometry * occgeo = Ng_OCC_Load_STEP(AName.c_str());

	ng_mesh_param_surf_.grading = 0.1;
	ng_mesh_param_surf_.fineness = 1;
	ng_mesh_param_surf_.closeedgeenable = 1;
	ng_mesh_param_surf_.elementsperedge=10;
	ng_mesh_param_surf_.elementspercurve=10;
	ng_mesh_param_surf_.maxh=2;
	// nécessaire car les maillages par défaut OCC y sont générés
	Ng_OCC_SetLocalMeshSize   (occgeo, ng_mesh_,&ng_mesh_param_surf_);

	Ng_OCC_GenerateEdgeMesh(occgeo, ng_mesh_,&ng_mesh_param_surf_);

	Ng_OCC_GenerateSurfaceMesh(occgeo, ng_mesh_,&ng_mesh_param_surf_);

	buildGMDSOutput(ASurfMesh);

	finalizeNetgen();

}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMesh(const std::string& AName,
		 gmds::IGMesh& ASurfMesh,
		 std::vector<Node>& ANodes,
		 std::vector<std::vector<Node> >& AEdges,
		 const double AQuality)

{
	initializeNetgen(AQuality);;

	Ng_OCC_Geometry * occgeo = Ng_OCC_Load_STEP(AName.c_str());

	ng_mesh_param_surf_.grading = 0.1;
	ng_mesh_param_surf_.fineness = 1;
	ng_mesh_param_surf_.closeedgeenable = 1;
	ng_mesh_param_surf_.elementsperedge=10;
	ng_mesh_param_surf_.elementspercurve=10;
	ng_mesh_param_surf_.maxh=2;
	// nécessaire car les maillages par défaut OCC y sont générés
	Ng_OCC_SetLocalMeshSize   (occgeo, ng_mesh_,&ng_mesh_param_surf_);


	std::vector<double> nodes_xyz;

	for(unsigned int i=0;i<ANodes.size();i++){
		Node n = ANodes[i];
		nodes_xyz.push_back(n.X());
		nodes_xyz.push_back(n.Y());
		nodes_xyz.push_back(n.Z());
	}

	std::vector<std::vector<double> > edges_xyz;
	for(unsigned int i=0;i<AEdges.size();i++){
		std::vector<Node> edge_nodes = AEdges[i];
		std::vector<double> xyz;
		for (unsigned int j=0;j<edge_nodes.size();j++){
			Node n = edge_nodes[j];
			xyz.push_back(n.X());
			xyz.push_back(n.Y());
			xyz.push_back(n.Z());
		}
		edges_xyz.push_back(xyz);
	}

    Ng_OCC_InitEdges(occgeo, ng_mesh_,
    				 nodes_xyz, edges_xyz,&ng_mesh_param_surf_);


	Ng_OCC_GenerateSurfaceMesh(occgeo, ng_mesh_,&ng_mesh_param_surf_);

	buildGMDSOutput(ASurfMesh);

	finalizeNetgen();

}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMesh(TopoDS_Face& AFace,
		 gmds::IGMesh& ASurfMesh,
		 std::vector<Node>& ANodes,
		 std::vector<std::vector<Node> >& AEdges,
		 std::map<TCellID,TCellID>& ANodeMap,
		  NetgenParams& AParams)

{
	std::set<TCellID> init_nodes;
    //----------------------------------------------
    //---------------------- OLD -------------------
    //----------------------------------------------
    // déclaration d'une géométrie OCC-netgen
	initializeNetgen(AParams);

	Ng_OCC_Geometry * occgeo = Ng_OCC_Load_Shape(AFace);

	// nécessaire car les maillages par défaut OCC y sont générés
	Ng_OCC_SetLocalMeshSize   (occgeo, ng_mesh_,&ng_mesh_param_surf_);


    // Déclaration d'un maillage Netgen




//	// nécessaire car les maillages par défaut OCC y sont générés
//	Ng_OCC_SetLocalMeshSize   (occgeo, ng_mesh_,&ng_mesh_param_surf_);

	std::vector<double> nodes_xyz;
	for(unsigned int i=0;i<ANodes.size();i++){
		Node n = ANodes[i];
		nodes_xyz.push_back(n.X());
		nodes_xyz.push_back(n.Y());
		nodes_xyz.push_back(n.Z());
		init_nodes.insert(n.getID());
	}

	std::vector<std::vector<double> > edges_xyz;
	for(unsigned int i=0;i<AEdges.size();i++){
		std::vector<Node> edge_nodes = AEdges[i];
		std::vector<double> xyz;
		for (unsigned int j=0;j<edge_nodes.size();j++){
			Node n = edge_nodes[j];
			xyz.push_back(n.X());
			xyz.push_back(n.Y());
			xyz.push_back(n.Z());
			init_nodes.insert(n.getID());
		}
		edges_xyz.push_back(xyz);
	}


    Ng_OCC_InitEdges(occgeo, ng_mesh_,
    				 nodes_xyz, edges_xyz,&ng_mesh_param_surf_);
//  Ng_OCC_GenerateEdgeMesh (occgeo,ng_mesh_,&ng_mesh_param_surf_);
//
//
	Ng_OCC_GenerateSurfaceMesh(occgeo, ng_mesh_,&ng_mesh_param_surf_);

//	Ng_OCC_MyGenerateSurfaceMesh(AFace, ng_mesh_,&ng_mesh_param_surf_);
    //----------------------------------------------
    //----------------------------------------------
	buildGMDSOutput(ASurfMesh, init_nodes, ANodeMap);

	finalizeNetgen();

}
/*----------------------------------------------------------------------------*/

void NetgenFacade::generateTriMeshWithOCC(
		std::vector<gmds::math::Point >& APnts,
			 	 	 	 	 gmds::IGMesh& ASurfMesh,
							 const double AQuality ){

	initializeNetgen(AQuality);
	// FACE CONSTRUITE
    /** les coordonnées en X Y, Z des points */
	int nbPoints = APnts.size();


    for(unsigned int i=0;i<nbPoints;i++)
    {
    	double coords[3];
    	coords[0] = APnts[i].X();
    	coords[1] = APnts[i].Y();
    	coords[2] = APnts[i].Z();
    	Ng_AddPoint(ng_mesh_,coords);
    }
    int pi[3];
    pi[0] = 1;
    pi[1] = 2;
    pi[2] = 3;
    Ng_AddSurfaceElement (ng_mesh_, NG_TRIG, pi);


    pi[0] = 1;
    pi[1] = 3;
    pi[2] = 4;
    Ng_AddSurfaceElement (ng_mesh_, NG_TRIG, pi);

	std::cout<<"OOC FACE CREATION"<<"\n";
	TopoDS_Face f = buildOCCSurfaceFromPoints(APnts);
	std::cout<<"NETGEN MESHING"<<"\n";
	meshing2D(f);
	std::cout<<"GMDS WRITING"<<"\n";
	//MAILLAGE AVEC NETGEN
	buildGMDSOutput(ASurfMesh);

	finalizeNetgen();
}

/*----------------------------------------------------------------------------*/

TopoDS_Face NetgenFacade::buildOCCSurfaceFromPoints(
			std::vector<gmds::math::Point >& APnts)
{
    /** les coordonnées en X Y, Z des points */
	int nbPoints = APnts.size();
    double *x, *y, *z;
    x = new double[nbPoints];
    y = new double[nbPoints];
    z = new double[nbPoints];
    for(unsigned int i=0;i<nbPoints;i++)
    {
    	x[i] = APnts[i].X();
    	y[i] = APnts[i].Y();
    	z[i] = APnts[i].Z();
    }

    // PART 1 == CREATION D'UNE LIGNE OUVERTE BORDANT LA FACE
    gp_Pnt  aPnt1;
    gp_Pnt  aPnt2;
    bool wireVide = true;
    Handle(Geom_Curve) result_curve;
    std::vector<gp_Pnt> cleaned_points;

    for(int i=0; i<nbPoints - 1;++i)
    {
        aPnt1 = gp_Pnt(x[i  ],y[i  ],z[i  ]);
        aPnt2 = gp_Pnt(x[i+1],y[i+1],z[i+1]);
        if (aPnt1.IsEqual(aPnt2,1e-12))
        {
     //       std::cout<<"Make Polyline: null length segment detected in MakePolyline!"<<"\n";
            continue;
        }
        else
            cleaned_points.push_back(aPnt1);
    }

    cleaned_points.push_back(gp_Pnt(x[nbPoints-1],
                                    y[nbPoints-1],
                                    z[nbPoints-1]) );

    for(int i=0; i<cleaned_points.size()-1;++i)
    {
        aPnt1 = cleaned_points[i];
        aPnt2 = cleaned_points[i+1];
        wireVide=false;

        GC_MakeSegment seg=GC_MakeSegment(aPnt1,  aPnt2);

        if (i==0)
            result_curve = seg.Value();
        else
        {

            Handle(Geom_Curve) tmp;
            Standard_Real sta1 = result_curve->FirstParameter();
            Standard_Real sta2 = (seg.Value())->FirstParameter();
            Standard_Real end1 = result_curve->LastParameter();
            Standard_Real end2 = (seg.Value())->LastParameter();
            Handle(Geom_Curve) tmp2 = (Handle(Geom_Curve))(seg.Value());
            Standard_Boolean b1,b2;
            bool succeed = ShapeConstruct::JoinCurves(result_curve,
                                       tmp2,
                                       TopAbs_FORWARD,TopAbs_FORWARD,
                                       sta1,end1,sta2,end2,
                                       tmp,
                                       b1,b2);

            if (!succeed)
            {
 //               std::cout<<"OCeane Make Polyline: impossible to build polyline "<<"\n";
                exit(0);
            }
            result_curve = tmp;
        }
    }

    if (wireVide)
    {
//        cout<<"\n"<<"Erreur : Tous les points de la ligne brisee sont confondus"<<"\n";
//        std::cout << "### Ligne brisee " << "\n";
//        for(int i=0; i<nbPoints; ++i)
//            std::cout << "Point " << i << " "<<x[i] << " , " <<y[i] << " , " <<z[i] << "\n";
//        std::cout << "\n";
        exit(0);
    }

    BRepBuilderAPI_MakeEdge mkEdge(result_curve);
    TopoDS_Edge boundary_edge = mkEdge.Edge();

    // PART 2 == CREATION DU CONTOUR (WIRE)

    Handle(ShapeExtend_WireData) data = new ShapeExtend_WireData();
    data->Add(boundary_edge);

    ShapeFix_Wire builder;
    builder.Load(data);
    builder.Perform();
    TopoDS_Wire w = builder.Wire();//WireAPIMake();

    // PART 3 == CREATION DE LA SURFACE (OU FACE OCC)

    TopoDS_Face aFace;
    BRepBuilderAPI_MakeFace MF(w);
    aFace   = MF.Face();

    return aFace;
}
/*----------------------------------------------------------------------------*/

void NetgenFacade::meshing2D(TopoDS_Face & AFace)
{
    std::cout<<"Maillage triangulaire NETGEN"<<"\n";

//    MoniTool_TimerSentry time("meshing");


    Ng_OCC_Geometry * occgeo = Ng_OCC_Load_Shape(AFace);

   	Ng_OCC_GenerateSurfaceMesh(occgeo, ng_mesh_,&ng_mesh_param_surf_);
   	Ng_OCC_Uniform_Refinement (occgeo, ng_mesh_);
    // déclaration d'une géométrie OCC-netgen
//    Ng_OCC_Geometry occgeo;
//
//    occgeo.occdeflection = 0.01;
//
//
//    TopoDS_Shape shape_loc = AFace;
//    occgeo.shape = shape_loc;
//    occgeo.changed = 1;
//
//    occgeo.BuildFMap();
//
//    // nécessaire car les maillages par défaut OCC y sont générés
//    occgeo.BuildVisualizationMesh();
//
//    // Déclaration d'un maillage Netgen
//    Mesh *maill = new Mesh();
//
//    maill->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
//
//    ng_mesh_param_surf_.maxh =2;
//    ng_mesh_param_surf_.fineness = 0.5;
//    ng_mesh_param_surf_.secondorder=0;
//    //ng_mesh_param_surf_.
//
////    netgen::mparam.quad = 0;
////    netgen::mparam.grading = m_grading;
////    netgen::mparam.curvaturesafety = m_curve;
////    netgen::mparam.segmentsperedge = m_edge;
////    netgen::mparam.uselocalh = 1;
////    netgen::mparam.checkoverlap = 0;
//
//    OCCGenerateMesh( occgeo, maill,
//                             MESHCONST_ANALYSE,
//                             MESHCONST_OPTSURFACE, "maillage");
//
//
//    /* permet de remettre à 0 le timer chaque fois que l'on rentre dans la
//     * méthode. Sinon le temps affiché cumule tous les appels antérieurs à la
//     * méthode. */
////    time.Timer()->ClearTimers();
////    time.Destroy();
}
