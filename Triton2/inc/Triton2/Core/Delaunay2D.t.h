/*----------------------------------------------------------------------------*/
/*
 * Delaunay2D.t.h
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
Delaunay2D<TNum, TGeomKernel>::OrientedBoundary::
OrientedBoundary(const std::vector<gmds::math::Point >& APoints)
:points_(APoints)
{
	/* TODO
	 * Check that we don't have twice the same point
	 */
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
const std::vector<gmds::math::Point >&
Delaunay2D<TNum, TGeomKernel>::OrientedBoundary::getOrderedPoints() const

{
	return points_;
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::initVariables()
{
//	init_mark_ = mesh_.getNewMark();
//	init_N2F_node_variable = mesh_.template newVariable<std::list<gmds::Face> >(gmds::GMDS_NODE,"N2F");
//	/* not use for a simple triangulation but for mesh algorithms where inner
//	 * points must be added */
//	size_field = mesh_.template newVariable<TNum>(gmds::GMDS_NODE,"size");
//
//	/* for each boundary node, we assign a size depending on the incident edges (i.e. distance
//	 * to the prev. and next nodes stored in init_nodes
//	 */
//	for(unsigned int i=0;i<init_nodes_.size();i++){
//		gmds::Node n = init_nodes_[i];
//		gmds::Node prev, next;
//		if(i==0)
//			prev = init_nodes_[init_nodes_.size()-1];
//		else
//			prev = init_nodes_[i-1];
//		if(i==init_nodes_.size()-1)
//			next = init_nodes_[0];
//		else
//			next = init_nodes_[i+1];
//
//		gmds::math::Point p_cur  = n.getPoint();
//		gmds::math::Point p_prev = prev.getPoint();
//		gmds::math::Point p_next = next.getPoint();
//
//		TNum d1 = p_cur.distance(p_prev);
//		TNum d2 = p_cur.distance(p_next);
//
//		(*size_field)[n.getID()] = .5*(d1+d2);
//	}
//
//	adimensional_size_field_ = mesh_.template newVariable<TNum>(gmds::GMDS_FACE,"adimensional size");
//
//	circumCenter_field_ = mesh_.template newVariable<gmds::math::Point >(gmds::GMDS_FACE,"circumCenter");
//
//	delete_mark_ = mesh_.getNewMark();
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
Delaunay2D<TNum,TGeomKernel>::
Delaunay2D(const std::vector<gmds::math::Point >& APoints)
: mesh_(gmds::IGMesh(gmds::DIM2|gmds::N|gmds::F|gmds::F2N|gmds::F2F))

{
//	/* TODO
//	 * Check that we don't have twice the same point
//	 */
//
//	/* we directly create the mesh nodes */
//	typename std::vector<gmds::math::Point >::const_iterator it = APoints.begin();
//	for(;it!=APoints.end();it++)
//	{
//		gmds::math::Point p = *it;
//		init_nodes_.push_back(mesh_.newNode(p.X(),p.Y()));
//		init_points_.push_back(p);
//	}
//
//	for(int i=0;i<init_nodes_.size();i++)
//		mesh_.mark(init_nodes_[i],init_mark_);
//
//	initVariables();
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
Delaunay2D<TNum,TGeomKernel>::
Delaunay2D(const OrientedBoundary& ABoundary)
: mesh_(gmds::IGMesh(gmds::DIM2|gmds::N|gmds::F|gmds::F2N|gmds::F2F))
{

	//	initVariables();
//	addBoundary(ABoundary);
//
//	for(int i=0;i<init_nodes_.size();i++)
//		mesh_.mark(init_nodes_[i],init_mark_);
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
Delaunay2D<TNum, TGeomKernel>::
Delaunay2D(const OrientedBoundary& AOutBoundary,
		const OrientedBoundary* AInBoundary, const int ANbIn)
		: mesh_(gmds::IGMesh(gmds::DIM2|gmds::N|gmds::F|gmds::F2N|gmds::F2F))
		  {
	//	initVariables();
	//	addBoundary(AOutBoundary);
	//	for(unsigned int i=0;i<ANbIn;i++)
	//		addBoundary(AInBoundary[i]);
//
//	for(int i=0;i<init_nodes_.size();i++)
//		mesh_.mark(init_nodes_[i],init_mark_);
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
Delaunay2D<TNum, TGeomKernel>::~Delaunay2D(){
//	mesh_.unmarkAll(init_mark_);
//	mesh_.freeMark(init_mark_);
//	mesh_.unmarkAll(delete_mark_);
//	mesh_.freeMark(delete_mark_);
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::
addBoundary(const OrientedBoundary& ABoundary){
	const std::vector<gmds::math::Point >&  boundary_ref = ABoundary.getOrderedPoints();
	int nb_nodes_before = init_nodes_.size();
	/* TODO
	 * Check that we don't have twice the same segment or null lenght segment
	 */
	/* we directly create the mesh nodes */
	typename std::vector<gmds::math::Point >::const_iterator it=boundary_ref.begin();
	for(;it!=boundary_ref.end();it++)
	{
		gmds::math::Point p = *it;
		init_nodes_.push_back(mesh_.newNode(p.X(),p.Y()));
		init_points_.push_back(p);
	}
	/* we also keep in mind the pair of mesh nodes that will define mesh edges
	 *
	 */
	for(int i=nb_nodes_before+1;i<init_nodes_.size();i++){
		init_segments_.push_back(std::pair<gmds::Node,gmds::Node>(init_nodes_[i-1],init_nodes_[i]));
		init_boundary_association_[init_nodes_[i-1].getID()].push_back(init_nodes_[i].getID());
		init_boundary_association_[init_nodes_[i].getID()].push_back(init_nodes_[i-1].getID());
		TNum length = init_points_[i-1].distance(init_points_[i]);
		(*size_field)[init_nodes_[i-1].getID()]=(*size_field)[init_nodes_[i-1].getID()]+length;
		(*size_field)[init_nodes_[i].getID()]=(*size_field)[init_nodes_[i].getID()]+length;
	}
	// it remains to add the last segment
	TNum length = init_points_[nb_nodes_before].distance(init_points_.back());
	init_segments_.push_back(std::pair<gmds::Node,gmds::Node>(init_nodes_.back(),init_nodes_[nb_nodes_before]));
	init_boundary_association_[init_nodes_.back().getID()].push_back(init_nodes_[nb_nodes_before].getID());
	init_boundary_association_[init_nodes_[nb_nodes_before].getID()].push_back(init_nodes_.back().getID());
	(*size_field)[init_nodes_.back().getID()]=(*size_field)[init_nodes_.back().getID()]+length;
	(*size_field)[init_nodes_[nb_nodes_before].getID()]=(*size_field)[init_nodes_[nb_nodes_before].getID()]+length;

	for(int i=nb_nodes_before;i<init_nodes_.size();i++)
		(*size_field)[init_nodes_[i].getID()]=(*size_field)[init_nodes_[i].getID()]/2.0;

}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum, TGeomKernel>::isInCircumCircle(gmds::Face& AFace,
					  gmds::math::Point& p)
{
	std::vector<gmds::Node> nodes = AFace.get<gmds::Node>();
	if(nodes.size()!=3)
		return false;

	gmds::math::Point p1,p2,p3;

	p1.setXYZ(nodes[0].X(), nodes[0].Y(),0.0);
	p2.setXYZ(nodes[1].X(), nodes[1].Y(),0.0);
	p3.setXYZ(nodes[2].X(), nodes[2].Y(),0.0);
	TNum res = TGeomKernel<TNum>::isInCircumCircle(p1,p2,p3,p);
	return (res>0.0 && !TGeomKernel<TNum>::compare(res,0));
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum,TGeomKernel>::isInRelaxedCircumCircle(gmds::Face& AFace,
					  gmds::math::Point& p)
{
	std::vector<gmds::Node> nodes = AFace.get<gmds::Node>();
	if(nodes.size()!=3)
		return false;

	gmds::math::Point p1,p2,p3;

	TNum coord1[2]={nodes[0].X(),nodes[0].Y()};
	TNum coord2[2]={nodes[1].X(),nodes[1].Y()};
	TNum coord3[2]={nodes[2].X(),nodes[2].Y()};
	p1.setXYZ(nodes[0].X(), nodes[0].Y(),0.0);
	p2.setXYZ(nodes[1].X(), nodes[1].Y(),0.0);
	p3.setXYZ(nodes[2].X(), nodes[2].Y(),0.0);
	TNum res = TGeomKernel<TNum>::isInCircumCircle(p1,p2,p3,p);
	return (res>=0.0);
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::buildCavity(gmds::Face& AFace, gmds::math::Point& p,
								std::vector<gmds::Face>& cavity,
								const bool& relax)
{
//	gmds::Face current_face = AFace;
//	std::vector<gmds::Face> neighbours = current_face.get<gmds::Face>();
////	std::cout<<"nb voisins = "<<neighbours.size()<<std::endl;
//
//	cavity.clear();
//
//	cavity.push_back(current_face);
//
//	int cavity_mark = mesh_.getNewMark();
//	mesh_.mark(current_face,cavity_mark);
//
//	/* we look for the faces composing the cavity */
//	if(!relax){
////		std::cout<<"not relax"<<std::endl;
//		for(unsigned int i=0;i<neighbours.size();i++){
//			gmds::Face * next = neighbours[i];
//			if(next!=0 && !mesh_.isMarked(next,cavity_mark) &&
//				isInCircumCircle(next,p)){
//				cavity.push_back(next);
//				mesh_.mark(next,cavity_mark);
//				std::vector<gmds::Face*> next_neighbours = next.getFaces();
//				for(unsigned int j=0;j<next_neighbours.size();j++){
//					if(!mesh_.isMarked(next_neighbours[j],cavity_mark))
//						neighbours.push_back(next_neighbours[j]);
//				}
//			}
//		}
//	}
//	else{
////		std::cout<<" relax"<<std::endl;
//
//
//		for(unsigned int i=0;i<neighbours.size();i++){
//			gmds::Face * next = neighbours[i];
//			if(next!=0 && !mesh_.isMarked(next,cavity_mark) &&
//				isInRelaxedCircumCircle(next,p)){
//				cavity.push_back(next);
//				mesh_.mark(next,cavity_mark);
//				std::vector<gmds::Face*> next_neighbours = next.getFaces();
//				// add new triangles to check
//				for(unsigned int j=0;j<next_neighbours.size();j++){
//					if(!mesh_.isMarked(next_neighbours[j],cavity_mark))
//						neighbours.push_back(next_neighbours[j]);
//				}
//			}
//		}
//	}
//
//	for(unsigned int i=0;i<cavity.size();i++)
//		mesh_.unmark(cavity[i],cavity_mark);
//	mesh_.freeMark(cavity_mark);
//
//#ifdef _DEBUG
////	std::cout<<"Cavity Size: "<<cavity.size()<<"\n";
//#endif //_DEBUG
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::splitCavity(gmds::Node& ANodeToInsert,
								std::vector<gmds::Face>& ACavity)
{

//	std::map<gmds::Node,std::vector<gmds::Face> > node2faces;
//	std::map<gmds::Face,std::vector<gmds::Face> >	tmp_F2F;
//
//	//=========================================================
//	//PASS 1 - all faces are marked
//	//=========================================================
//	int cavity_mark = mesh_.getNewMark<gmds::Face>();
//	for(unsigned int i=0;i<ACavity.size();i++){
//		mesh_.mark(ACavity[i],cavity_mark);
//	}
//
//	//=========================================================
//	//PASS 2 - we build new faces
//	//=========================================================
//	for(unsigned int i=0;i<ACavity.size();i++){
//		Face* current_face = ACavity[i];
//		std::vector<Face*> adj_faces = current_face.getAllFaces();
//		std::vector<Node*> adj_nodes = current_face.getNodes();
//
//		for(unsigned int j=0;j<adj_faces.size();j++){
//			Face *adj_f = adj_faces[j];
//
//			// the adjacent face belongs to the cavity too
//			if(adj_f!=0 && mesh_.isMarked(adj_f,cavity_mark))
//				continue;
//
//			//we get the nodes of the edge shared by current_face and adj_f
//			Node *n1, *n2;
//			if(j==0){
//				n1 = adj_nodes[1];
//				n2 = adj_nodes[2];
//			}
//			else if(j==1){
//				n1 = adj_nodes[2];
//				n2 = adj_nodes[0];
//			}
//			else { //i==2
//				n1 = adj_nodes[0];
//				n2 = adj_nodes[1];
//			}
//
//			// a new face is created
//			Face * newFace = mesh_.newTriangle(ANodeToInsert,n1,n2);
//			node2faces[n1].push_back(newFace);
//			node2faces[n2].push_back(newFace);
//			//F2F connectivity
//			tmp_F2F[newFace].resize(3);
//			if(adj_f!=0)
//				adj_f->replaceFace(current_face,newFace);
//
//			tmp_F2F[newFace][0]=adj_f;
//		}
//	}
//
//
//	//=========================================================
//	//PASS 3 - connectivities between new faces must be built
//	//=========================================================
//	/* F2F connectivities between new faces must be build */
//	std::map<Node*,std::vector<Face*> >::iterator it_N2F=node2faces.begin();
//
//	for(;it_N2F!=node2faces.end();it_N2F++){
//			std::vector<Face*> faces = it_N2F->second;
//			if (faces.size()!=2)
//				continue;
//			Node *n  = it_N2F->first;
//			Face *f1 = faces[0];
//			Face *f2 = faces[1];
//
//			std::vector<Node*> f1_nodes = f1.getNodes();
//			std::vector<Node*> f2_nodes = f2.getNodes();
//
//			if(n==f1_nodes[1])
//				tmp_F2F[f1][2]=f2;
//			else if(n==f1_nodes[2])
//				tmp_F2F[f1][1]=f2;
//			else
//				throw GMDSException("Error in the cavity filling");
//			if(n==f2_nodes[1])
//				tmp_F2F[f2][2]=f1;
//			else if(n==f2_nodes[2])
//				tmp_F2F[f2][1]=f1;
//			else
//				throw GMDSException("Error in the cavity filling");
//	}
//
//	std::map<Face*,std::vector<Face*> >::iterator it_F2F=tmp_F2F.begin();
//	for(;it_F2F!=tmp_F2F.end();it_F2F++){
//		std::vector<Face*> ff = it_F2F->second;
//		std::vector<gmds::id> ff_id;
//		ff_id.resize(3);
//		for(unsigned int k=0;k<ff.size();k++){
//			if(ff[k])
//				ff_id[k]=ff[k].getID();
//			else
//				ff_id[k]=gmds::NullID;
//		}
//#ifdef _DEBUG
////		std::cout<<(it_F2F->first).getID()<<" -> ";
////		std::vector<Node*> nn = (it_F2F->first).getNodes();
////		std::cout<<ff_id[0]<<" "<<ff_id[1]<<" "<<ff_id[2]<<"\n";
////		std::cout<<"     -> nodes ";
////		for(unsigned int k=0; k<nn.size();k++){
////			std::cout<<"("<<nn[k].X()<<", "<<nn[k].Y()<<") ";
////		}
////		std::cout<<"\n";
//#endif//_DEBUG
//		(it_F2F->first)->setFaces(ff_id);
//	}
//
//	//=========================================================
//	//PASS 4 - connectivities between new faces must be built
//	//=========================================================
//	for(unsigned int i=0;i<ACavity.size();i++)
//	{
////		Face* fi = ACavity[i];
////		std::vector<Node*> local_nodes;
////		fi.getNodes(local_nodes);
////		for(unsigned int j=0;j<local_nodes.size();j++)
////			(*init_N2F_node_variable)[local_nodes[j].getID()].remove(fi);
//
//		mesh_.unmark(ACavity[i],cavity_mark);
//	}
//	mesh_.freeMark(cavity_mark);
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::splitCavityOld(gmds::Node& ANodeToInsert,
								std::vector<gmds::Face>& ACavity)
{

//	std::map<Node,std::vector<Face> > node2faces;
//	std::map<Face,std::vector<Face> >	tmp_F2F;
//
//	int cavity_mark = mesh_.getNewMark();
//	for(unsigned int i=0;i<ACavity.size();i++){
//		mesh_.mark(ACavity[i],cavity_mark);
//	}
//	for(unsigned int i=0;i<ACavity.size();i++){
//		Face* current_face = ACavity[i];
//		std::vector<Face*> adj_faces = current_face.getAllFaces();
//		std::vector<gmds::id> adj_faces_id = current_face.getAllFaceIDs();
//		std::vector<Node*> adj_nodes = current_face.getNodes();
//
//		for(unsigned int j=0;j<adj_faces.size();j++){
//			Face *adj_f = adj_faces[j];
//			// the adj face belongs to the cavity too
//			if(adj_f!=0 && mesh_.isMarked(adj_f,cavity_mark))
//				continue;
//
//			Node *n1, *n2;
//			if(j==0){
//				n1 = adj_nodes[1];
//				n2 = adj_nodes[2];
//			}
//			else if(j==1){
//				n1 = adj_nodes[2];
//				n2 = adj_nodes[0];
//			}
//			else { //i==2
//				n1 = adj_nodes[0];
//				n2 = adj_nodes[1];
//			}
//			Face * newFace = mesh_.newTriangle(ANodeToInsert,n1,n2);
//			(*init_N2F_node_variable)[ANodeToInsert.getID()].push_back(newFace);
//			(*init_N2F_node_variable)[n1.getID()].push_back(newFace);
//			(*init_N2F_node_variable)[n2.getID()].push_back(newFace);
//			node2faces[n1].push_back(newFace);
//			node2faces[n2].push_back(newFace);
//			//F2F connectivity
//			tmp_F2F[newFace].resize(3);
//			if(adj_f!=0)
//				adj_f->replaceFace(current_face,newFace);
//
//			tmp_F2F[newFace][0]=adj_f;
//		}
//	}
//
//	/* F2F connectivities between new faces must be build */
//	std::map<Node*,std::vector<Face*> >::iterator it_N2F=node2faces.begin();
//
//	for(;it_N2F!=node2faces.end();it_N2F++){
//			std::vector<Face*> faces = it_N2F->second;
//			if (faces.size()!=2)
//				continue;
//			Node *n  = it_N2F->first;
//			Face *f1 = faces[0];
//			Face *f2 = faces[1];
//
//			std::vector<Node*> f1_nodes = f1.getNodes();
//			std::vector<Node*> f2_nodes = f2.getNodes();
//
//			if(n==f1_nodes[1])
//				tmp_F2F[f1][2]=f2;
//			else if(n==f1_nodes[2])
//				tmp_F2F[f1][1]=f2;
//			else
//				throw GMDSException("Error in the cavity filling");
//			if(n==f2_nodes[1])
//				tmp_F2F[f2][2]=f1;
//			else if(n==f2_nodes[2])
//				tmp_F2F[f2][1]=f1;
//			else
//				throw GMDSException("Error in the cavity filling");
//	}
//
//	std::map<Face*,std::vector<Face*> >::iterator it_F2F=tmp_F2F.begin();
//
//	for(;it_F2F!=tmp_F2F.end();it_F2F++){
//		std::vector<Face*> ff = it_F2F->second;
//		std::vector<gmds::id> ff_id;
//		ff_id.resize(3);
//		for(unsigned int k=0;k<ff.size();k++){
//			if(ff[k])
//				ff_id[k]=ff[k].getID();
//			else
//				ff_id[k]=gmds::NullID;
//		}
//#ifdef _DEBUG
////		std::cout<<(it_F2F->first).getID()<<" -> ";
////		std::vector<Node*> nn = (it_F2F->first).getNodes();
////		std::cout<<ff_id[0]<<" "<<ff_id[1]<<" "<<ff_id[2]<<"\n";
////		std::cout<<"     -> nodes ";
////		for(unsigned int k=0; k<nn.size();k++){
////			std::cout<<"("<<nn[k].X()<<", "<<nn[k].Y()<<") ";
////		}
////		std::cout<<"\n";
//#endif//_DEBUG
//		(it_F2F->first)->setFaces(ff_id);
//	}
//	for(unsigned int i=0;i<ACavity.size();i++)
//	{
//		Face* fi = ACavity[i];
//		std::vector<Node*> local_nodes;
//		fi.getNodes(local_nodes);
//		for(unsigned int j=0;j<local_nodes.size();j++)
//			(*init_N2F_node_variable)[local_nodes[j].getID()].remove(fi);
//		mesh_.deleteFace(fi);
//	}
//	//TODO expensive
//	mesh_.unmarkAll(cavity_mark);
//	mesh_.freeMark(cavity_mark);
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
gmds::Face  Delaunay2D<TNum, TGeomKernel>::findFirstTriangle(gmds::math::Point& p)
{
	gmds::IGMesh::face_iterator it = mesh_.faces_begin();
	for(;!it.isDone();it.next()){
		gmds::Face f = it.value();
		if(isInCircumCircle(f,p))
			return f;
	}
}
/*--------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
gmds::Face  Delaunay2D<TNum, TGeomKernel>::
findRelaxedFirstTriangle(gmds::math::Point& p)
{
	gmds::IGMesh::face_iterator it = mesh_.faces_begin();
	for(;!it.isDone();it.next()){
		gmds::Face f = it.value();
		if(isInRelaxedCircumCircle(f,p))
			return f;
	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::displayMeshInfo(){
//	std::cout<<"##################################################"<<"\n";
//	gmds::IGMesh::face_iterator it = mesh_.faces_begin();
//	for(;!it.isDone();it.next()){
//			gmds::Face f = it.value();
//		std::cout<<"Face "<<f.getID()<<"\n";
//		std::vector<gmds::TCellID> nodeIDs = f.getAllNodeIDs();
//		std::vector<gmds::Node> nodes = f.getAllNodes();
//		std::vector<gmds::TCellID> faceIDs = f.getAllFaceIDs();
//		std::cout<<"	- node ids: ";
//		for(unsigned int i=0; i<nodeIDs.size();i++)
//			std::cout<<nodeIDs[i]<<" ";
//		std::cout<<"\n";
//		std::cout<<"	- node loc: ";
//		for(unsigned int i=0; i<nodes.size();i++)
//			std::cout<<"("<<nodes[i].X()<<", "<<nodes[i].Y()<<", "<<nodes[i].Z()<<") ";
//		std::cout<<"\n";
//		std::cout<<"	- face ids: ";
//		for(unsigned int i=0; i<faceIDs.size();i++)
//			std::cout<<faceIDs[i]<<" ";
//		std::cout<<"\n";
//	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::insertNode(gmds::Node& ANode)
{
	gmds::math::Point p(ANode.X(),ANode.Y());
	bool relax = false;
	gmds::Face f = findFirstTriangle(p);
	if(f.getID()==gmds::NullID){
		//std::cout<<"NO FIRST TRIANGLE!!!"<<"\n";
		relax=true;
		f =findRelaxedFirstTriangle(p);
	}
	if(f.getID()==gmds::NullID){
			std::cout<<"NO RELAXED FIRST TRIANGLE!!!"<<"\n";
	}
	std::vector<gmds::Face> cavity;
	buildCavity(f,p,cavity, relax);
	splitCavityOld(ANode,cavity);
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::triangulate()
{
	gmds::math::Point pmin,pmax, p_xmin_ymax, p_xmax_ymin;
	TGeomKernel<TNum>::getBoundingBox(init_points_,pmin,pmax);

	TNum diameter = sqrt((pmax.X()-pmin.X()) * (pmax.X()-pmin.X()) +
			(pmax.Y()-pmin.Y()) * (pmax.Y()-pmin.Y()));
	TNum coord[2]={pmin.X()-diameter,pmax.Y()+diameter};
	p_xmin_ymax.setXYZ(coord[0],coord[1],0.0);
	coord[0] = pmax.X()+diameter;
	coord[1] = pmin.Y()-diameter;
	p_xmax_ymin.setXYZ(coord[0],coord[1],0.0);
	coord[0] = pmin.X()-diameter;
	coord[1] = pmin.Y()-diameter;
	pmin.setXYZ(coord[0],coord[1],0.0);
	coord[0] = pmax.X()+diameter;
	coord[1] = pmax.Y()+diameter;
	pmax.setXYZ(coord[0],coord[1],0.0);
	// bounding box creation
	initTriangulation(pmin,p_xmax_ymin,pmax,p_xmin_ymax);

	// triangulation of the init_nodes_ set
	for(unsigned int i=0;i<init_nodes_.size();i++){
		gmds::Node n = init_nodes_[i];
		insertNode(n);
//		static int iNode=1;
//		std::ostringstream str;
//		str<<"insert_node_"<<iNode++;
//
//		gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(mesh_);
//		writer.write(str.str(),N|F);

	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::triangulateConstrained()
{
	/* after getting a triangulation of the set of points we need to enforce
	 * the segment preservation. To do that we need to work on mesh nodes and
	 * no more on points. Thus for each inner point, we keep in memory in
	 * the corresponding nodes in init_nodes_.
	 * This task is done in the triangulate() operation.
	 */
	triangulate();

	constraintBoundary();

	removeOuterTriangles();
}

/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::setADimensionalSize(gmds::Face* AFace)
{
//	TNum adim = 0;
//	std::vector<gmds::Node*> nodes;
//	AFace.getNodes(nodes);
//	gmds::math::Point p0(nodes[0].X(), nodes[0].Y());
//	gmds::math::Point p1(nodes[1].X(), nodes[1].Y());
//	gmds::math::Point p2(nodes[2].X(), nodes[2].Y());
//	gmds::math::Point circumCenter;
//	GeomToolKit<TNum>::computeCircumCenter(p0,p1,p2,circumCenter);
//
//
//    const TNum dx = nodes[0].X() - circumCenter.X();
//    const TNum dy = nodes[0].Y() - circumCenter.Y();
//
//
//    TNum circum_radius = sqrt((dx * dx + dy * dy).toDouble());
//
//    (*adimensional_size_field_)[AFace.getID()]= circum_radius;
//    (*circumCenter_field_)[AFace.getID()]= circumCenter;
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::initADimensionalSize()
{
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next())
//		setADimensionalSize(it->currentItem());
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::mesh(const TAlgoStrategy AStrategy)
{
	triangulate();
	std::cout<<"Initial triangulation: \n";
	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;
#ifdef _DEBUG
	gmds::LimaWriter<gmds::IGMesh> writer(mesh_);
	writer.write("after_triangulation",gmds::N|gmds::F);
#endif //_DEBUG

	constraintBoundary();

	std::cout<<"Boundary constraint: \n";
	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;

#ifdef _DEBUG
	writer.write("after_boundary_reconstruction",gmds::N|gmds::F);
#endif //_DEBUG

	//removeOuterTriangles();
	//removeOuterTrianglesWithRays();

	std::cout<<"Outer triangles removal: \n";
	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;

#ifdef _DEBUG
	writer.write("after_outer_triangles_removing",gmds::N|gmds::F);
#endif //_DEBUG

	initADimensionalSize();
	if(AStrategy==DELAUNAY_STRATEGY)
		middleEdgeInsertionOld();
	else if(AStrategy== EDGE_LENGTH_STRATEGY)
		edgeLengthBasedInsertion();

		std::cout<<"Inner node insertion: \n";
	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;

#ifdef _DEBUG
	writer.write("after_inner_nodes_insertion",gmds::N|gmds::F);
#endif //_DEBUG

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum, TGeomKernel>::meshNew(const TAlgoStrategy AStrategy)
{
//	triangulate();
//	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
//	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;
//#ifdef _DEBUG
//	gmds::VTKWriter<IGMesh> writer(mesh_);
//	writer.write("after_triangulation",N|F);
//#endif //_DEBUG
//
//	constraintBoundary();
//	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
//	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;
//
//#ifdef _DEBUG
//	writer.write("after_boundary_reconstruction",N|F);
//#endif //_DEBUG
//
//	removeOuterTriangles();
//	//removeOuterTrianglesWithRays();
//	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
//	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;
//
//#ifdef _DEBUG
//	writer.write("after_outer_triangles_removing",N|F);
//#endif //_DEBUG
//
//	initADimensionalSize();
//
//	//middleEdgeInsertion();
//
//	std::cout<<"\t nb faces = "<<mesh_.getNbFaces()<<std::endl;
//	std::cout<<"\t nb nodes = "<<mesh_.getNbNodes()<<std::endl;
//
//#ifdef _DEBUG
//	writer.write("after_inner_nodes_insertion",N|F);
//#endif //_DEBUG
//
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum,TGeomKernel>::
isBoundaryEdge(gmds::Node& AN1, gmds::Node& AN2){
//	if(mesh_.isMarked(AN1,init_mark_) && mesh_.isMarked(AN2,init_mark_))
//	{
//		std::vector<id> segmentsj =
//				init_boundary_association_[AN1.getID()];
//		if( segmentsj[0]==AN2.getID() ||segmentsj[1]==AN2.getID())
//			return true;
//	}
	return false;
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::middleEdgeInsertion()
{
//	/* simple insertion process based on the size field stored in the mesh node
//	 * All edges are split until getting the right size.
//	 * The size is defined on the boundary nodes by as the middle of the length
//	 * of the adjacent edge.
//	 * An edge is split if its length > the length of at least one of its node
//	 * When we split it, a  new node is inserted along the edge according to
//	 * a linear ponderation.
//	 */
//	/* we work on faces as we do not have edges in the model */
//	std::set<gmds::Face> sf;
//
//	/* Initialisation of the list of faces to replace */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face *f = it->currentItem();
//		std::vector<gmds::Node*> nodes;
//		f.getNodes(nodes);
//		gmds::math::Point p0(nodes[0].X(),nodes[0].Y());
//		gmds::math::Point p1(nodes[1].X(),nodes[1].Y());
//		gmds::math::Point p2(nodes[2].X(),nodes[2].Y());
//		TNum w0 = (*size_field)[nodes[0].getID()];
//		TNum w1 = (*size_field)[nodes[1].getID()];
//		TNum w2 = (*size_field)[nodes[2].getID()];
//		TNum face_size = (w0+w1+w2)/3.0;
//		TNum radius= (*adimensional_size_field_)[f.getID()];
//
//		if (face_size > 0.01*radius)
//			sf.insert(f);
//	}
//
//	while(!sf.empty()){
//		Face* f=0;
//		std::cout<<"NB BAD FACES: "<<sf.size()<<std::endl;
//		// we look for the worst face to remove
//		sf.clear();
//		it = mesh_.faces_begin();
//		for(;!it->isDone();it->next()){
//			gmds::Face *f = it->currentItem();
//			if(mesh_.isMarked(f,delete_mark_))
//				continue;
//
//			std::vector<gmds::Node*> nodes;
//			f.getNodes(nodes);
//			gmds::math::Point p0(nodes[0].X(),nodes[0].Y());
//			gmds::math::Point p1(nodes[1].X(),nodes[1].Y());
//			gmds::math::Point p2(nodes[2].X(),nodes[2].Y());
//			TNum w0 = (*size_field)[nodes[0].getID()];
//			TNum w1 = (*size_field)[nodes[1].getID()];
//			TNum w2 = (*size_field)[nodes[2].getID()];
//			TNum face_size = (w0+w1+w2)/3.0;
//			TNum radius= (*adimensional_size_field_)[f.getID()];
//
//			if (face_size > 0.01*radius)
//				sf.insert(f);
//		}
//		do{
//			f= *(sf.begin());
//			if(f!=0 && mesh_.isMarked(f,delete_mark_)){
//				std::cout<<"\t DEL"<<std::endl;
//				sf.erase(sf.begin());
//			}
//		}
//		while(f==0 && !sf.empty());
//
//		if(f!=0){
//			gmds::math::Point circumCenter = (*circumCenter_field_)[f.getID()];
//
//			std::cout<<"=== "<<circumCenter<<" ==="<<std::endl;
//
//			std::vector<gmds::Face*> cavity;
//
//			buildCavity(f,circumCenter,cavity, false);
//			std::cout<<"-> cavity size: "<<cavity.size()<<std::endl;
//
//			for(unsigned int i=0;i<cavity.size();i++){
//				mesh_.mark(cavity[i],delete_mark_);
//			}
//			gmds::Node* n_insert = mesh_.newNode(circumCenter.X(),circumCenter.Y());
//			std::vector<gmds::Node*> nodes;
//			f.getNodes(nodes);
//			gmds::math::Point p0(nodes[0].X(),nodes[0].Y());
//			gmds::math::Point p1(nodes[1].X(),nodes[1].Y());
//			gmds::math::Point p2(nodes[2].X(),nodes[2].Y());
//			TNum w0 = (*size_field)[nodes[0].getID()];
//			TNum w1 = (*size_field)[nodes[1].getID()];
//			TNum w2 = (*size_field)[nodes[2].getID()];
//			TNum face_size = (w0+w1+w2)/3.0;
//			(*size_field)[n_insert.getID()]=face_size;
//
//
//			splitCavity(n_insert,cavity);
//
//			for(unsigned int i=0;i<cavity.size();i++){
//				mesh_.deleteFace(cavity[i]);
//			}
//
//			static int i=1;
//			std::ostringstream str;
//			str<<"inner_node_"<<i++;
//
//			gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(mesh_);
//			writer.write(str.str(),N|F);
//
//		}
//		else
//			sf.clear();
//
//	}
//
//	it = mesh_.faces_begin();
//	int idel = 0;
//	for(;!it->isDone();it->next()){
//		gmds::Face *f = it->currentItem();
//		if(mesh_.isMarked(f,delete_mark_)){
//			mesh_.deleteFace(f);
//			idel++;
//		}
//	}
//
//	std::cout<<"NB DELETED = "<<idel<<std::endl;
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::middleEdgeInsertionOld()
{
//	/* simple insertion process based on the size field stored in the mesh node
//	 * All edges are split until getting the right size.
//	 * The size is defined on the boundary nodes by as the middle of the length
//	 * of the adjacent edge.
//	 * An edge is split if its length > the length of at least one of its node
//	 * When we split it, a  new node is inserted along the edge according to
//	 * a linear ponderation.
//	 */
//
//	/* we work on faces as we do not have edges in the model */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face *f = it->currentItem();
//		std::vector<gmds::Node*> nodes;
//		f.getNodes(nodes);
//		gmds::math::Point p0(nodes[0].X(),nodes[0].Y());
//		gmds::math::Point p1(nodes[1].X(),nodes[1].Y());
//		gmds::math::Point p2(nodes[2].X(),nodes[2].Y());
//		TNum w0 = (*size_field)[nodes[0].getID()];
//		TNum w1 = (*size_field)[nodes[1].getID()];
//		TNum w2 = (*size_field)[nodes[2].getID()];
//		/* we check the length of each edge comparing to the size specified on
//		 * its extrimities
//		 */
//
//		TNum length01 = p0.distance(p1);
//		TNum length12 = p1.distance(p2);
//		TNum length02 = p0.distance(p2);
//
//		gmds::math::Point p_insert;
//		TNum w_insert;
//		bool insert = true;
//		if(length01>w0 || length01>w1){
//			if((w0!=0.0) && (w1!=0.0) && !isBoundaryEdge(nodes[0],nodes[1]))
//			{
//				p_insert = w0/(w0+w1)*p0+w1/(w0+w1)*p1;
//				TNum length_i = p0.distance(p_insert);
//				w_insert = (length_i/length01)*w0+(1.0-length_i/length01)*w1;
//			}
//		}
//		else if(length02>w0 || length02>w2){
//			if((w0!=0.0) && (w2!=0.0)&& !isBoundaryEdge(nodes[0],nodes[2])){
//				p_insert = w0/(w0+w2)*p0+w2/(w0+w2)*p2;
//				TNum length_i = p0.distance(p_insert);
//				w_insert = (length_i/length02)*w0+(1.0-length_i/length02)*w2;
//			}
//		}
//		else if(length12>w1 || length12>w2){
//			if((w1!=0.0) && (w2!=0.0)&& !isBoundaryEdge(nodes[1],nodes[2])){
//				p_insert = w1/(w1+w2)*p1+w2/(w1+w2)*p2;
//				TNum length_i = p1.distance(p_insert);
//				w_insert = (length_i/length12)*w1+(1.0-length_i/length12)*w2;
//			}
//		}
//		else
//			insert =false;
//
//		if(insert){
//			std::vector<gmds::Face*> cavity;
//			buildCavity(f,p_insert,cavity, false);
//			if(cavity.size()!=1){
//				gmds::Node* n_insert = mesh_.newNode(p_insert.X(),p_insert.Y());
//				(*size_field)[n_insert.getID()]=w_insert;
//				splitCavityOld(n_insert,cavity);
//			static int i=1;
//			std::ostringstream str;
//			str<<"inner_node_"<<i++;
//
//			gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(mesh_);
//			writer.write(str.str(),N|F);
//
//			it = mesh_.faces_begin();
//			}
//		}
//
//	}

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::edgeLengthBasedInsertion()
{
//	/*
//	 * Each edge is stored by knowing 4 gmds::id(s) :
//	 * 	- its endings nodes
//	 * 	- ids adjacent faces. If there is only one adj. face, the fourth id is
//	 * 	  equal to ???
//	 */
//	typedef gmds::id edgeDS[4];
//
//	std::vector<edgeDS> wrongEdges;
//
//	/*
//	 * ALGORITHM STEPS:
//	 * 	1- We get all the edges that do not fit the expected size (too long or
//	 * 	too small with epsilon) -> add into wrongEdges
//	 * 	2 - While wrongEdges is not empty
//	 * 		 - if lenght(e)> expected(size)
//	 * 		 		we split it in two parts and cut the 2 adj. triangles
//	 * 		 		4 edges are created. Their lenghts are checked to be inserted
//	 * 		 		in wrongEdges or not
//	 * 		 - if lenght(e)< expected(size)
//	 * 		 		e is collapsed in the middle point of e. Again all modified
//	 * 		 		curves are checked to be inserted in wrongEdges or not
//	 *
//	 */
//
//	/* we work on faces as we do not have edges in the model */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face *f = it->currentItem();
//		std::vector<gmds::Node*> nodes;
//		f.getNodes(nodes);
//		gmds::math::Point p0(nodes[0].X(),nodes[0].Y());
//		gmds::math::Point p1(nodes[1].X(),nodes[1].Y());
//		gmds::math::Point p2(nodes[2].X(),nodes[2].Y());
//		TNum w0 = (*size_field)[nodes[0].getID()];
//		TNum w1 = (*size_field)[nodes[1].getID()];
//		TNum w2 = (*size_field)[nodes[2].getID()];
//		/* we check the length of each edge comparing to the size specified on
//		 * its extrimities
//		 */
//
//		TNum length01 = p0.distance(p1);
//		TNum length12 = p1.distance(p2);
//		TNum length02 = p0.distance(p2);
//
//		gmds::math::Point p_insert;
//		TNum w_insert;
//		bool insert = true;
//		if(length01>w0 || length01>w1){
//			if((w0!=0.0) && (w1!=0.0) && !isBoundaryEdge(nodes[0],nodes[1]))
//			{
//				p_insert = w0/(w0+w1)*p0+w1/(w0+w1)*p1;
//				TNum length_i = p0.distance(p_insert);
//				w_insert = (length_i/length01)*w0+(1.0-length_i/length01)*w1;
//			}
//		}
//		else if(length02>w0 || length02>w2){
//			if((w0!=0.0) && (w2!=0.0)&& !isBoundaryEdge(nodes[0],nodes[2])){
//				p_insert = w0/(w0+w2)*p0+w2/(w0+w2)*p2;
//				TNum length_i = p0.distance(p_insert);
//				w_insert = (length_i/length02)*w0+(1.0-length_i/length02)*w2;
//			}
//		}
//		else if(length12>w1 || length12>w2){
//			if((w1!=0.0) && (w2!=0.0)&& !isBoundaryEdge(nodes[1],nodes[2])){
//				p_insert = w1/(w1+w2)*p1+w2/(w1+w2)*p2;
//				TNum length_i = p1.distance(p_insert);
//				w_insert = (length_i/length12)*w1+(1.0-length_i/length12)*w2;
//			}
//		}
//		else
//			insert =false;
//
//		if(insert){
//			std::vector<gmds::Face*> cavity;
//			buildCavity(f,p_insert,cavity, false);
//			if(cavity.size()!=1){
//				gmds::Node* n_insert = mesh_.newNode(p_insert.X(),p_insert.Y());
//				(*size_field)[n_insert.getID()]=w_insert;
//				splitCavity(n_insert,cavity);
////			static int i=1;
////			std::ostringstream str;
////			str<<"inner_node_"<<i++;
////
////			gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(mesh_);
////			writer.write(str.str(),N|F);
//
//			it = mesh_.faces_begin();
//			}
//		}
//
//	}

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::constraintBoundary(){
//	std::vector<std::pair<Node*, Node* > >::iterator it = init_segments_.begin();
//	for(;it!=init_segments_.end();it++){
//		Node* n1 = it->first;
//		Node* n2 = it->second;
//		if(!checkEdgeExistence(n1,n2))
//			constraintToEdge(n1,n2);
//	}

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::removeOuterTriangles()
{
//
//	gmds::Variable<int>* color = mesh_.template newVariable<int>(GMDS_FACE,"color");
//	int current_color = 1;
//
//	bool finish=false;
//
//	while(!finish){
//		Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//		bool find_first = false;
//		std::vector<Face*> faces;
//
//		if(current_color==1){
//			//get a triangle on the boundary
//			for(;!it->isDone() && !find_first;it->next()){
//				if((it->currentItem()).getFaces().size()!=3){
//					find_first = true;
//					(*color)[(it->currentItem()).getID()]=current_color;
//					faces.push_back(it->currentItem());
//				}
//			}
//		}
//		else{
//			//get a triangle whose color is still 0
//			for(;!it->isDone() && !find_first;it->next()){
//				if((*color)[(it->currentItem()).getID()]==0){
//					find_first = true;
//					(*color)[(it->currentItem()).getID()]=current_color;
//					faces.push_back(it->currentItem());
//				}
//			}
//		}
//		if(find_first==false)
//			finish=true;
//
//		while (!faces.empty()){
//			Face *current_face = faces.back();
//			faces.pop_back();
//
//			std::vector<Node*> current_face_nodes;
//			current_face.getNodes(current_face_nodes);
//			std::vector<Face*> neighbours;
//			current_face.getAllFaces(neighbours);
//			for(int i=0;i<3;i++){
//				//foreach face we check if it is in the same domain
//				if(neighbours[i]==0)
//					continue; //no face along the ith edge
//				if((*color)[neighbours[i].getID()]==current_color)
//					continue; //this face is already in the domain
//				int j = (i+1)%3;
//				int k = (i+2)%3;
//				// the edge is not a boundary edge
//				if(!mesh_.isMarked(current_face_nodes[j],init_mark_) ||
//						!mesh_.isMarked(current_face_nodes[k],init_mark_)){
//					(*color)[neighbours[i].getID()]=current_color;
//					faces.push_back(neighbours[i]);
//					continue;
//				}
//				//the edge could be a boundary edge
//				std::vector<id> segmentsj =
//						init_boundary_association_[current_face_nodes[j].getID()];
//				if( segmentsj[0]!=current_face_nodes[k].getID() &&
//					segmentsj[1]!=current_face_nodes[k].getID())
//				{
//					faces.push_back(neighbours[i]);
//					(*color)[neighbours[i].getID()]=current_color;
//				}
//			}
//		}
//		current_color++;
//	}
//
//	/* Every triangles being adjacent to n1, n2, n3 or n4 can be removed */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face* current_face = it->currentItem();
//		if((*color)[current_face.getID()]!=2)
//			mesh_.deleteFace(current_face);
//
//	}
//
//
//	mesh_.deleteNode(bounding_box_nodes_[0]);
//	mesh_.deleteNode(bounding_box_nodes_[1]);
//	mesh_.deleteNode(bounding_box_nodes_[2]);
//	mesh_.deleteNode(bounding_box_nodes_[3]);
//	bounding_box_nodes_.clear();
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::removeOuterTrianglesNew(){

//	gmds::Variable<int>* color = mesh_.template newVariable<int>(GMDS_FACE,"color");
//	int current_color = 1;
//
//
//	Node *n = bounding_box_nodes_[0];
//
//	bool finish=false;
//
//
//	while(!finish){
//		Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//		bool find_first = false;
//		std::vector<Face*> faces;
//
//		if(current_color==1){
//			//get a triangle on the boundary
//			for(;!it->isDone() && !find_first;it->next()){
//				if((it->currentItem()).getFaces().size()!=3){
//					find_first = true;
//					(*color)[(it->currentItem()).getID()]=current_color;
//					faces.push_back(it->currentItem());
//				}
//			}
//		}
//		else{
//			//get a triangle whose color is still 0
//			for(;!it->isDone() && !find_first;it->next()){
//				if((*color)[(it->currentItem()).getID()]==0){
//					find_first = true;
//					(*color)[(it->currentItem()).getID()]=current_color;
//					faces.push_back(it->currentItem());
//				}
//			}
//		}
//		if(find_first==false)
//			finish=true;
//
//		while (!faces.empty()){
//			Face *current_face = faces.back();
//			faces.pop_back();
//
//			std::vector<Node*> current_face_nodes;
//			current_face.getNodes(current_face_nodes);
//			std::vector<Face*> neighbours;
//			current_face.getAllFaces(neighbours);
//			for(int i=0;i<3;i++){
//				//foreach face we check if it is in the same domain
//				if(neighbours[i]==0)
//					continue; //no face along the ith edge
//				if((*color)[neighbours[i].getID()]==current_color)
//					continue; //this face is already in the domain
//				int j = (i+1)%3;
//				int k = (i+2)%3;
//				// the edge is not a boundary edge
//				if(!mesh_.isMarked(current_face_nodes[j],init_mark_) ||
//						!mesh_.isMarked(current_face_nodes[k],init_mark_)){
//					(*color)[neighbours[i].getID()]=current_color;
//					faces.push_back(neighbours[i]);
//					continue;
//				}
//				//the edge could be a boundary edge
//				std::vector<id> segmentsj =
//						init_boundary_association_[current_face_nodes[j].getID()];
//				if( segmentsj[0]!=current_face_nodes[k].getID() &&
//					segmentsj[1]!=current_face_nodes[k].getID())
//				{
//					faces.push_back(neighbours[i]);
//					(*color)[neighbours[i].getID()]=current_color;
//				}
//			}
//		}
//		current_color++;
//	}
//
//	/* Every triangles being adjacent to n1, n2, n3 or n4 can be removed */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face* current_face = it->currentItem();
//		if((*color)[current_face.getID()]!=2)
//			mesh_.deleteFace(current_face);
//
//	}
//
//
//	mesh_.deleteNode(bounding_box_nodes_[0]);
//	mesh_.deleteNode(bounding_box_nodes_[1]);
//	mesh_.deleteNode(bounding_box_nodes_[2]);
//	mesh_.deleteNode(bounding_box_nodes_[3]);
//	bounding_box_nodes_.clear();
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::removeOuterTrianglesWithRays(){

//	// the fourth first nodes are those of the extra bounding box
//	Node* n1 = bounding_box_nodes_[0];
//	Node* n2 = bounding_box_nodes_[1];
//	Node* n3 = bounding_box_nodes_[2];
//	Node* n4 = bounding_box_nodes_[3];
//
//	/* Every triangles being adjacent to n1, n2, n3 or n4 can be removed */
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//	for(;!it->isDone();it->next()){
//		gmds::Face* current_face = it->currentItem();
//		std::vector<Node*> nodes_of_face;
//		current_face.getNodes(nodes_of_face);
//		bool toDelete=false;
//		for(unsigned int i=0;i<nodes_of_face.size() && !toDelete;i++){
//			gmds::Node* ni =nodes_of_face[i];
//			if (ni==n1 || ni==n2 ||ni==n3 || ni==n4)
//				toDelete=true;
//		}
//		if(toDelete){
//			std::vector<gmds::Face*> faces_to_update = current_face.getFaces();
//			id current_id = current_face.getID();
//			for(unsigned int i=0;i<faces_to_update.size();i++) {
//				faces_to_update[i]->removeFace(current_id);
//			}
//			mesh_.deleteFace(current_face);
//		}
//	}
//
//	/* Alle the remaining triangles are built on the boundary nodes. However
//	 * some triangles are inside and other outside of the domain to triangulate.
//	 * To find the outside ones, we throw infinite rays from the center of mass of each
//	 * triangle in order to get the number of intersection of this ray with the
//	 * boundary of the domain to mesh.
//	 */
//
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator itf = mesh_.faces_begin();
//	for(;!itf->isDone();itf->next()){
//		gmds::Face* current_face = itf->currentItem();
//
//		std::vector<Node*> nodes_of_face;
//		current_face.getNodes(nodes_of_face);
//		TNum x= 	(nodes_of_face[0].X() + nodes_of_face[1].X() +
//					nodes_of_face[2].X())/3.0;
//		TNum y= 	(nodes_of_face[0].Y() + nodes_of_face[1].Y() +
//							nodes_of_face[2].Y())/3.0;
//
//		if(isOutOfTheDomain(x,y)){
//			std::vector<gmds::Face*> faces_to_update = current_face.getFaces();
//			id current_id = current_face.getID();
//			for(unsigned int i=0;i<faces_to_update.size();i++) {
//				faces_to_update[i]->removeFace(current_id);
//			}
//			mesh_.deleteFace(current_face);
//		}
//	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum,TGeomKernel>::isOutOfTheDomain(const TNum& AX, const TNum& AY){
//	TNum max_ray = 1000000.0;
//	GEPETO::Segment<2,TNum> ray1(gmds::math::Point(AX,AY),
//						gmds::math::Point(AX+max_ray,AY));
//	GEPETO::Segment<2,TNum> ray2(gmds::math::Point(AX,AY),
//						gmds::math::Point(AX,AY+max_ray));
//
//	int nb_intersections1=0, nb_intersections2=0;
//
//	for(std::vector<std::pair<Node*, Node* > >::iterator it=init_segments_.begin();
//			it!=init_segments_.end();it++){
//		Node* n1 = it->first;
//		Node* n2 = it->second;
//		gmds::math::Point p1(n1.X(),n1.Y());
//		gmds::math::Point p2(n2.X(),n2.Y());
//		GEPETO::Segment<2,TNum> s(p1,p2);
//
//		if(s.intersect(ray1)==GEPETO::GEOM_YES){
//			nb_intersections1++;
//		}
//		if(s.intersect(ray2)==GEPETO::GEOM_YES)
//		{	nb_intersections2++;
//
//		}
//	}
//
//	if (nb_intersections1%2==nb_intersections2%2)
//		return (nb_intersections1%2==0);
//	else{
//		//we throw a third ray
//
//		GEPETO::Segment<2,TNum> ray3(gmds::math::Point(AX,AY),
//							gmds::math::Point(AX,AY-max_ray));
//
//		int nb_intersections3=0;
//
//		for(std::vector<std::pair<Node*, Node* > >::iterator it=init_segments_.begin();
//				it!=init_segments_.end();it++){
//			Node* n1 = it->first;
//			Node* n2 = it->second;
//			gmds::math::Point p1(n1.X(),n1.Y());
//			gmds::math::Point p2(n2.X(),n2.Y());
//			GEPETO::Segment<2,TNum> s(p1,p2);
//			if(s.intersect(ray3)==GEPETO::GEOM_YES)
//				nb_intersections3++;
//		}
//		return (nb_intersections3%2==0);
//	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum,TGeomKernel>::checkEdgeExistence(gmds::Node& AN1, gmds::Node& AN2){
//	std::list<Face*> faces_AN1 = (*init_N2F_node_variable)[AN1.getID()];
//
//	for(std::list<Face*>::iterator it=faces_AN1.begin();it!=faces_AN1.end();it++)
//	{
//		Face* current_face = *it;
//		std::vector<Node*> nodes_of_face;
//		current_face.getNodes(nodes_of_face);
//		for(unsigned int i=0;i<nodes_of_face.size();i++)
//			if(nodes_of_face[i]==AN2)
//				return true;
//	}
	return false;

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::initTriangulation(
		gmds::math::Point& AP_xmin_ymin,
		gmds::math::Point& AP_xmax_ymin,
		gmds::math::Point& AP_xmax_ymax,
		gmds::math::Point& AP_xmin_ymax)
{
//	gmds::Node n1 = mesh_.newNode(AP_xmin_ymin.X(),AP_xmin_ymin.Y(),0.0);
//	gmds::Node n2 = mesh_.newNode(AP_xmax_ymin.X(),AP_xmax_ymin.Y(),0.0);
//	gmds::Node n3 = mesh_.newNode(AP_xmax_ymax.X(),AP_xmax_ymax.Y(),0.0);
//	gmds::Node n4 = mesh_.newNode(AP_xmin_ymax.X(),AP_xmin_ymax.Y(),0.0);
//
//	bounding_box_nodes_.push_back(n1);
//	bounding_box_nodes_.push_back(n2);
//	bounding_box_nodes_.push_back(n3);
//	bounding_box_nodes_.push_back(n4);
//
//	gmds::Face t1 = mesh_.newTriangle(n1,n2,n3);
//	(*init_N2F_node_variable)[n1.getID()].push_back(t1);
//	(*init_N2F_node_variable)[n2.getID()].push_back(t1);
//	(*init_N2F_node_variable)[n3.getID()].push_back(t1);
//	gmds::Face t2 = mesh_.newTriangle(n3,n4,n1);
//	(*init_N2F_node_variable)[n3.getID()].push_back(t2);
//	(*init_N2F_node_variable)[n4.getID()].push_back(t2);
//	(*init_N2F_node_variable)[n1.getID()].push_back(t2);
//	std::vector<gmds::TCellID> f2f;
//	f2f.resize(3);
//	f2f[0]=gmds::NullID;
//	f2f[1]=t2.getID();
//	f2f[2]=gmds::NullID;
//	t1.set<gmds::Face>(f2f);
//
//	f2f[1]=t1.getID();
//	t2.set<gmds::Face>(f2f);
//
//	gmds::IGMeshDoctor doc(&mesh_);
//	doc.orientFaces();
//

}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
bool Delaunay2D<TNum,TGeomKernel>::
swap(gmds::Face& AF1, gmds::Face& AF2, const bool& AB)
{
//	std::vector<gmds::id> f1_face_ids, f2_face_ids;
//	AF1.getAllFaceIDs(f1_face_ids);
//	AF2.getAllFaceIDs(f2_face_ids);
//
//	std::vector<gmds::Face*> f1_faces, f2_faces;
//	AF1.getAllFaces(f1_faces);
//	AF2.getAllFaces(f2_faces);
//
//	//displayMeshInfo();
//	int f2_from_f1, f1_from_f2;
//
//	if(f1_face_ids[0]==AF2.getID())
//		f2_from_f1 = 0;
//	else if(f1_face_ids[1]==AF2.getID())
//		f2_from_f1 = 1;
//	else
//		f2_from_f1 = 2;
//
//	if(f2_face_ids[0]==AF1.getID())
//		f1_from_f2=0;
//	else if(f2_face_ids[1]==AF1.getID())
//		f1_from_f2=1;
//	else
//		f1_from_f2 = 2;
//
//	//before swapping, we check if the swapping is possible (non convex area)
//	std::vector<gmds::id> f1_node_ids, f2_node_ids;
//	AF1.getAllNodeIDs(f1_node_ids);
//	AF2.getAllNodeIDs(f2_node_ids);
//
//	if(!AB){
//		std::vector<gmds::Node*> f1_nodes, f2_nodes;
//		AF1.getAllNodes(f1_nodes);
//		AF2.getAllNodes(f2_nodes);
//
//
//		gmds::math::Point p1(	f1_nodes[f2_from_f1].X(),
//				f1_nodes[f2_from_f1].Y());
//		gmds::math::Point p2(	f2_nodes[f1_from_f2].X(),
//				f2_nodes[f1_from_f2].Y());
//		gmds::math::Point p3(	f1_nodes[(f2_from_f1+1)%3].X(),
//				f1_nodes[(f2_from_f1+1)%3].Y());
//		gmds::math::Point p4(	f1_nodes[(f2_from_f1+2)%3].X(),
//				f1_nodes[(f2_from_f1+2)%3].Y());
//
//		GEPETO::EGeomPredicateResult r1 =p3.isStrictlyOnLeft(p1,p2); // TODO Not strict
//		GEPETO::EGeomPredicateResult r2 =p4.isStrictlyOnLeft(p1,p2); //idem
//		GEPETO::Segment<2,TNum> seg(p1,p2);
//		if (r1==r2 ||seg.isIn(p3) ||seg.isIn(p4)  )
//			return false;
//	}
//	//update the F2N connection
//	std::vector<gmds::id> new_f1_node_ids;
//	new_f1_node_ids.resize(3);
//	new_f1_node_ids[0] = f1_node_ids[(f2_from_f1+2)%3];
//	new_f1_node_ids[1] = f1_node_ids[f2_from_f1];
//	new_f1_node_ids[2] = f2_node_ids[f1_from_f2];
//
//	std::vector<gmds::id> new_f2_node_ids;
//	new_f2_node_ids.resize(3);
//	new_f2_node_ids[0] = f2_node_ids[(f1_from_f2+2)%3];
//	new_f2_node_ids[1] = f2_node_ids[f1_from_f2];
//	new_f2_node_ids[2] = f1_node_ids[f2_from_f1];
//
//	AF1->setNodes(new_f1_node_ids);
//	AF2->setNodes(new_f2_node_ids);
//
//	//update the F2F connection
//	std::vector<gmds::id> new_f1_face_ids;
//
//	new_f1_face_ids.resize(3);
//	new_f1_face_ids[0] = f1_face_ids[f2_from_f1];
//
//	new_f1_face_ids[1] = f2_face_ids[(f1_from_f2+2)%3];
//	// if the face exists, it must known its new neighbour
//	if(f2_face_ids[(f1_from_f2+2)%3]!=gmds::NullID)
//		f2_faces[(f1_from_f2+2)%3]->replaceFace(AF2,AF1);
//	new_f1_face_ids[2] = f1_face_ids[(f2_from_f1+1)%3];
//
//	std::vector<gmds::id> new_f2_face_ids;
//	new_f2_face_ids.resize(3);
//	new_f2_face_ids[0] = f2_face_ids[f1_from_f2];
//	new_f2_face_ids[1] = f1_face_ids[(f2_from_f1+2)%3];
//	// if the face exists, it must known its new neighbour
//	if(f1_face_ids[(f2_from_f1+2)%3]!=gmds::NullID)
//		f1_faces[(f2_from_f1+2)%3]->replaceFace(AF1,AF2);
//	new_f2_face_ids[2] = f2_face_ids[(f1_from_f2+1)%3];
//
//	AF1->setFaces(new_f1_face_ids);
//	AF2->setFaces(new_f2_face_ids);

	//displayMeshInfo();
//	static int i=1;
//	std::ostringstream str;
//	str<<"swap_"<<i++;
//
//	gmds::VTKWriter<DIM2|N|F|F2N|F2F> writer(mesh_);
//	writer.write(str.str(),N|F);

	return true;
}

/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::
getTriangleStrip(gmds::Node& AN1, gmds::Node& AN2,std::vector<gmds::Face>& AStrip)
{
//	AStrip.clear();
//	gmds::math::Point p1(AN1.X(), AN1.Y());
//	gmds::math::Point p2(AN2.X(), AN2.Y());
//	GEPETO::Segment<2,TNum> segment(p1,p2);
//
//	Mesh<DIM2|N|F|F2N|F2F>::faces_iterator it = mesh_.faces_begin();
//
//	for(;!it->isDone();it->next()){
//		Face* f = it->currentItem();
//		std::vector<gmds::Node*> nodes;
//		f.getNodes(nodes);
//		// we only consider triangles
//		if (nodes.size()!=3)
//			continue;
//
//		gmds::math::Point f_p0(nodes[0].X(), nodes[0].Y());
//		gmds::math::Point f_p1(nodes[1].X(), nodes[1].Y());
//		gmds::math::Point f_p2(nodes[2].X(), nodes[2].Y());
//		GEPETO::Triangle<2,TNum> tri(f_p0,f_p1,f_p2);
//		if (tri.intersect(segment,true)==GEPETO::GEOM_YES)
//				AStrip.push_back(f);
//	}
}
/*----------------------------------------------------------------------------*/
template<typename TNum, template<typename T> class TGeomKernel>
void Delaunay2D<TNum,TGeomKernel>::
constraintToEdge(gmds::Node& AN1, gmds::Node& AN2)
{
//	std::vector<gmds::Face*> strip;
// 	getTriangleStrip(AN1,AN2,strip);
// 	if(strip.size()==0)
// 		return;
//
//
//	gmds::math::Point p1(AN1.X(), AN1.Y());
//	gmds::math::Point p2(AN2.X(), AN2.Y());
//	GEPETO::Segment<2,TNum> segment(p1,p2);
//
//
// 	int nbTriangles = strip.size();
// 	gmds::Face* left_face=0, *right_face=0;
// 	int strip_mark = mesh_.getNewMark();
// 	int done_mark  = mesh_.getNewMark();
//
// 	/* we get the triangles containing AN1 and AN2 and we mark all the
// 	 * triangles of the strip with strip_mark */
// 	int id1 = AN1.getID(), id2 = AN2.getID();
// 	for(unsigned int i=0;i<nbTriangles;i++)
// 	{
// 		mesh_.mark(strip[i],strip_mark);
// 		std::vector<gmds::id> ids;
// 		strip[i].getNodeIDs(ids);
// 		if(left_face==0 && (ids[0]==id1 || ids[1]==id1 || ids[2]==id1) )
// 			left_face = strip[i];
// 		else if(right_face==0 && (ids[0]==id2 || ids[1]==id2 || ids[2]==id2) )
// 			right_face = strip[i];
// 	}
//
// 	// alternative swapping from left and right are performed until building edge
// 	// [AN1,AN2]
// 	bool edge12_done = false;
// 	bool left = true;
// 	while(!edge12_done){
// 		gmds::Face* current_face = (left)?left_face:right_face;
// 		std::vector<gmds::Face*> adj_faces;
// 		current_face.getFaces(adj_faces);
// 		gmds::Face* next_face=0;
// 		int i=0;
// 		while(next_face==0 && i<adj_faces.size()){
// 			gmds::Face* adj_face_i = adj_faces[i];
// 			if(mesh_.isMarked(adj_face_i,strip_mark)&& !mesh_.isMarked(adj_face_i,done_mark))
// 				next_face = adj_face_i;
// 			i++;
// 		}
//
// 		if (next_face==0 && i==adj_faces.size())
// 			edge12_done=true;
// 		else if(swap(current_face,next_face)){
// 			/* we have to check which face between current and next face is still
// 			 * intersected by [AN1,AN2] */
// 			std::vector<gmds::Node*> nodes;
// 			current_face.getNodes(nodes);
//
// 			gmds::math::Point f_p0(nodes[0].X(), nodes[0].Y());
// 			gmds::math::Point f_p1(nodes[1].X(), nodes[1].Y());
// 			gmds::math::Point f_p2(nodes[2].X(), nodes[2].Y());
// 			GEPETO::Triangle<2,TNum> tri(f_p0,f_p1,f_p2);
// 			if (segment.intersect(tri,true)==GEPETO::GEOM_YES)
//			{
// 				if (left)
// 					left_face = current_face;
// 				else
// 					right_face = current_face;
//
// 				mesh_.mark(next_face,done_mark);
//			}
// 			else
//			{
// 				if (left)
// 					left_face = next_face;
// 				else
// 					right_face = next_face;
// 				mesh_.mark(current_face,done_mark);
//			}
// 		}
// 		left=!left;
// 	}
//
// 	// we clean the strip_mark
// 	for(unsigned int i=0;i<nbTriangles;i++){
// 		mesh_.unmark(strip[i],strip_mark);
// 		mesh_.unmark(strip[i],done_mark );
// 	}
// 	mesh_.freeMark(strip_mark);
// 	mesh_.freeMark(done_mark );

}
/*----------------------------------------------------------------------------*/

