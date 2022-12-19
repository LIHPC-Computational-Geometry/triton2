/*--------------------------------------------------------------------------*/
/*
 * GeomToolKit.t.h
 *
 *  Created on: 19 janv. 2011
 *      Author: ledouxf
 */
/*--------------------------------------------------------------------------*/
template<typename T>
bool GeomToolKit<T>::compare(const T& a, const T&b,
						  const T& tolerance)
{
	return (a>b-tolerance && a<b+tolerance);
}
/*--------------------------------------------------------------------------*/
template<typename T>
void GeomToolKit<T>::getBoundingBox(
		const std::vector<gmds::math::Point >& APnts,
		gmds::math::Point& APmin,
		gmds::math::Point& APmax,
		const T& ATol)
{
	T xmin, xmax, ymin, ymax;
	xmin = xmax = APnts[0].X();
	ymin = ymax = APnts[0].Y();

	const unsigned int nbPnts = APnts.size();
	for(unsigned int i=1;i<nbPnts;i++){
		if(APnts[i].X()<xmin)
			xmin = APnts[i].X();
		if (APnts[i].X()>xmax)
			xmax = APnts[i].X();

		if(APnts[i].Y()<ymin)
			ymin = APnts[i].Y();
		if (APnts[i].Y()>ymax)
			ymax = APnts[i].Y();
	}

	T min[2] = {xmin-ATol,ymin-ATol};
	T max[2] = {xmax+ATol,ymax+ATol};
	APmin.setXYZ(min[0],min[1],0.0);
	APmax.setXYZ(max[0], max[1], 0.0);
}
/*--------------------------------------------------------------------------*/
template<typename T>
T GeomToolKit<T>::isInCircumCircle(
		gmds::math::Point& p1, gmds::math::Point& p2,
		gmds::math::Point& p3, gmds::math::Point& n)
{
	gmds::math::Matrix<4,4> mat;
	T nx = n.X();
	T ny = n.Y();
	T p1x = p1.X();
	T p1y = p1.Y();
	T p2x = p2.X();
	T p2y = p2.Y();
	T p3x = p3.X();
	T p3y = p3.Y();

	T val[4][4] ={ {p1x,p1y,p1x*p1x+p1y*p1y,1},
						{p2x,p2y,p2x*p2x+p2y*p2y,1},
						{p3x,p3y,p3x*p3x+p3y*p3y,1},
						{nx ,ny ,nx*nx+ny*ny    ,1}
						};
	mat.set(val);
#ifdef _DEBUG
//
//	bool i = mat.det()>0;
////	if(i!=0){
//		std::cout<<"("<<nx<<", "<<ny<<") in ("<<p1x<<", "<<p1y<<") ("
//				<<p2x<<", "<<p2y<<") "<<"("<<p3x<<", "<<p3y<<")"<<"\n";
//		std::cout<<"    "<<i<<" with det "<<mat.det()<<"\n";
//
////	}
#endif //_DEBUG
	return mat.det();

}
/*--------------------------------------------------------------------------*/
template<typename T>
void GeomToolKit<T>::computeCircumCenter(
		const gmds::math::Point& p1, const gmds::math::Point& p2,
		const gmds::math::Point& p3, gmds::math::Point& center){
	T d, a1, a2, a3;

	const T x1 = p1.X();
	const T x2 = p2.X();
	const T x3 = p3.X();
	const T y1 = p1.Y();
	const T y2 = p2.Y();
	const T y3 = p3.Y();

	d = 2. * (T)(y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2));
	if(d == 0.0) {
		// Msg::Warning("Colinear points in circum circle computation");
		center.setXYZ(-99999., -99999.,0);
		return ;
	}

	a1 = x1 * x1 + y1 * y1;
	a2 = x2 * x2 + y2 * y2;
	a3 = x3 * x3 + y3 * y3;
	center.setXYZ((T)((a1 * (y3 - y2) + a2 * (y1 - y3) + a3 * (y2 - y1)) / d),
				 (T)((a1 * (x2 - x3) + a2 * (x3 - x1) + a3 * (x1 - x2)) / d),
				 (T)(0.0));
}


/*--------------------------------------------------------------------------*/

