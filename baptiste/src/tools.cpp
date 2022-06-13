#include <gmds/baptiste/tools.h>
#include "gmds/io/IGMeshIOService.h"
#include <gmds/io/VTKReader.h>
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/baptiste/RLBlockSet.h"

#include "gmds/igalgo/r2d.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/math/Quadrilateral.h>
#include <gmds/math/Triangle.h>

using namespace gmds;

std::vector<double> gmds::LinearSpacedArray(double a, double b, std::size_t N)
{
	double h = (b - a) / static_cast<double>(N);
	std::vector<double> v(N+1);
	std::vector<double>::iterator x;
	double val;
	for (x = v.begin(), val = a; x != v.end(); ++x, val += h)
	{
		*x = val;
	}
	return v;
}

void gmds::cloneMesh(const Mesh &originalMesh, Mesh &newMesh)
{
	if (originalMesh.getDim() == 3 or newMesh.getDim() == 3)
	{
		throw GMDSException("Dimension must be 2");
	}

	for (int nodeID:originalMesh.nodes())
	{
		newMesh.newNode(originalMesh.get<Node>(nodeID).point());
	}
	std::set<TCellID> invalidIndexes;
	for (int faceID = 0; faceID <= originalMesh.getMaxLocalID(2); faceID++)
	{
		if (originalMesh.has<Face>(faceID))
		{
			newMesh.newFace(originalMesh.get<Face>(faceID).get<Node>());
		}
		else
		{
			newMesh.newFace({0, 0, 0, 0});
			invalidIndexes.insert(faceID);
		}
	}
	for (int faceID:invalidIndexes)
	{
		newMesh.deleteFace(faceID);
	}
}

void gmds::cloneBlockSet(const RLBlockSet &originalBlockSet, RLBlockSet &newBlockSet)
{
	cloneMesh(originalBlockSet.m_mesh, newBlockSet.m_mesh);
	newBlockSet.xSize = originalBlockSet.xSize;
	newBlockSet.ySize = originalBlockSet.ySize;
}


Mesh gmds::readMesh(std::string filename)
{
	Mesh mesh = Mesh(MeshModel(DIM2|F|N|F2N));
	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
	return mesh;
}

std::vector<Action *> gmds::getActionsVector()
{
	std::vector<Action *> actions;
	for (int iV = 0; iV <= 1; iV++)
	{
		for (int iAxis = 0; iAxis <= 1; iAxis++) {
			std::vector<int> ranges = {-2, -1, 1, 2};
			for (int iRange : ranges)
			{
				ActionEdit e = ActionEdit(iV, iAxis, iRange);
				actions.push_back(&e);
			}
		}
	}
	ActionDelete d = ActionDelete();
	//actions.push_back(&d);
	return actions;
}

double gmds::getMeshArea(Mesh &mesh)
{
	double res = 0;
	for (int faceID : mesh.faces())
	{
		Face face = mesh.get<Face>(faceID);
		res += face.area();
	}
	return res;
}

void gmds::volfraccomputation_2d_reverse(const gmds::Mesh *AMesh, gmds::Mesh *AImprintMesh, gmds::Variable<double> *AVolFrac)
{
	// check validity of the inputs
	bool valid_input = true;
	std::string msg("volfraccomputation_2d ");

	gmds::MeshModel model = AMesh->getModel();
	if (!model.has(gmds::F) || !model.has(gmds::F2N) || !model.has(gmds::DIM2))
	{
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AMesh->getNbFaces() != AMesh->getNbQuadrilaterals())
	{
		msg += std::string("AMesh should have only quads");
		valid_input = false;
	}
	gmds::MeshModel modelImprint = AImprintMesh->getModel();
	if (!modelImprint.has(gmds::F) || !modelImprint.has(gmds::F2N) || !model.has(gmds::DIM2))
	{
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}

	if(AImprintMesh->getNbFaces() != AImprintMesh->getNbTriangles())
	{
		msg += std::string("AImprintMesh should have only triangles");
		valid_input = false;
	}

	// check mesh orientation
	// TODO
	for(auto f_id: AMesh->faces())
	{
		gmds::Face f = AMesh->get<Face>(f_id);
		std::vector<gmds::Node> n = f.get<gmds::Node>();
		AVolFrac->set(f_id,0);

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		//std::cout<<"SJ Face "<<f_id<<"("<<n[0].id()<<", "<<n[1].id()<<", "<<n[2].id()<<", "<<n[3].id()<<": "<<sj<<std::endl;
		if(sj < 0)
		{
			msg += std::string("AMesh has a bad cell.");
			valid_input = false;
			break;
		}
	}

	for (auto f_id: AMesh->faces())
	{

		gmds::Face f = AMesh->get<Face>(f_id);

		r2d_rvec2 verticesQuad[4];
		std::vector<gmds::Node> n_quad = f.get<gmds::Node>();

		verticesQuad[0].x = n_quad[0].X();
		verticesQuad[0].y = n_quad[0].Y();
		verticesQuad[1].x = n_quad[1].X();
		verticesQuad[1].y = n_quad[1].Y();
		verticesQuad[2].x = n_quad[2].X();
		verticesQuad[2].y = n_quad[2].Y();
		verticesQuad[3].x = n_quad[3].X();
		verticesQuad[3].y = n_quad[3].Y();

		r2d_int numverts = 4;
		r2d_plane planes[4];
		r2d_poly_faces_from_verts(planes,  verticesQuad, numverts);

		for(auto tri_id: AImprintMesh->faces())
		{
			gmds::Face tri = AImprintMesh->get<Face>(tri_id);
			std::vector<gmds::Node> n_tri = tri.get<gmds::Node>();

			gmds::TCoord xyz[3][3];

			xyz[0][0] = n_tri[0].X();
			xyz[0][1] = n_tri[0].Y();
			xyz[0][2] = n_tri[0].Z();
			xyz[1][0] = n_tri[1].X();
			xyz[1][1] = n_tri[1].Y();
			xyz[1][2] = n_tri[1].Z();
			xyz[2][0] = n_tri[2].X();
			xyz[2][1] = n_tri[2].Y();
			xyz[2][2] = n_tri[2].Z();

			gmds::math::Point pt0(xyz[0][0], xyz[0][1], xyz[0][2]);
			gmds::math::Point pt1(xyz[1][0], xyz[1][1], xyz[1][2]);
			gmds::math::Point pt2(xyz[2][0], xyz[2][1], xyz[2][2]);
			gmds::math::Triangle t(pt0, pt1, pt2);

			const double volsurf = t.area();

			r2d_poly poly;
			r2d_rvec2 verts[3];

			verts[0].x = xyz[0][0];
			verts[0].y = xyz[0][1];
			verts[1].x = xyz[1][0];
			verts[1].y = xyz[1][1];
			verts[2].x = xyz[2][0];
			verts[2].y = xyz[2][1];

			r2d_init_poly(&poly, verts, 3);

			r2d_clip(&poly, planes,4);
			r2d_int POLY_ORDER = 2;
			r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
			r2d_reduce(&poly, om, POLY_ORDER);

			double vf_tri = AVolFrac->value(tri_id);

			std::cout<<"Id Face : "<<tri_id<<std::endl;
			std::cout<<"valeur om : "<<om<<std::endl;
			std::cout<<"valeur om[0] : "<<om[0]<<std::endl;
			std::cout<<"valeur IoU : "<<volsurf<<std::endl;
			std::cout<<"old value : "<<vf_tri<<std::endl;

			AVolFrac->set(tri_id, vf_tri + om[0]/volsurf);

			std::cout<<"valeur volfrac : "<<AVolFrac->value(tri_id)<<std::endl;

		}
	}
}

void gmds::anotherVolFrac(const gmds::Mesh *AMesh, gmds::Mesh *AImprintMesh, gmds::Variable<double> *AVolFrac)
{
	// check validity of the inputs
	bool valid_input = true;
	std::string msg("volfraccomputation_2d ");

	gmds::MeshModel model = AMesh->getModel();
	if (!model.has(gmds::F) || !model.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}
	if(AMesh->getNbFaces() != AMesh->getNbQuadrilaterals()) {
		msg += std::string("AMesh should have only quads");
		valid_input = false;
	}
	gmds::MeshModel modelImprint = AImprintMesh->getModel();
	if (!modelImprint.has(gmds::F) || !modelImprint.has(gmds::F2N) || !model.has(gmds::DIM2)) {
		msg += std::string("bad model for AMesh");
		valid_input = false;
	}
	if(AImprintMesh->getNbFaces() != AImprintMesh->getNbTriangles()) {
		msg += std::string("AImprintMesh should have only triangles");
		valid_input = false;
	}

	// check mesh orientation
	// TODO
	for(auto f_id: AMesh->faces()) {
		AVolFrac->set(f_id, 0);
		gmds::Face f = AMesh->get<Face>(f_id);
		std::vector<gmds::Node> n = f.get<gmds::Node>();

		gmds::math::Quadrilateral quad(n[0].point(), n[1].point(), n[2].point(), n[3].point());
		double sj = quad.computeScaledJacobian2D();
		if(sj < 0) {
			msg += std::string("AMesh has a bad cell.");
			valid_input = false;
			break;
		}
	}

	//	for(auto f_id: AImprintMesh->faces()) {
	//		gmds::Face f = AImprintMesh->get<Face>(f_id);
	//		std::vector<gmds::Node> n = f.get<gmds::Node>();
	//
	//		gmds::math::Triangle tri(n[0].point(), n[1].point(), n[2].point());
	//		double sj = tri.computeScaledJacobian2D();
	//		if(sj < 0) {
	//			msg += std::string("AImprintMesh has a bad cell.");
	//			valid_input = false;
	//			break;
	//		}
	//	}

	if(!valid_input) {
		throw gmds::GMDSException(msg);
	}

	for(auto tri_id: AImprintMesh->faces()) {
		gmds::Face tri = AImprintMesh->get<Face>(tri_id);
		std::vector<gmds::Node> n_tri = tri.get<gmds::Node>();

		r2d_rvec2 vertices[3];
		vertices[0].x = n_tri[0].X();
		vertices[0].y = n_tri[0].Y();
		vertices[1].x = n_tri[2].X();
		vertices[1].y = n_tri[2].Y();
		vertices[2].x = n_tri[1].X();
		vertices[2].y = n_tri[1].Y();
		r2d_int numverts = 3;
		r2d_plane planes[3];
		r2d_poly_faces_from_verts(planes,  vertices, numverts);

		for(auto f_id: AMesh->faces()) {
			gmds::Face f = AMesh->get<Face>(f_id);

			std::vector<gmds::Node> n_quad = f.get<gmds::Node>();
			gmds::TCoord xyz[4][3];
			xyz[0][0] = n_quad[0].X();
			xyz[0][1] = n_quad[0].Y();
			xyz[0][2] = n_quad[0].Z();
			xyz[1][0] = n_quad[1].X();
			xyz[1][1] = n_quad[1].Y();
			xyz[1][2] = n_quad[1].Z();
			xyz[2][0] = n_quad[2].X();
			xyz[2][1] = n_quad[2].Y();
			xyz[2][2] = n_quad[2].Z();
			xyz[3][0] = n_quad[3].X();
			xyz[3][1] = n_quad[3].Y();
			xyz[3][2] = n_quad[3].Z();

			gmds::math::Point pt0(vertices[0].x, vertices[0].y);
			gmds::math::Point pt1(vertices[1].x, vertices[1].y);
			gmds::math::Point pt2(vertices[2].x, vertices[2].y);
			gmds::math::Triangle t(pt0, pt1, pt2);


			const double volsurf = t.area();

			r2d_poly poly;
			r2d_rvec2 verts[4];

			verts[0].x = xyz[0][0];
			verts[0].y = xyz[0][1];
			verts[1].x = xyz[1][0];
			verts[1].y = xyz[1][1];
			verts[2].x = xyz[2][0];
			verts[2].y = xyz[2][1];
			verts[3].x = xyz[3][0];
			verts[3].y = xyz[3][1];

			r2d_init_poly(&poly, verts, 4);

			r2d_clip(&poly, planes, 3);
			r2d_int POLY_ORDER = 2;
			r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
			r2d_reduce(&poly, om, POLY_ORDER);

			double vf = AVolFrac->value(tri_id);

			std::cout<<"Id Face : "<<tri_id<<std::endl;
			std::cout<<"valeur om : "<<om<<std::endl;
			std::cout<<"valeur om[0] : "<<om[0]<<std::endl;
			std::cout<<"valeur IoU : "<<volsurf<<std::endl;
			std::cout<<"old value : "<<vf<<std::endl;

			//			std::cout<<"f_id "<<f_id<<" volsurf "<<volsurf<<" vf "<<vf<<" om "<<om[0]/volsurf<<std::endl;
			// add but do not forget to divide by cell area
			AVolFrac->set(tri_id, vf + om[0]/volsurf);
		}
	}
}
