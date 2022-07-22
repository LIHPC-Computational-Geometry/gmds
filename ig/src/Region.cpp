/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * Region.cpp
 *
 *  Created on: 20 may 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Region.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <iomanip>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Cell.h>
#include <gmds/ig/RegionContainer.h>
#include <gmds/math/Hexahedron.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Pyramid.h>
#include <gmds/math/Prism3.h>
#include <gmds/math/Triangle.h>
#include <gmds/math/Quadrilateral.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    Region::
    Region()
    : Cell(0,GMDS_TETRA,NullID),m_regions_container(0),m_type_id(NullID)
    {
    }
    /*----------------------------------------------------------------------------*/
    Region::
    Region(Mesh* AMesh, const ECellType AType, const TCellID& AID)
    : Cell(AMesh,AType,AID) //mesh, type and id are filled in
    {
        //============================================
        // we keep a reference on the face container
        if(AMesh!=0){
            m_regions_container = AMesh->m_regions_container;
            m_type_id = m_regions_container->getTypeID(AID);
        }
        else {
            m_regions_container = 0;
            m_type_id = NullID;
        }
        
    }
    /*----------------------------------------------------------------------------*/
    Region::
    Region(const Region& AReg)
    : Cell(AReg.m_owner,AReg.m_type,AReg.m_id)
    {
        if(m_owner!=0)
            m_regions_container = m_owner->m_regions_container;
        else
            m_regions_container = 0;
        
        m_type_id = AReg.m_type_id;
    }
    /*----------------------------------------------------------------------------*/
    Region::~Region()
    {}
    /*----------------------------------------------------------------------------*/
    bool Region::operator==(const Region& ARegion) const
    {
        return (m_owner == ARegion.m_owner && m_id==ARegion.m_id);
    }
    /*----------------------------------------------------------------------------*/
    bool Region::operator!=(const Region& ARegion) const
    {
        return (!(*this == ARegion));
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::nbNodes() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val= 4;
        else if (m_type==GMDS_HEX)
            val= 8;
        else if (m_type==GMDS_PYRAMID)
            val= 5;
        else if (m_type==GMDS_PRISM3)
            val= 6;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2N)[m_type_id];
            val= cells.size();
        }
        return val;
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::nbEdges() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=6;
        else if (m_type==GMDS_HEX)
            val=12;
        else if (m_type==GMDS_PYRAMID)
            val=8;
        else if (m_type==GMDS_PRISM3)
            val=9;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2E)[m_type_id];
            val=cells.size();
        }
        return val;
        
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::nbFaces() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=4;
        else if (m_type==GMDS_HEX)
            val=6;
        else if (m_type==GMDS_PYRAMID || m_type==GMDS_PRISM3)
            val=5;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2F)[m_type_id];
            val=cells.size();
        }
        return val;
        
    }
    /*----------------------------------------------------------------------------*/
    TInt Region::nbRegions() const
    {
        TInt val=0;
        if (m_type==GMDS_TETRA)
            val=4;
        else if (m_type==GMDS_PYRAMID || m_type==GMDS_PRISM3)
            val=5;
        else if (m_type==GMDS_HEX)
            val=6;
        else if (m_type==GMDS_POLYHEDRA){
            TabCellID<size_undef> cells = (*m_regions_container->m_P2R)[m_type_id];
            val=cells.size();
        }
        return val;
    }
    /*----------------------------------------------------------------------------*/
    math::Point
    Region::center() const
    {
        TCoord p_coords[3] = {0.0,0.0,0.0};
        
        std::vector<Node> nodes = this->get<Node>();
        unsigned  nb_nodes = nodes.size();
        
        for(unsigned int i=0; i<nb_nodes; i++)
        {
            Node n = nodes[i];
            p_coords[0] += n.X();
            p_coords[1] += n.Y();
            p_coords[2] += n.Z();
        }
        
        p_coords[0] = p_coords[0] / nodes.size();
        p_coords[1] = p_coords[1] / nodes.size();
        p_coords[2] = p_coords[2] / nodes.size();
        
        math::Point p(p_coords[0],p_coords[1],p_coords[2]);
        
        return p;
    }
    
    /*----------------------------------------------------------------------------*/
    TCoord Region::volume() const
    {
        TCoord vol=0.0;
        std::vector<Node> n = this->get<Node>();
        if(this->type()==GMDS_TETRA){
            vol = math::Tetrahedron(n[0].point(),
                                    n[1].point(),
                                    n[2].point(),
                                    n[3].point()).getVolume();
        }
        else if(this->type()==GMDS_HEX){
            vol = math::Hexahedron(n[0].point(),
                                   n[1].point(),
                                   n[2].point(),
                                   n[3].point(),
                                   n[4].point(),
                                   n[5].point(),
                                   n[6].point(),
                                   n[7].point()).getVolume();
	}
        else if(this->type()==GMDS_PYRAMID){
            vol = math::Pyramid(n[0].point(),
                                n[1].point(),
                                n[2].point(),
                                n[3].point(),
                                n[4].point()).getVolume();
        }
        else if(this->type()==GMDS_PRISM3){
            vol = math::Prism3(n[0].point(),
                               n[1].point(),
                               n[2].point(),
                               n[3].point(),
                               n[4].point(),
                               n[5].point()).getVolume();
        }
	else {
	  throw GMDSException("Region::volume can not be computed for this value type.");
	}
        
        return vol;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<math::Point>
    Region::computeNGLLPoints(int ADegree) const
    {
        //	throw GMDSException("Region::computeNGLLPoints not implemented yet.");
        
        std::vector<math::Point> points;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
            {
                points.resize(27);
                
                // we value the 8 corners of the hex
                gmds::math::Point p[8];
                for(unsigned int iNode=0; iNode<nodes.size(); iNode++)
                {
                    p[iNode] = nodes[iNode].point();
                }
                // Now we build an inner grid of  27 points
                
                // we build intermediate nodes to create our 27 inner points
                // middle of edge
                gmds::math::Point p01, p12, p23, p03,
                p45, p56, p67, p47, p04, p15, p26, p37;
                
                p01 = (p[0]+p[1])*0.5;
                p12 = (p[1]+p[2])*0.5;
                p23 = (p[2]+p[3])*0.5;
                p03 = (p[0]+p[3])*0.5;
                p45 = (p[4]+p[5])*0.5;
                p56 = (p[5]+p[6])*0.5;
                p67 = (p[6]+p[7])*0.5;
                p47 = (p[4]+p[7])*0.5;
                p04 = (p[0]+p[4])*0.5;
                p15 = (p[1]+p[5])*0.5;
                p26 = (p[2]+p[6])*0.5;
                p37 = (p[3]+p[7])*0.5;																	  //face center
                gmds::math::Point p0123, p4567, p0154, p1265, p2376, p0374;
                p0123 = (p[0]+p[1]+p[2]+p[3])*0.25;
                p4567 = (p[4]+p[5]+p[6]+p[7])*0.25;
                p0154 = (p[0]+p[1]+p[5]+p[4])*0.25;
                p1265 = (p[1]+p[2]+p[6]+p[5])*0.25;
                p2376 = (p[2]+p[3]+p[7]+p[6])*0.25;
                p0374 = (p[0]+p[3]+p[7]+p[4])*0.25;
                //hex center
                gmds::math::Point bary = (p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7])*0.125;
                points[0] = bary;
                points[1] = (bary+p[0])*0.5;
                points[2] = (bary+p[1])*0.5;
                points[3] = (bary+p[2])*0.5;
                points[4] = (bary+p[3])*0.5;
                points[5] = (bary+p[4])*0.5;
                points[6] = (bary+p[5])*0.5;
                points[7] = (bary+p[6])*0.5;
                points[8] = (bary+p[7])*0.5;
                points[9] = (bary+p01)*0.5;
                points[10]= (bary+p12)*0.5;
                points[11]= (bary+p23)*0.5;
                points[12]= (bary+p03)*0.5;
                points[13]= (bary+p45)*0.5;
                points[14]= (bary+p56)*0.5;
                points[15]= (bary+p67)*0.5;
                points[16]= (bary+p47)*0.5;
                points[17]= (bary+p04)*0.5;
                points[18]= (bary+p15)*0.5;
                points[19]= (bary+p26)*0.5;
                points[20]= (bary+p37)*0.5;
                points[21]= (bary+p0123)*0.5;
                points[22]= (bary+p4567)*0.5;
                points[23]= (bary+p0154)*0.5;
                points[24]= (bary+p1265)*0.5;
                points[25]= (bary+p2376)*0.5;
                points[26]= (bary+p0374)*0.5;
            }
                break;
            case GMDS_TETRA:
                for(unsigned int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].point());
                }
                break;
            case GMDS_PYRAMID:
                for(unsigned int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].point());
                }
                break;
            case GMDS_PRISM3:
                for(unsigned int iNode=0; iNode<nodes.size(); iNode++)
                {
                    points.push_back(nodes[iNode].point());
                }
                break;
            default:
                throw GMDSException("Region::computeNGLLPoints not implemented yet for this value type.");
                break;
        }
        
        return points;
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeScaledJacobian() const
    {
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                                     nodes[4].point(), nodes[5].point(), nodes[6].point(), nodes[7].point());
                return hex.computeScaledJacobian();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
                return tet.computeScaledJacobian();
            }
                break;
            case GMDS_PYRAMID:
	    {
		math::Pyramid pyr(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                          nodes[4].point());
		return pyr.computeScaledJacobian();
	    }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].point(), nodes[1].point(), nodes[2].point(),
                                   nodes[3].point(), nodes[4].point(), nodes[5].point());
                return prism.computeScaledJacobian();
            }
                break;
            default:
                throw GMDSException("Region::computeScaledJacobian not implemented yet for this value type.");
                break;
        }
return 0;
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeNormalizedScaledJacobian() const
    {
        std::vector<Node> nodes = this->get<Node>();

        switch(this->type()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                                     nodes[4].point(), nodes[5].point(), nodes[6].point(), nodes[7].point());
                return hex.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
                return tet.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_PYRAMID:
            {
                math::Pyramid pyr(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                                  nodes[4].point());
                return pyr.computeNormalizedScaledJacobian();
            }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].point(), nodes[1].point(), nodes[2].point(),
                                   nodes[3].point(), nodes[4].point(), nodes[5].point());
                return prism.computeNormalizedScaledJacobian();
            }
                break;
            default:
                throw GMDSException("Region::computeNormalizedScaledJacobian not implemented yet for this value type.");
                break;
        }
        return 0;
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::computeMeanRatio() const
    {
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
            {
                math::Hexahedron hex(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                                     nodes[4].point(), nodes[5].point(), nodes[6].point(), nodes[7].point());
                return hex.computeMeanRatio();
            }
                break;
            case GMDS_TETRA:
            {
                math::Tetrahedron tet(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point());
                return tet.computeMeanRatio();
            }
                break;
            case GMDS_PYRAMID:
	    {
		math::Pyramid pyr(nodes[0].point(), nodes[1].point(), nodes[2].point(), nodes[3].point(),
                          nodes[4].point());
		return pyr.computeMeanRatio();
	    }
                break;
            case GMDS_PRISM3:
            {
                math::Prism3 prism(nodes[0].point(), nodes[1].point(), nodes[2].point(),
                                   nodes[3].point(), nodes[4].point(), nodes[5].point());
                return prism.computeMeanRatio();
            }
                break;
            default:
                throw GMDSException("Region::computeMeanRatio not implemented yet for this value type.");
                break;
        }
        return 0;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<std::vector<Node> >
    Region::getOrderedNodesFaces() const
    {
        std::vector<std::vector<Node> > orderedNodesFaces;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
                orderedNodesFaces.resize(6);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[7];
                orderedNodesFaces[1][1] = nodes[4];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[6];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[5];
                orderedNodesFaces[2][1] = nodes[4];
                orderedNodesFaces[2][2] = nodes[0];
                orderedNodesFaces[2][3] = nodes[1];
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[7];
                orderedNodesFaces[3][1] = nodes[6];
                orderedNodesFaces[3][2] = nodes[2];
                orderedNodesFaces[3][3] = nodes[3];
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[4];
                orderedNodesFaces[4][1] = nodes[7];
                orderedNodesFaces[4][2] = nodes[3];
                orderedNodesFaces[4][3] = nodes[0];
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[6];
                orderedNodesFaces[5][1] = nodes[5];
                orderedNodesFaces[5][2] = nodes[1];
                orderedNodesFaces[5][3] = nodes[2];
                break;
            case GMDS_TETRA:
                orderedNodesFaces.resize(4);
                orderedNodesFaces[0].resize(3);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[2];
                orderedNodesFaces[0][2] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[3];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[0];
                orderedNodesFaces[3][2] = nodes[3];
                break;
            case GMDS_PYRAMID:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[4];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[4];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[3];
                orderedNodesFaces[3][2] = nodes[4];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[0];
                orderedNodesFaces[4][2] = nodes[4];
                break;
            case GMDS_PRISM3:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[4];
                orderedNodesFaces[0][3] = nodes[3];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[1];
                orderedNodesFaces[1][1] = nodes[2];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[4];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[2];
                orderedNodesFaces[2][1] = nodes[0];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[2][3] = nodes[5];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[0];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[1];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[4];
                orderedNodesFaces[4][2] = nodes[5];
                break;
            default:
                throw GMDSException("Region::getOrderedNodesFaces not implemented for this region type");
                break;
        }
        
        return orderedNodesFaces;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<std::vector<TCellID> >
    Region::getOrderedNodesFacesIDs() const
    {
        std::vector<std::vector<TCellID> > orderedNodesFaces;
        
        std::vector<TCellID> nodes = this->getIDs<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
                orderedNodesFaces.resize(6);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[7];
                orderedNodesFaces[1][1] = nodes[4];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[6];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[5];
                orderedNodesFaces[2][1] = nodes[4];
                orderedNodesFaces[2][2] = nodes[0];
                orderedNodesFaces[2][3] = nodes[1];
                orderedNodesFaces[3].resize(4);
                orderedNodesFaces[3][0] = nodes[7];
                orderedNodesFaces[3][1] = nodes[6];
                orderedNodesFaces[3][2] = nodes[2];
                orderedNodesFaces[3][3] = nodes[3];
                orderedNodesFaces[4].resize(4);
                orderedNodesFaces[4][0] = nodes[4];
                orderedNodesFaces[4][1] = nodes[7];
                orderedNodesFaces[4][2] = nodes[3];
                orderedNodesFaces[4][3] = nodes[0];
                orderedNodesFaces[5].resize(4);
                orderedNodesFaces[5][0] = nodes[6];
                orderedNodesFaces[5][1] = nodes[5];
                orderedNodesFaces[5][2] = nodes[1];
                orderedNodesFaces[5][3] = nodes[2];
                break;
            case GMDS_TETRA:
                orderedNodesFaces.resize(4);
                orderedNodesFaces[0].resize(3);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[2];
                orderedNodesFaces[0][2] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[3];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[0];
                orderedNodesFaces[3][2] = nodes[3];
                break;
            case GMDS_PYRAMID:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[3];
                orderedNodesFaces[0][2] = nodes[2];
                orderedNodesFaces[0][3] = nodes[1];
                orderedNodesFaces[1].resize(3);
                orderedNodesFaces[1][0] = nodes[0];
                orderedNodesFaces[1][1] = nodes[1];
                orderedNodesFaces[1][2] = nodes[4];
                orderedNodesFaces[2].resize(3);
                orderedNodesFaces[2][0] = nodes[1];
                orderedNodesFaces[2][1] = nodes[2];
                orderedNodesFaces[2][2] = nodes[4];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[2];
                orderedNodesFaces[3][1] = nodes[3];
                orderedNodesFaces[3][2] = nodes[4];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[0];
                orderedNodesFaces[4][2] = nodes[4];
                break;
            case GMDS_PRISM3:
                orderedNodesFaces.resize(5);
                orderedNodesFaces[0].resize(4);
                orderedNodesFaces[0][0] = nodes[0];
                orderedNodesFaces[0][1] = nodes[1];
                orderedNodesFaces[0][2] = nodes[4];
                orderedNodesFaces[0][3] = nodes[3];
                orderedNodesFaces[1].resize(4);
                orderedNodesFaces[1][0] = nodes[1];
                orderedNodesFaces[1][1] = nodes[2];
                orderedNodesFaces[1][2] = nodes[5];
                orderedNodesFaces[1][3] = nodes[4];
                orderedNodesFaces[2].resize(4);
                orderedNodesFaces[2][0] = nodes[2];
                orderedNodesFaces[2][1] = nodes[0];
                orderedNodesFaces[2][2] = nodes[3];
                orderedNodesFaces[2][3] = nodes[5];
                orderedNodesFaces[3].resize(3);
                orderedNodesFaces[3][0] = nodes[0];
                orderedNodesFaces[3][1] = nodes[2];
                orderedNodesFaces[3][2] = nodes[1];
                orderedNodesFaces[4].resize(3);
                orderedNodesFaces[4][0] = nodes[3];
                orderedNodesFaces[4][1] = nodes[4];
                orderedNodesFaces[4][2] = nodes[5];
                break;
            default:
                throw GMDSException("Region::getOrderedNodesFacesIDs not implemented for this region type");
                break;
        }
        
        return orderedNodesFaces;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<std::vector<TCellID> >
    Region::getNodesEdgesIDs() const
    {
        std::vector<std::vector<TCellID> > nodesEdges;

        std::vector<TCellID> nodes = this->getAllIDs<Node>();

        switch(this->type()) {
            case GMDS_HEX:
                nodesEdges.resize(12);
                nodesEdges[ 0].resize(2);
                nodesEdges[ 0][0] = nodes[0];
                nodesEdges[ 0][1] = nodes[1];
                nodesEdges[ 1].resize(2);
                nodesEdges[ 1][0] = nodes[1];
                nodesEdges[ 1][1] = nodes[2];
                nodesEdges[ 2].resize(2);
                nodesEdges[ 2][0] = nodes[2];
                nodesEdges[ 2][1] = nodes[3];
                nodesEdges[ 3].resize(2);
                nodesEdges[ 3][0] = nodes[3];
                nodesEdges[ 3][1] = nodes[0];
                nodesEdges[ 4].resize(2);
                nodesEdges[ 4][0] = nodes[4];
                nodesEdges[ 4][1] = nodes[5];
                nodesEdges[ 5].resize(2);
                nodesEdges[ 5][0] = nodes[5];
                nodesEdges[ 5][1] = nodes[6];
                nodesEdges[ 6].resize(2);
                nodesEdges[ 6][0] = nodes[6];
                nodesEdges[ 6][1] = nodes[7];
                nodesEdges[ 7].resize(2);
                nodesEdges[ 7][0] = nodes[7];
                nodesEdges[ 7][1] = nodes[0];
                nodesEdges[ 8].resize(2);
                nodesEdges[ 8][0] = nodes[0];
                nodesEdges[ 8][1] = nodes[4];
                nodesEdges[ 9].resize(2);
                nodesEdges[ 9][0] = nodes[1];
                nodesEdges[ 9][1] = nodes[5];
                nodesEdges[10].resize(2);
                nodesEdges[10][0] = nodes[2];
                nodesEdges[10][1] = nodes[6];
                nodesEdges[11].resize(2);
                nodesEdges[11][0] = nodes[3];
                nodesEdges[11][1] = nodes[7];
                break;
            case GMDS_TETRA:
                nodesEdges.resize(6);
                nodesEdges[0].resize(2);
                nodesEdges[0][0] = nodes[0];
                nodesEdges[0][1] = nodes[1];
                nodesEdges[1].resize(2);
                nodesEdges[1][0] = nodes[1];
                nodesEdges[1][1] = nodes[2];
                nodesEdges[2].resize(2);
                nodesEdges[2][0] = nodes[2];
                nodesEdges[2][1] = nodes[0];
                nodesEdges[3].resize(2);
                nodesEdges[3][0] = nodes[0];
                nodesEdges[3][1] = nodes[3];
                nodesEdges[4].resize(2);
                nodesEdges[4][0] = nodes[1];
                nodesEdges[4][1] = nodes[3];
                nodesEdges[5].resize(2);
                nodesEdges[5][0] = nodes[2];
                nodesEdges[5][1] = nodes[3];
                break;
            case GMDS_PYRAMID:
                nodesEdges.resize(8);
                nodesEdges[0].resize(2);
                nodesEdges[0][0] = nodes[0];
                nodesEdges[0][1] = nodes[1];
                nodesEdges[1].resize(2);
                nodesEdges[1][0] = nodes[1];
                nodesEdges[1][1] = nodes[2];
                nodesEdges[2].resize(2);
                nodesEdges[2][0] = nodes[2];
                nodesEdges[2][1] = nodes[3];
                nodesEdges[3].resize(2);
                nodesEdges[3][0] = nodes[3];
                nodesEdges[3][1] = nodes[0];
                nodesEdges[4].resize(2);
                nodesEdges[4][0] = nodes[0];
                nodesEdges[4][1] = nodes[4];
                nodesEdges[5].resize(2);
                nodesEdges[5][0] = nodes[1];
                nodesEdges[5][1] = nodes[4];
                nodesEdges[6].resize(2);
                nodesEdges[6][0] = nodes[2];
                nodesEdges[6][1] = nodes[4];
                nodesEdges[7].resize(2);
                nodesEdges[7][0] = nodes[3];
                nodesEdges[7][1] = nodes[4];
                break;
            case GMDS_PRISM3:
                nodesEdges.resize(9);
                nodesEdges[0].resize(2);
                nodesEdges[0][0] = nodes[0];
                nodesEdges[0][1] = nodes[1];
                nodesEdges[1].resize(2);
                nodesEdges[1][0] = nodes[1];
                nodesEdges[1][1] = nodes[2];
                nodesEdges[2].resize(2);
                nodesEdges[2][0] = nodes[2];
                nodesEdges[2][1] = nodes[0];
                nodesEdges[3].resize(2);
                nodesEdges[3][0] = nodes[3];
                nodesEdges[3][1] = nodes[4];
                nodesEdges[4].resize(2);
                nodesEdges[4][0] = nodes[4];
                nodesEdges[4][1] = nodes[5];
                nodesEdges[5].resize(2);
                nodesEdges[5][0] = nodes[5];
                nodesEdges[5][1] = nodes[3];
                nodesEdges[6].resize(2);
                nodesEdges[6][0] = nodes[0];
                nodesEdges[6][1] = nodes[3];
                nodesEdges[7].resize(2);
                nodesEdges[7][0] = nodes[1];
                nodesEdges[7][1] = nodes[4];
                nodesEdges[8].resize(2);
                nodesEdges[8][0] = nodes[2];
                nodesEdges[8][1] = nodes[5];
                break;
            default:
                throw GMDSException("Region::getNodesEdgesIDs not implemented for this region type");
                break;
        }

        return nodesEdges;
    }
    /*----------------------------------------------------------------------------*/
    bool
    Region::isFaceOrientedOutward(std::vector<Node> ANodes) const
    {
      std::vector<TCellID> ids(ANodes.size());
      for(size_t i=0; i<ANodes.size(); i++) {
	ids[i] = ANodes[i].id();
      }

      try{
	return this->isFaceOrientedOutward(ids);
      } catch(GMDSException& e) {
	throw e;
      }
      
    }
    /*----------------------------------------------------------------------------*/
    bool
    Region::isFaceOrientedOutward(std::vector<TCellID> AIDs) const
    {
        std::vector<std::vector<TCellID> > orderedNodesFaces = this->getOrderedNodesFacesIDs();
        
        // first find the face
        for(unsigned int iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

            if(AIDs.size() == orderedNodesFaces[iFace].size()) {
                
                unsigned int nbNodesMatched = 0;
                
                for(unsigned int iNode1=0; iNode1<orderedNodesFaces[iFace].size(); iNode1++) {
                    for(unsigned int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                        if(AIDs[iNode2] == orderedNodesFaces[iFace][iNode1]) {
                            nbNodesMatched++;
                        }
                    }
                }
                if(nbNodesMatched == AIDs.size()) {
                    // face is found, now check the orientation
                    if(orderedNodesFaces[iFace].size() < 3) {
                        throw GMDSException("Region::isFaceOrientedOutward face with less than 3 nodes.");
                    }
                    TCellID firstNode = orderedNodesFaces[iFace][0];
                    TCellID secondNode = orderedNodesFaces[iFace][1];
                    
                    for(unsigned int iNode2=0; iNode2<AIDs.size(); iNode2++) {
                        if(AIDs[iNode2] == firstNode) {
                            if(AIDs[(iNode2+1)%AIDs.size()] == secondNode) {
                                return true;
                            } else  if (AIDs[(AIDs.size()+iNode2-1)%AIDs.size()] == secondNode) {
                                return false;
                            } else {
                                throw GMDSException("Region::isFaceOrientedOutward jumbled face.");
                            }
                        }
                    }
                }
            }
        }
        
        throw GMDSException("Region::isFaceOrientedOutward face not found.");
    }
   /*----------------------------------------------------------------------------*/
    std::vector<VirtualFace>
    Region::getFakeFaces() const
    {
	std::vector<std::vector<Node> > orderedNodesFaces = this->getOrderedNodesFaces();
	std::vector<VirtualFace> fakeFaces;

	for(size_t iFace=0; iFace<orderedNodesFaces.size(); iFace++) {

		std::vector<gmds::TCellID> ids(orderedNodesFaces[iFace].size());
		for(size_t iNode=0; iNode<ids.size(); iNode++) {
			ids[iNode] = orderedNodesFaces[iFace][iNode].id();
		}

		fakeFaces.push_back(VirtualFace(ids));
	}

        return fakeFaces;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<VirtualEdge>
    Region::getFakeEdges() const
    {
        std::vector<VirtualEdge> fakeEdges;
        
        std::vector<Node> nodes = this->get<Node>();
        
        switch(this->type()) {
            case GMDS_HEX:
                fakeEdges.resize(12);
                fakeEdges[ 0] = VirtualEdge(nodes[0].id(), nodes[1].id());
                fakeEdges[ 1] = VirtualEdge(nodes[1].id(), nodes[2].id());
                fakeEdges[ 2] = VirtualEdge(nodes[2].id(), nodes[3].id());
                fakeEdges[ 3] = VirtualEdge(nodes[3].id(), nodes[0].id());
                fakeEdges[ 4] = VirtualEdge(nodes[4].id(), nodes[5].id());
                fakeEdges[ 5] = VirtualEdge(nodes[5].id(), nodes[6].id());
                fakeEdges[ 6] = VirtualEdge(nodes[6].id(), nodes[7].id());
                fakeEdges[ 7] = VirtualEdge(nodes[7].id(), nodes[4].id());
                fakeEdges[ 8] = VirtualEdge(nodes[0].id(), nodes[4].id());
                fakeEdges[ 9] = VirtualEdge(nodes[1].id(), nodes[5].id());
                fakeEdges[10] = VirtualEdge(nodes[2].id(), nodes[6].id());
                fakeEdges[11] = VirtualEdge(nodes[3].id(), nodes[7].id());
                break;
            case GMDS_TETRA:
                fakeEdges.resize(6);
		fakeEdges[ 0] = VirtualEdge(nodes[0].id(), nodes[1].id());
		fakeEdges[ 1] = VirtualEdge(nodes[1].id(), nodes[2].id());
		fakeEdges[ 2] = VirtualEdge(nodes[2].id(), nodes[0].id());
		fakeEdges[ 3] = VirtualEdge(nodes[0].id(), nodes[3].id());
		fakeEdges[ 4] = VirtualEdge(nodes[1].id(), nodes[3].id());
		fakeEdges[ 5] = VirtualEdge(nodes[2].id(), nodes[3].id());
                break;
            case GMDS_PYRAMID:
                fakeEdges.resize(8);
		fakeEdges[ 0] = VirtualEdge(nodes[0].id(), nodes[1].id());
                fakeEdges[ 1] = VirtualEdge(nodes[1].id(), nodes[2].id());
                fakeEdges[ 2] = VirtualEdge(nodes[2].id(), nodes[3].id());
                fakeEdges[ 3] = VirtualEdge(nodes[3].id(), nodes[0].id());
                fakeEdges[ 4] = VirtualEdge(nodes[0].id(), nodes[4].id());
                fakeEdges[ 5] = VirtualEdge(nodes[1].id(), nodes[4].id());
                fakeEdges[ 6] = VirtualEdge(nodes[2].id(), nodes[4].id());
                fakeEdges[ 7] = VirtualEdge(nodes[3].id(), nodes[4].id());
                break;
            case GMDS_PRISM3:
                fakeEdges.resize(9);
		fakeEdges[ 0] = VirtualEdge(nodes[0].id(), nodes[1].id());
                fakeEdges[ 1] = VirtualEdge(nodes[1].id(), nodes[2].id());
                fakeEdges[ 2] = VirtualEdge(nodes[2].id(), nodes[0].id());
                fakeEdges[ 3] = VirtualEdge(nodes[3].id(), nodes[4].id());
                fakeEdges[ 4] = VirtualEdge(nodes[4].id(), nodes[5].id());
                fakeEdges[ 5] = VirtualEdge(nodes[5].id(), nodes[3].id());
                fakeEdges[ 6] = VirtualEdge(nodes[0].id(), nodes[3].id());
                fakeEdges[ 7] = VirtualEdge(nodes[1].id(), nodes[4].id());
                fakeEdges[ 8] = VirtualEdge(nodes[2].id(), nodes[5].id());
                break;
            default:
                throw GMDSException("Region::getFakeEdges not implemented for this region type");
                break;
        }
        
        return fakeEdges;
    }
    /*----------------------------------------------------------------------------*/
    double
    Region::findTheNorm(double ATargetVolume, std::set<VirtualFace> freeFaces) const
    {
        if(ATargetVolume <= 0.) {
            throw GMDSException("Region::findTheNorm cannot be called with a negative or null target volume");
        }

        if(freeFaces.empty()) {
            throw GMDSException("Region::findTheNorm cannot be called without any free face");
        }

        // TODO: choose good constants
        const double stoppingCriteria = 1.e-8;
        const int numIter = 100;

        double current_volume = this->volume();

        if(std::abs(current_volume - ATargetVolume)/ATargetVolume < stoppingCriteria) {
            return 0.;
        }

        // identify the free nodes and their associated faces
        std::vector<Node> nodes = this->get<Node>();
        std::map<TCellID, int> nodeID2localID;
        std::vector<gmds::math::Point> points(8);
        for(unsigned int i=0; i<nodes.size(); i++) {
            nodeID2localID[nodes[i].id()] = i;
            points[i] = nodes[i].point();
        }

        std::vector<bool> isNodeFree(nodes.size(),false);
        std::vector<gmds::math::Vector3d> displacement(nodes.size(), gmds::math::Vector3d({0.,0.,0.}));
        std::vector<VirtualFace> ffs = this->getFakeFaces();
        for(auto ff: ffs) {
            if(freeFaces.find(ff) != freeFaces.end()) {
                std::vector<TCellID> faceNodeIDs = ff.node_ids();
                std::vector<gmds::math::Point> facePoints;
                for(auto n: faceNodeIDs) {
                    facePoints.push_back(points[nodeID2localID[n]]);
                }
                switch (facePoints.size()) {
                    case 3 : {
                        gmds::math::Triangle tri(facePoints[0], facePoints[1], facePoints[2]);
                        gmds::math::Vector3d displ(tri.getNormal());
                        displ.normalize();
                        if(!this->isFaceOrientedOutward(faceNodeIDs)) {
                            displ = (-1.) * displ;
                        }
                        for(auto n: faceNodeIDs) {
                            displacement[nodeID2localID[n]] = displacement[nodeID2localID[n]] + displ;
                        }
                    }
                        break;
                    case 4 : {
                        gmds::math::Quadrilateral quad(facePoints[0], facePoints[1], facePoints[2], facePoints[3]);
                        gmds::math::Vector3d displ(quad.getNormal());
                        displ.normalize();
                        if(!this->isFaceOrientedOutward(faceNodeIDs)) {
                            displ = (-1.) * displ;
                        }
                        for(auto n: faceNodeIDs) {
                            displacement[nodeID2localID[n]] = displacement[nodeID2localID[n]] + displ;
                        }
                    }
                        break;
                    default :
                        break;
                }

            }
        }

        double adjustlength_under = (current_volume < ATargetVolume) ? 0. : -2.*std::cbrt(ATargetVolume);
        double adjustlength_over  = (current_volume < ATargetVolume) ? 2.*std::cbrt(ATargetVolume) : 0.;

        // augment adjustlength_under till the adjusted volume is greater than ATargetVolume
        {
            double under_volume = HUGE_VALF;

            while (under_volume > ATargetVolume) {

                gmds::math::Hexahedron hex(
                        points[0] + adjustlength_under*displacement[0],
                        points[1] + adjustlength_under*displacement[1],
                        points[2] + adjustlength_under*displacement[2],
                        points[3] + adjustlength_under*displacement[3],
                        points[4] + adjustlength_under*displacement[4],
                        points[5] + adjustlength_under*displacement[5],
                        points[6] + adjustlength_under*displacement[6],
                        points[7] + adjustlength_under*displacement[7]
                );

                under_volume = hex.getVolume();

                if(under_volume > ATargetVolume) {
                    break;
                } else {
                    adjustlength_under *= 2.;
                }

            }
        }

        // augment adjustlength_over till the adjusted volume is greater than ATargetVolume
        {
            //std::cout<<"finding adjustlength_over "<<std::endl;

            double over_volume = - HUGE_VALF;

            while (over_volume < ATargetVolume) {

                gmds::math::Hexahedron hex(
                        points[0] + adjustlength_over*displacement[0],
                        points[1] + adjustlength_over*displacement[1],
                        points[2] + adjustlength_over*displacement[2],
                        points[3] + adjustlength_over*displacement[3],
                        points[4] + adjustlength_over*displacement[4],
                        points[5] + adjustlength_over*displacement[5],
                        points[6] + adjustlength_over*displacement[6],
                        points[7] + adjustlength_over*displacement[7]
                );

                //std::cout<<over_volume<<" "<<ATargetVolume<<" "<<adjustlength_over<<std::endl;

                over_volume = hex.getVolume();

                if(over_volume > ATargetVolume) {
                    break;
                } else {
                    adjustlength_over *= 2.;
                }

            }
        }

        double next_volume = HUGE_VALF;
        int iter = 0;

        while (iter < numIter) {

            double adjustlength_current = (adjustlength_over + adjustlength_under) / 2.;
            gmds::math::Hexahedron hex(
                points[0] + adjustlength_current*displacement[0],
                points[1] + adjustlength_current*displacement[1],
                points[2] + adjustlength_current*displacement[2],
                points[3] + adjustlength_current*displacement[3],
                points[4] + adjustlength_current*displacement[4],
                points[5] + adjustlength_current*displacement[5],
                points[6] + adjustlength_current*displacement[6],
                points[7] + adjustlength_current*displacement[7]
            );
            //std::cout<<hex<<std::endl;

            next_volume = hex.getVolume();

            //std::cout<<"iter "<<iter<<" "<<adjustlength_under<<" "<<adjustlength_over<<" "<<next_volume<<std::endl;

            if(std::abs(next_volume - ATargetVolume)/ATargetVolume < stoppingCriteria) {
                return adjustlength_current;
            }

            if(std::abs(next_volume - current_volume)/ ATargetVolume < stoppingCriteria) {
                return adjustlength_current;
            }

            current_volume = next_volume;
            if(next_volume < ATargetVolume){
                adjustlength_under = adjustlength_current;
            } else {
                adjustlength_over = adjustlength_current;
            }
            iter++;
        }

        /*
        std::cerr<<std::setprecision(16)<<std::endl;
        std::cerr<<this->volume()<<" "<<ATargetVolume<<" "<<current_volume<<" "<<next_volume<<" "<<adjustlength_under<<" "<<adjustlength_over<<std::endl;
        std::cerr<<"gmds::Node n0 = mesh.newNode"<<points[0]<<std::endl;
        std::cerr<<"gmds::Node n1 = mesh.newNode"<<points[1]<<std::endl;
        std::cerr<<"gmds::Node n2 = mesh.newNode"<<points[2]<<std::endl;
        std::cerr<<"gmds::Node n3 = mesh.newNode"<<points[3]<<std::endl;
        std::cerr<<"gmds::Node n4 = mesh.newNode"<<points[4]<<std::endl;
        std::cerr<<"gmds::Node n5 = mesh.newNode"<<points[5]<<std::endl;
        std::cerr<<"gmds::Node n6 = mesh.newNode"<<points[6]<<std::endl;
        std::cerr<<"gmds::Node n7 = mesh.newNode"<<points[7]<<std::endl;
        for(auto v: displacement) {
            std::cerr<<v<<std::endl;
        }
        std::cerr<<std::endl;
         */
        throw GMDSException("Region::findTheNorm could not find the adjustment length.");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Node>& ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Edge>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Face>& ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGet(std::vector<Region>& ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].values(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].values(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_regions_container->buildRegion(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetNodeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetEdgeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetFaceIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetRegionIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].values(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].values(ACells);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Node>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_nodes_container->buildNode(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Edge>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_edges_container->buildEdge(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Face>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_owner->m_faces_container->buildFace(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAll(std::vector<Region>& ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        std::vector<TCellID> cellIDs;
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].allValues(cellIDs);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].allValues(cellIDs);
        }
        else
            throw GMDSException("Not yet implemented");
        
        ACells.resize(cellIDs.size());
        for(unsigned int i=0;i<cellIDs.size();i++){
            ACells[i] =m_regions_container->buildRegion(cellIDs[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllNodeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");

        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllEdgeIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllFaceIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateGetAllRegionIDs(std::vector<TCellID>&   ACells) const
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
            
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].allValues(ACells);
        }   
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].allValues(ACells);
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2R)[m_type_id].allValues(ACells);
        }
        else
            throw GMDSException("Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetNodeIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");

        unsigned  nb_cells = nbNodes();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_H2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2N)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2N)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetEdgeIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");


        unsigned  nb_cells = nbEdges();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_H2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_PY2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            (*m_regions_container->m_PR2E)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2E)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetFaceIDs(const std::vector<TCellID>&   ACells)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");

        unsigned int nb_cells = nbFaces();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size()){
                for(auto c_id:ACells){
                    Face f = m_owner->get<Face>(c_id);
                    std::cout<<"Face: "<<c_id<<std::endl;
                    std::vector<Node> fn =f.get<Node>();
                    for(auto ni:fn){
                        std::cout << "\t Node " << ni.id() << ": " << ni.point() << "\n";
                    }
                }
                throw GMDSException("Invalid number of adj. entities");
            }
            (*m_regions_container->m_H2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2F)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            (*m_regions_container->m_P2F)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateSetRegionIDs(const std::vector<TCellID>& ACells)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");


        unsigned int nb_cells = nbRegions();
        if(m_type==GMDS_TETRA){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_T2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_HEX){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_H2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PYRAMID){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PY2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_PRISM3){
            if(nb_cells!=ACells.size())
                throw GMDSException("Invalid number of adj. entities");
            
            (*m_regions_container->m_PR2R)[m_type_id]=ACells;
        }
        else if (m_type==GMDS_POLYHEDRA){
            
            (*m_regions_container->m_P2R)[m_type_id]=ACells;
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2E adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionAdd(TCellID AID)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].add(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].add(AID);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].add(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2N adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2F adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionRemove(TCellID AID)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].del(AID);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].del(AID);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].del(AID);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateNodeReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2N))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2N)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2N)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2N)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateEdgeReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2E))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2E)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2E)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2E)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateFaceReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2F))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2F)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2F)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2F)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    void Region::delegateRegionReplace(TCellID AID1, TCellID AID2)
    {
        if(!m_owner->m_model.has(R2R))
            throw GMDSException("R2R adjacency is not supported by the mesh model");
        
        if(m_type==GMDS_TETRA){
            (*m_regions_container->m_T2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_HEX){
            (*m_regions_container->m_H2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PYRAMID){
            (*m_regions_container->m_PY2R)[m_type_id].replace(AID1, AID2);
        }
        else if (m_type==GMDS_PRISM3){
            (*m_regions_container->m_PR2R)[m_type_id].replace(AID1, AID2);
        }
        else {
            (*m_regions_container->m_P2R)[m_type_id].replace(AID1, AID2);
        }
    }
    /*----------------------------------------------------------------------------*/
    std::ostream & operator << (std::ostream & AStream, const Region & AR)
    {
        AStream<<"Region "<< AR.id();;
        return AStream;
    }
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
