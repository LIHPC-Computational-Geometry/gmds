#include <gtest/gtest.h>
#include <gmds/math/BezierTriangle.h>

using namespace gmds::math;

TEST(BezierTriangleClass, ConstructeurParDefaut)
{
    // Test du constructeur par défaut
    BezierTriangle bezierTriangle;
    ASSERT_NO_THROW(bezierTriangle.operator()(0.0, 0.0));
    ASSERT_NO_THROW(bezierTriangle.normal(0.0, 0.0));
    ASSERT_NO_THROW(bezierTriangle.geomInfo(0.0, 0.0, Point(), Vector3d(), Vector3d(), Vector3d()));
    ASSERT_NO_THROW(bezierTriangle.getDiscretization(10));
}

TEST(BezierTriangleClass, ConstructeurAvecPointsDeControleEtNormales)
{
    // Test du constructeur avec points de contrôle et normales
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 0.0, 0.0);
    Point p3(0.0, 1.0, 0.0);

    Vector3d n1(0.0, 0.0, 1.0);
    Vector3d n2(0.0, 0.0, 1.0);
    Vector3d n3(0.0, 0.0, 1.0);

    BezierTriangle bezierTriangle(p1, p2, p3, n1, n2, n3);

    ASSERT_NO_THROW(bezierTriangle.operator()(0.5, 0.5));
    ASSERT_NO_THROW(bezierTriangle.normal(0.5, 0.5));
    ASSERT_NO_THROW(bezierTriangle.geomInfo(0.5, 0.5, Point(), Vector3d(), Vector3d(), Vector3d()));
    ASSERT_NO_THROW(bezierTriangle.getDiscretization(10));
}

TEST(BezierTriangleClass, EvaluationOperateur)
{
    // Test de l'évaluation de l'opérateur()
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 0.0, 0.0);
    Point p3(0.0, 1.0, 0.0);

    Vector3d n1(0.0, 0.0, 1.0);
    Vector3d n2(0.0, 0.0, 1.0);
    Vector3d n3(0.0, 0.0, 1.0);

    BezierTriangle bezierTriangle(p1, p2, p3, n1, n2, n3);

    // Test de l'évaluation du point à différentes coordonnées (u, v)
    ASSERT_NO_THROW(bezierTriangle.operator()(0.0, 0.0));
    ASSERT_NO_THROW(bezierTriangle.operator()(0.5, 0.5));
    ASSERT_NO_THROW(bezierTriangle.operator()(1.0, 1.0));
}

TEST(BezierTriangleClass, CalculNormal)
{
    // Test du calcul de la normale à différentes coordonnées (u, v)
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 0.0, 0.0);
    Point p3(0.0, 1.0, 0.0);

    Vector3d n1(0.0, 0.0, 1.0);
    Vector3d n2(0.0, 0.0, 1.0);
    Vector3d n3(0.0, 0.0, 1.0);

    BezierTriangle bezierTriangle(p1, p2, p3, n1, n2, n3);

    ASSERT_NO_THROW(bezierTriangle.normal(0.0, 0.0));
    ASSERT_NO_THROW(bezierTriangle.normal(0.5, 0.5));
    ASSERT_NO_THROW(bezierTriangle.normal(1.0, 1.0));
}

TEST(BezierTriangleClass, InfoGeometrie)
{
    // Test de la méthode geomInfo à différentes coordonnées (u, v)
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 0.0, 0.0);
    Point p3(0.0, 1.0, 0.0);

    Vector3d n1(0.0, 0.0, 1.0);
    Vector3d n2(0.0, 0.0, 1.0);
    Vector3d n3(0.0, 0.0, 1.0);

    BezierTriangle bezierTriangle(p1, p2, p3, n1, n2, n3);

    Point pointResult;
    Vector3d normalResult;
    Vector3d duResult;
    Vector3d dvResult;

    ASSERT_NO_THROW(bezierTriangle.geomInfo(0.0, 0.0, pointResult, normalResult, duResult, dvResult));
    ASSERT_NO_THROW(bezierTriangle.geomInfo(0.5, 0.5, pointResult, normalResult, duResult, dvResult));
    ASSERT_NO_THROW(bezierTriangle.geomInfo(1.0, 1.0, pointResult, normalResult, duResult, dvResult));
}

TEST(BezierTriangleClass, Discretisation)
{
    // Test de la discrétisation à différentes résolutions
    Point p1(0.0, 0.0, 0.0);
    Point p2(1.0, 0.0, 0.0);
    Point p3(0.0, 1.0, 0.0);

    Vector3d n1(0.0, 0.0, 1.0);
    Vector3d n2(0.0, 0.0, 1.0);
    Vector3d n3(0.0, 0.0, 1.0);

    BezierTriangle bezierTriangle(p1, p2, p3, n1, n2, n3);

    ASSERT_NO_THROW(bezierTriangle.getDiscretization(10));
    ASSERT_NO_THROW(bezierTriangle.getDiscretization(20));
    ASSERT_NO_THROW(bezierTriangle.getDiscretization(30));
}