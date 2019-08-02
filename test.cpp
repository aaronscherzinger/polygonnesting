#include <algorithm>
#include <iostream>

#define POLYGONNESTING_PRINT_DEBUG
#include "polygonnesting.h"

///////
// Example Implementation of VERTEX_TYPE, POLYGON_TYPE, VALUE_TYPE = float, and functors
///////

// simple 2D vertex
struct Vertex2D
{
    float x;
    float y;
};

using Polygon = std::vector<Vertex2D>;
using PolygonSet = std::vector<Polygon*>;

/// Computes the winding order of the vertices of a polygon
PolygonNesting<Polygon, Vertex2D>::VertexOrder ComputeVertexOrder(const Polygon& p) {
    float area = 0;
    //begin with edge to the first vertex and then iterate through all edges
    Vertex2D startVertex = p.back();
    Vertex2D endVertex;
    for (Polygon::const_iterator i = p.begin(); i != p.end(); ++i) {
        endVertex = *i;
        area += startVertex.x * endVertex.y - startVertex.y * endVertex.x;
        startVertex = endVertex;
    }

    return (area > 0) ? PolygonNesting<Polygon, Vertex2D>::VertexOrder::CCW : PolygonNesting<Polygon, Vertex2D>::VertexOrder::CW;
}

// functors for example implementation
auto getVertexOrder = [](const Polygon* p) -> PolygonNesting<Polygon, Vertex2D>::VertexOrder { return ComputeVertexOrder(*p); };
auto getNumVertices = [](const Polygon* p) -> size_t { return p->size(); };
auto getVertex = [](const Polygon* p, size_t i) -> const Vertex2D& { return (*p)[i]; };
auto getX = [](const Vertex2D& v) -> float { return v.x; };
auto getY = [](const Vertex2D& v) -> float { return v.y; };

///////
// End Example implementation
///////

void PrintPolygonSet(const PolygonSet& polySet)
{
    for (size_t i = 0; i < polySet.size(); ++i)
    {
        std::cout << "Polygon " << i << ":" << std::endl;
        for (auto& v : *polySet[i])
        {
            std::cout << "[" << v.x << " , " << v.y << "] "; 
        }
        std::cout << std::endl << std::endl;
    }
}

void PrintParents(const std::vector<size_t>& parents)
{
    std::cout << std::endl;
    for (size_t i = 0; i < parents.size(); ++i)
    {
        std::cout << "Parent(" << i << ") = ";
        if (parents[i] == static_cast<size_t>(-1))
        {
            std::cout << "none" << std::endl;
        }
        else
        {
            std::cout << parents[i] << std::endl;
        }    
    }
    std::cout << std::endl;
}

#define TEST_CONFIGURATION_1
#define TEST_CONFIGURATION_2
#define TEST_CONFIGURATION_3
#define TEST_CONFIGURATION_4
#define TEST_CONFIGURATION_5

int main()
{
    PolygonNesting<Polygon, Vertex2D> polygonNesting(getVertexOrder, getNumVertices, getVertex, getX, getY);

    std::vector<size_t> parentsResult;
#ifdef TEST_CONFIGURATION_1
    // test 1: two nested triangle in CCW orientation
    std::cout << std::endl << " TEST 1A " << std::endl << "-----" << std::endl;
    Polygon p1 = { Vertex2D{ 0.f, 0.f }, Vertex2D{2.f, 0.f}, Vertex2D{1.f, 2.f} }; 
    Polygon p2 = { Vertex2D{ 0.5f, 0.5f }, Vertex2D{1.5f, 0.5f}, Vertex2D{1.f, 1.5f} };

    PolygonSet polySet1 = { &p1, &p2 };

    PrintPolygonSet(polySet1);

    polygonNesting.Clear();
    for (auto p : polySet1)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents();

    // the parent of the second polygon should be the first
    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == 0);

    PrintParents(parentsResult);

    // test 2: invert orientation of both polygons to CW
    std::cout << std::endl << " TEST 1B " << std::endl << "-----" << std::endl;
    std::reverse(p1.begin(), p1.end());
    std::reverse(p2.begin(), p2.end());

    PrintPolygonSet(polySet1);

    polygonNesting.Clear();
    for (auto p : polySet1)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents();

    // the parent of the second polygon should be the first
    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == 0);

    PrintParents(parentsResult);

    // test 3: outer orientation is CCW, inner is CW 
    std::cout << std::endl << " TEST 1C " << std::endl << "-----" << std::endl;
    std::reverse(p1.begin(), p1.end());

    PrintPolygonSet(polySet1);

    polygonNesting.Clear();
    for (auto p : polySet1)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents();

    // the parent of the second polygon should be the first
    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == 0);

    PrintParents(parentsResult);
#endif

#ifdef TEST_CONFIGURATION_2
    // test 4: Figure 2.1 from technical report 
    std::cout << std::endl << " TEST 2 " << std::endl << "-----" << std::endl;
    Polygon p3 = { Vertex2D{ 0.f, 7.f }, Vertex2D{6.5f, 11.f}, Vertex2D{10.f, 11.5f}, Vertex2D{11.5f, 14.f}, Vertex2D{14.5f, 12.f}, Vertex2D{13.8f, 10.9f}, Vertex2D{12.f, 10.5f}, Vertex2D{10.2f, 8.5f}, Vertex2D{12.5f, 5.f}, Vertex2D{11.f, 2.f}, Vertex2D{7.5f, 0.8f}, Vertex2D{4.3f, 1.f}, Vertex2D{1.f, 3.f} }; 

    PolygonSet polySet2 = { &p3 };

    PrintPolygonSet(polySet2);

    polygonNesting.Clear();
    for (auto p : polySet2)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents();

    assert(parentsResult[0] == static_cast<size_t>(-1));

    PrintParents(parentsResult);
#endif

#ifdef TEST_CONFIGURATION_3
    // Figure 1.1 from technical report 
    std::cout << std::endl << " TEST 3 " << std::endl << "-----" << std::endl;
    Polygon p4 = { Vertex2D{ 22.f, 14.f }, Vertex2D{ 21.f, 11.f }, Vertex2D{ 18.75f, 11.5f }, Vertex2D{ 20.f, 8.f }, Vertex2D{ 23.f, 8.5f }, Vertex2D{ 24.f, 11.f } };
    Polygon p5 = { Vertex2D{ 10.f, 0.f }, Vertex2D{ 19.f, 4.f }, Vertex2D{ 18.5f, 9.f }, Vertex2D{ 11.5f, 10.f }, Vertex2D{ 11.f, 14.f }, Vertex2D{ 3.f, 12.f }, Vertex2D{ 2.f, 7.f }, Vertex2D{ 4.f, 3.f } };
    Polygon p6 = { Vertex2D{ 10.f, 12.f }, Vertex2D{ 12.f, 6.f }, Vertex2D{ 11.f, 2.f }, Vertex2D{ 5.f, 4.f }, Vertex2D{ 4.f, 7.f }, Vertex2D{ 6.f, 11.f } }; 
    Polygon p7 = { Vertex2D{ 6.f, 9.f }, Vertex2D{ 7.f, 7.f }, Vertex2D{ 9.5f, 8.5f }, Vertex2D{ 7.5f, 9.5f }, Vertex2D{ 7.f, 8.5f } };
    Polygon p8 = { Vertex2D{ 11.f, 7.5f }, Vertex2D{ 8.5f, 6.5f }, Vertex2D{ 11.f, 5.5f } };

    PolygonSet polySet3 = { &p4, &p5, &p6, &p7, &p8 };

    PrintPolygonSet(polySet3);

    polygonNesting.Clear();
    for (auto p : polySet3)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents(); 

    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == static_cast<size_t>(-1));
    assert(parentsResult[2] == 1);
    assert(parentsResult[3] == 2);
    assert(parentsResult[4] == 2);

    PrintParents(parentsResult);
#endif

#ifdef TEST_CONFIGURATION_4
    // Figure 3.2 from technical report 
    std::cout << std::endl << " TEST 4 " << std::endl << "-----" << std::endl;
    Polygon p9 = { Vertex2D{ 12.f, -5.f }, Vertex2D{ 12.f, -1.f }, Vertex2D{ 9.f, 2.f }, Vertex2D{ 5.f, 6.f }, Vertex2D{ -2.f, 5.f }, Vertex2D{ -5.f, 4.f }, Vertex2D{ -6.f, -5.f }, Vertex2D{ -2.f, -8.f }, Vertex2D{ 7.f, -4.f } };
    Polygon p10 = { Vertex2D{ 0.f, 0.f }, Vertex2D{ 3.f, -1.f }, Vertex2D{ 4.f, -4.f }, Vertex2D{ 1.f, -6.f }, Vertex2D{ -1.f, -5.f }, Vertex2D{ 0.f, -4.f } };
    Polygon p11 = { Vertex2D{ 2.f, 2.f }, Vertex2D{ 0.f, 3.f }, Vertex2D{ 3.f, 4.f }, Vertex2D{ 6.f, 1.f }, Vertex2D{ 1.5f, 0.5f } };
    Polygon p12 = { Vertex2D{ 2.f, -4.5 }, Vertex2D{ 2.f, -3.f }, Vertex2D{ 0.5f, -4.5f } }; 

    PolygonSet polySet4 = { &p9, &p10, &p11, &p12 };

    PrintPolygonSet(polySet4);

    polygonNesting.Clear();
    for (auto p : polySet4)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents();

    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == 0);
    assert(parentsResult[2] == 0);
    assert(parentsResult[3] == 1);

    PrintParents(parentsResult);
#endif

#ifdef TEST_CONFIGURATION_5
    // Two quads with collinear edges 
    std::cout << std::endl << " TEST 5 " << std::endl << "-----" << std::endl;
    Polygon p13 = { Vertex2D{ 6.f, 9.f }, Vertex2D{ 3.f, 9.f }, Vertex2D{ 0.f, 9.f }, Vertex2D{ 0.f, 6.f }, Vertex2D{ 0.f, 3.f }, Vertex2D{ 0.f, 0.f }, Vertex2D{ 3.f, 0.f },  Vertex2D{ 6.f, 0.f }, Vertex2D{ 9.f, 0.f }, Vertex2D{ 9.f, 3.f }, Vertex2D{ 9.f, 6.f }, Vertex2D{ 9.f, 9.f } };
    Polygon p14 = { Vertex2D{ 6.f, 3.f }, Vertex2D{ 3.f, 3.f }, Vertex2D{ 3.f, 6.f }, Vertex2D{ 6.f, 6.f } };

    PolygonSet polySet5 = { &p13, &p14 };

    PrintPolygonSet(polySet5);

    polygonNesting.Clear();
    for (auto p : polySet5)
    {
        polygonNesting.AddPolygon(p);
    }

    polygonNesting.ComputePolygonNesting();
    parentsResult = polygonNesting.GetParents(); 

    assert(parentsResult[0] == static_cast<size_t>(-1));
    assert(parentsResult[1] == 0);

    PrintParents(parentsResult);
#endif

    return 0;
}
