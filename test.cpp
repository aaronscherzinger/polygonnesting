#include <algorithm>
#include <iostream>

#define POLYGONNESTING_PRINT_DEBUG
#include "polygonnesting.h"

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

//#define TEST_CONFIGURATION_1
//#define TEST_CONFIGURATION_2
#define TEST_CONFIGURATION_3
#define TEST_CONFIGURATION_4
#define TEST_CONFIGURATION_5

//#TODO: once decided how the output should look like, add in asserts to check for the correct graph

int main()
{
#ifdef TEST_CONFIGURATION_1
    // test 1: two nested triangle in CCW orientation
    std::cout << std::endl << " TEST 1 " << std::endl << "-----" << std::endl;
    Polygon p1 = { Vertex2D{ 0.f, 0.f }, Vertex2D{2.f, 0.f}, Vertex2D{1.f, 2.f} }; 
    Polygon p2 = { Vertex2D{ 0.5f, 0.5f }, Vertex2D{1.5f, 0.5f}, Vertex2D{1.f, 1.5f} };

    PolygonSet polySet = { &p1, &p2 };

    PrintPolygonSet(polySet);

    PolygonNesting(polySet);

    // test 2: invert orientation of both polygons to CW
    std::cout << std::endl << " TEST 2 " << std::endl << "-----" << std::endl;
    std::reverse(p1.begin(), p1.end());
    std::reverse(p2.begin(), p2.end());

    PrintPolygonSet(polySet);

    PolygonNesting(polySet);

    // test 3: outer orientation is CCW, inner is CW 
    std::cout << std::endl << " TEST 3 " << std::endl << "-----" << std::endl;
    std::reverse(p1.begin(), p1.end());

    PrintPolygonSet(polySet);

    PolygonNesting(polySet);
#endif

#ifdef TEST_CONFIGURATION_2
    // test 4: Figure 2.1 from technical report 
    std::cout << std::endl << " TEST 4 " << std::endl << "-----" << std::endl;
    Polygon p1 = { Vertex2D{ 0.f, 7.f }, Vertex2D{6.5f, 11.f}, Vertex2D{10.f, 11.5f}, Vertex2D{11.5f, 14.f}, Vertex2D{14.5f, 12.f}, Vertex2D{13.8f, 10.9f}, Vertex2D{12.f, 10.5f}, Vertex2D{10.2f, 8.5f}, Vertex2D{12.5f, 5.f}, Vertex2D{11.f, 2.f}, Vertex2D{7.5f, 0.8f}, Vertex2D{4.3f, 1.f}, Vertex2D{1.f, 3.f} }; 

    PolygonSet polySet = { &p1 };

    PrintPolygonSet(polySet);

    PolygonNesting(polySet);
#endif

#ifdef TEST_CONFIGURATION_3
   // #TODO 
#endif

#ifdef TEST_CONFIGURATION_4
   // #TODO 
#endif

#ifdef TEST_CONFIGURATION_5
   // #TODO 
#endif

    return 0;
}
