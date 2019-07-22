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

int main()
{
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

    // #TODO: more interesting example - maybe one from the paper?

    return 0;
}
