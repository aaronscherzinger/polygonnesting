#include "polygonnesting.h"

int main()
{
    Polygon p1 = { Vertex2D{ 0.f, 0.f }, Vertex2D{2.f, 0.f}, Vertex2D{1.f, 2.f} }; 
    Polygon p2 = { Vertex2D{ 0.5f, 0.5f }, Vertex2D{1.5f, 0.5f}, Vertex2D{1.f, 1.5f} };

    PolygonSet polySet = { &p1, &p2 };

    PolygonNesting(polySet);

    return 0;
}
