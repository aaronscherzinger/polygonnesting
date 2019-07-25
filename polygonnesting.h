#include <vector>

#include <algorithm>
#include <assert.h>

#ifdef POLYGONNESTING_PRINT_DEBUG
#include <iostream>
#endif

// simple 2D vertex
// #TODO: replace this by a more flexible implementation allowing arbitrary data and functors for retrieving X and Y coordinates
struct Vertex2D
{
    float x;
    float y;
};

// #TODO: replae by more flexible implementation??
using Polygon = std::vector<Vertex2D>;
using PolygonSet = std::vector<Polygon*>;

enum class VertexOrder
{
    CCW,    // counter-clockwise orientation
    CW      // clock-wise orientation
};

struct Subchain
{
    size_t polygon;                 ///< index into the polygon set for the polygon containing this subchain
    std::vector<size_t> vertices;   ///< indices into the polygon's vertex list

    size_t currentEdge;             ///< index into the polygon's vertex list for the left vertex of the edge in the subchain currently intersected by the vertical sweep line
};

/// Computes the winding order of the vertices of a polygon (#TODO: do this when computing the polygon instead of afterwards)
VertexOrder ComputeVertexOrder(const Polygon& p) {
    float area = 0;
    //begin with edge to the first vertex and then iterate through all edges
    Vertex2D startVertex = p.back();
    Vertex2D endVertex;
    for (Polygon::const_iterator i = p.begin(); i != p.end(); ++i) {
        endVertex = *i;
        area += startVertex.x * endVertex.y - startVertex.y * endVertex.x;
        startVertex = endVertex;
    }

    return (area > 0) ? VertexOrder::CCW : VertexOrder::CW;
}

/// Computes the orientation of three consecutive vertices 
float OrientationTest(const Vertex2D& pa, const Vertex2D& pb, const Vertex2D& pc) {
    float acx = pa.x - pc.x;
    float bcx = pb.x - pc.x;
    float acy = pa.y - pc.y;
    float bcy = pb.y - pc.y;

    return acx * bcy - acy * bcx;
}

/// Checks if vertex pc lies to the left of the directed line through vertices pa and pb or on the line
bool LeftTurnCollinear(const Vertex2D& pa, const Vertex2D& pb, const Vertex2D& pc) {
    return (OrientationTest(pa, pb, pc) >= 0);
}

/// Checks if vertex pc lies to the right of the directed line through vertices pa and pb or on the line
bool RightTurnCollinear(const Vertex2D& pa, const Vertex2D& pb, const Vertex2D& pc) {
    return (OrientationTest(pa, pb, pc) <= 0);
}

// #TODO: ordering of points: x1<x2 or if x1==x2: y1>y2

// #TODO: compute intersection of edge with vertical line


/// returns the cyclic successor index for a vertex index
size_t succ(size_t vertexIndex, const Polygon& polygon)
{
    return (vertexIndex + 1) % polygon.size();
}

/// returns the cyclic predecessor index for a vertex index
size_t pred(size_t vertexIndex, const Polygon& polygon)
{
    assert(polygon.size() > 2);
    
    return (vertexIndex == 0) ? polygon.size() - 1 : vertexIndex - 1;
}

// #TODO: template this by: Polygon Type, Vertex Type, coordinate Type
// #TODO: member functions to add polygon (by pointers) and clear everything
// #TODO: set functors that get x and y from a vertex, get a specific vertex by index from a polygon, get the number of vertices from a polygon
// #TODO: maybe optional optionally allow to directly use point operator on vertices
// #TODO: how to handle succ, pred, getting the vertex order?

void PolygonNesting(const PolygonSet& polygonSet)
{
    constexpr size_t INVALID_INDEX = static_cast<size_t>(-1);

    assert(!polygonSet.empty());

    // step 1: break down each polygon into subchains

    // preliminary step: compute the vertex winding order for each polygon
    // #TODO: this step should be done prior to performing the polygon nesting algorithm if possible (e.g., during assembly of the polygon vertices)
    std::vector<VertexOrder> vertexOrder(polygonSet.size());
    for (size_t i = 0; i < polygonSet.size(); ++i)
    {
        vertexOrder[i] = ComputeVertexOrder(*(polygonSet[i]));
    }

    // find subchains for all polygons
    // #TODO: this step could be performed in parallel, but this might only be beneficial for a large number of polygons
    std::vector<Subchain> subchains;
    for (size_t i = 0; i < polygonSet.size(); ++i)
    {
        const Polygon& currentPolygon = *(polygonSet[i]);

        assert(currentPolygon.size() > 2);

        // depending on the winding order, we need the left turn or right turn test to check for reflex vertices
        auto convexityTest = (vertexOrder[i] == VertexOrder::CCW) ? &LeftTurnCollinear : &RightTurnCollinear;

        // we look for the leftmost vertex as it is a guaranteed starting point of a subchain
        size_t leftMostVertex = 0;   
        float minX = currentPolygon[0].x;     
        for (size_t currentVertex = 1; currentVertex < currentPolygon.size(); ++currentVertex)
        {
            if (currentPolygon[currentVertex].x < minX)
            {
                leftMostVertex = currentVertex;
                minX = currentPolygon[currentVertex].x;
            }
        }

        // traverse vertices
        // subchain ends: 
        // - when the convex polygonal line ends (i.e., at a reflex vertex)
        // - when the next vertex would break the convex chain (i.e., connecting the last to the first vertex would lead to a reflex angle)
        // - when the next vertex would break the subchain (i.e., x-monotony of the chain)
        size_t subChainEndVertex = currentPolygon.size();   // end vertex of the current subchain (in the beginning undefined)
        bool subChainEnded = false;                         // flag that is se when one of the terminating conditions for the current subchain is met
        bool increaseX = true;                              // for checking x-monotony (in the beginning, x increases as we start with the leftmost vertex)

        // we terminate once a subchain ends where we started
        size_t currentVertex = leftMostVertex;
        Subchain currentSubchain = { };
        currentSubchain.polygon = i;
        while(subChainEndVertex != leftMostVertex)
        {
            size_t nextVertex = succ(currentVertex, currentPolygon);
            
            // check conditions for ending subchain
            if (    // 1. we are back at the beginning
                    (currentVertex == leftMostVertex && subChainEndVertex != currentPolygon.size())
                    // 2. current vertex is a reflex vertex - convex polygonal line ends, therefore also the subchain
                    ||  (!convexityTest(currentPolygon[pred(currentVertex, currentPolygon)], currentPolygon[currentVertex], currentPolygon[nextVertex]))
                    // 3. next vertex breaks the convex chain - we need to terminate the polygonal chain
                    ||  ((currentSubchain.vertices.size() > 1) && !convexityTest(currentPolygon[nextVertex], currentPolygon[currentSubchain.vertices[0]], currentPolygon[succ(currentSubchain.vertices[0], currentPolygon)]))
                    // 4. next vertex breaks the monotony
                    ||  ((increaseX && currentPolygon[nextVertex].x < currentPolygon[currentVertex].x) || (!increaseX && currentPolygon[nextVertex].x > currentPolygon[currentVertex].x))
                )
            {
                subChainEnded = true;
                subChainEndVertex = currentVertex;
            }

            currentSubchain.vertices.push_back(currentVertex);

            if (subChainEnded)
            {
                // subchain ended - insert into list of subchains
                if (!increaseX)
                {
                    std::reverse(currentSubchain.vertices.begin(), currentSubchain.vertices.end());
                }
                subchains.push_back(currentSubchain);

                // set new subchain start vertex and continue
                currentSubchain = { };
                currentSubchain.polygon = i;
                currentSubchain.vertices.push_back(currentVertex);
                currentSubchain.currentEdge = INVALID_INDEX;
               
                increaseX = (currentPolygon[nextVertex].x > currentPolygon[currentVertex].x);                
                subChainEnded = false;
            }

            currentVertex = nextVertex;
        }
    }

#ifdef POLYGONNESTING_PRINT_DEBUG
    for (auto& s : subchains)
    {
        std::cout << "Subchain: " << " polygon " << s.polygon << ", vertices: ";
        for (auto i : s.vertices)
        {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
#endif
 
    // step 2: sort the endpoints of all subchains

    assert(subchains.size() > 1);

    // for each chain, insert its start and end into a list that is then sorted
    // we mark the start and end by the index into the subchain's vertex list
    using Endpoint = std::pair<size_t, size_t>; // first index: subchain, second index: subchain vertex index
    std::vector<Endpoint> sortedEndpoints(subchains.size() * 2);

    for (size_t i = 0; i < subchains.size(); ++i)
    {
        sortedEndpoints[2*i] = std::make_pair(i, 0);
        sortedEndpoints[2*i+1] = std::make_pair(i, subchains[i].vertices.size() - 1);
    }

    std::sort(sortedEndpoints.begin(), sortedEndpoints.end(), 
        // comparison functor: 1. by x coordinate, 2. by y coordinate, 3. by index (i.e., if it is an endpoint or a start point)        
        [&](Endpoint a, Endpoint b) 
        {
            size_t polyA = subchains[a.first].polygon;
            size_t polyB = subchains[b.first].polygon;
            size_t vertA = subchains[a.first].vertices[a.second];
            size_t vertB = subchains[b.first].vertices[b.second];

            if ((*(polygonSet[polyA]))[vertA].x != (*(polygonSet[polyB]))[vertB].x)
            {
                return ((*(polygonSet[polyA]))[vertA].x < (*(polygonSet[polyB]))[vertB].x);
            }
            else if ((*(polygonSet[polyA]))[vertA].y != (*(polygonSet[polyB]))[vertB].y)
            {
                return ((*(polygonSet[polyA]))[vertA].y > (*(polygonSet[polyB]))[vertB].y);
            }
            else
            {
                // at this point, the two endpoints connect two subchains of the same polygon
                // #TODO: in the case of two starting / ending subchains, insert the one above the other first
                return vertA > vertB;
            }
        });

#ifdef POLYGONNESTING_PRINT_DEBUG
    std::cout << std::endl << "Sorting result: " << std::endl;
    for (auto& e : sortedEndpoints)
    {
        size_t poly = subchains[e.first].polygon;
        size_t vert = subchains[e.first].vertices[e.second];
        float x = (*(polygonSet[poly]))[vert].x;
        float y = (*(polygonSet[poly]))[vert].y;
        std::cout << "[" << x << ", " << y << "] (" << e.first << ", " << e.second << ") ";
    }
    std::cout << std::endl;
#endif

    // step 3: setup the initial situation for the plane sweep

    assert(sortedEndpoints.size() > 2);

    std::vector<size_t> orderedSubchains;   // #TODO: approximate size?

    // insert the two first subchains - make sure to take into account the "above" relationship when inserting them
    size_t s1 = sortedEndpoints[0].first;
    size_t s2 = sortedEndpoints[1].first;

    assert(subchains[s1].polygon == subchains[s2].polygon);
    assert(subchains[s1].vertices.size() > 1);
    assert(subchains[s2].vertices.size() > 1);

    // #TODO: above should already be used during sorting, making this check unnecessary
    size_t p = subchains[s1].polygon;

    if ((*(polygonSet[p]))[subchains[s1].vertices[1]].y > (*(polygonSet[p]))[subchains[s2].vertices[1]].y)
    {
        orderedSubchains.push_back(s1);
        orderedSubchains.push_back(s2);
    }
    else
    {
        orderedSubchains.push_back(s2);
        orderedSubchains.push_back(s1);
    }

    // we remember the current edge where the sweep line stands by its start vertex
    // (INVALID_INDEX means that the subchain has not been encountered yet)
    subchains[s1].currentEdge = 0;
    subchains[s2].currentEdge = 0;

    // #TODO: we already know that the polygon of subchains s1, s2 does not have a parent as it has the leftmost vertex
    // #TODO: create a graph for the parent relationships? how to represent it in the most efficient / useful manner?

    // setp 4: perform plane sweep from left to right

    // #TODO
    for (size_t sweepLineIndex = 2; sweepLineIndex < sortedEndpoints.size(); ++sweepLineIndex)
    {

    }
}

 
