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

// struct for the list of subchains
struct Subchain
{
    size_t polygon;                 ///< index into the polygon set for the polygon containing this subchain
    std::vector<size_t> vertices;   ///< indices into the polygon's vertex list

    size_t currentEdge;             ///< index into the polygon's vertex list for the left vertex of the edge in the subchain currently intersected by the vertical sweep line
};

// struct for the ordered list of endpoints
// an endpoint corresponds to the vertex of a specific polygon and connects two of it's subchains
struct Endpoint
{
    // indices of the subchains in the global list of subchains
    size_t subchains[2];
    // indices into the subchains' vertex lists 
    size_t subchainVertexIndices[2];
    // index of the polygon
    size_t polygon;
    // index into the polygon's vertex list for the corresponding vertex
    size_t polygonVertexIndex; 
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

class SubchainCompareFunctor
{

    public:

        SubchainCompareFunctor(const PolygonSet& polygonSet, std::vector<Subchain>& subchains, float sweepLineCoord)
        : m_polygonSet(polygonSet)
        , m_subchains(subchains)
        , m_sweepLineCoord(sweepLineCoord)
        {
        }

        SubchainCompareFunctor() = delete;

        SubchainCompareFunctor(const SubchainCompareFunctor& other) = default;
        SubchainCompareFunctor& operator=(const SubchainCompareFunctor& other) = default;

        void SetSweepLineCoord(float s) { m_sweepLineCoord = s; }

        /// compares two subchains given by their index in the set of subchains with respect to an x-coordinate of the sweep line
        /// NOTE: this will increment the current edge of the subchains in order to progress along the subchain according to the sweep line x-coordinate
        /// Make sure that consecutive calls to this function with the same set of subchains have mononously increasing values of sweepLineCoord
        bool operator() (const size_t& a, const size_t& b) const
        {
            if (a == b)
            {
                return false;
            }

            auto& polygonA = *(m_polygonSet[m_subchains[a].polygon]);
            auto& polygonB = *(m_polygonSet[m_subchains[b].polygon]);

            // we need to compare the subchains regarding their y-coordinate of the intersection with the sweep line
            while(polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].x < m_sweepLineCoord)
            {
                m_subchains[a].currentEdge++;
            }
            assert(m_subchains[a].currentEdge < m_subchains[a].vertices.size());
            // note: vertices [currentEdge - 1] and [currentEdge] correspond to the edge intersecting the sweep line
            float yCoordA = polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].y;
            if (polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].x > m_sweepLineCoord)
            {
                assert(m_subchains[a].currentEdge > 0);
                // sweep line intersects the edge and not the vertex
                // we can compute the intersection by linear interpolation
                float x1 = polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge - 1]].x;
                float y1 = polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge - 1]].y;

                float x2 = polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].x;
                float y2 = polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].y;

                // #TODO: doesn't this go wrong for vertical edges?
                assert(x1 != x2);
                assert(m_sweepLineCoord >= x1);
                assert(m_sweepLineCoord <= x2);

                float ratio = (m_sweepLineCoord - x1) / (x2 - x1);
                yCoordA = (1.f - ratio) * y1 + ratio * y2;
            }

            while(polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].x < m_sweepLineCoord)
            {
                m_subchains[b].currentEdge++;
            }
            assert(m_subchains[b].currentEdge < m_subchains[b].vertices.size());
            // note: vertices [currentEdge - 1] and [currentEdge] correspond to the edge intersecting the sweep line
            float yCoordB = polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].y;
            if (polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].x > m_sweepLineCoord)
            {
                assert(m_subchains[b].currentEdge > 0);
                // sweep line intersects the edge and not the vertex
                // we can compute the intersection by linear interpolation
                float x1 = polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge - 1]].x;
                float y1 = polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge - 1]].y;

                float x2 = polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].x;
                float y2 = polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].y;

                // #TODO: doesn't this go wrong for vertical edges?
                assert(x1 != x2);
                assert(m_sweepLineCoord >= x1);
                assert(m_sweepLineCoord <= x2);

                float ratio = (m_sweepLineCoord - x1) / (x2 - x1);
                yCoordB = (1.f - ratio) * y1 + ratio * y2;
            }

            if (yCoordA == yCoordB)
            {
                // if this happens, a shared endpoint of two subchains has been met
                assert(polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge]].x == m_sweepLineCoord);
                assert(polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge]].x == m_sweepLineCoord);
                assert(m_subchains[a].polygon == m_subchains[b].polygon);

                // this means that one of the following cases has happened:
                // 1. we encounter the endpoint of one subchain that is the starting point of another one
                // 2. we encounter the startpoint of two subchains
                // 3. we encounter the endpoint of two subchains

                // case 1: we always put the endpoint of one chain before the starting point of the next as this does not affect the y-order
                if (m_subchains[a].currentEdge == 0 && m_subchains[b].currentEdge > 0)
                {
                    return true;
                }
                else if (m_subchains[b].currentEdge == 0 && m_subchains[a].currentEdge > 0)
                {
                    return false;
                } 
                else if (m_subchains[a].currentEdge == 0 && m_subchains[b].currentEdge == 0)
                {
                    // case 2: we put the subchain with the higher y-coordinate in the next vertex first
                    assert(m_subchains[a].vertices.size() > 1);
                    assert(m_subchains[b].vertices.size() > 1);

                    return (polygonA[m_subchains[a].vertices[1]].y > polygonB[m_subchains[b].vertices[1]].y);
                }
                else 
                {
                    // case 3: we put the subchains with the higher y-coordinate in the preceding vertex first
                    return (polygonA[m_subchains[a].vertices[m_subchains[a].currentEdge - 1]].y > polygonB[m_subchains[b].vertices[m_subchains[b].currentEdge - 1]].y);
                }
            }
            else
            {
                return (yCoordA > yCoordB);
            }
        }

    private:

        const PolygonSet& m_polygonSet;
        std::vector<Subchain>& m_subchains;
        float m_sweepLineCoord;
};

// #TODO: template this by: Polygon Type, Vertex Type, coordinate Type?
// #TODO: member functions to add polygon (by pointers) and clear everything
// #TODO: generic implementation, check paper, or else:
// #TODO: set functors that get x and y from a vertex, get a specific vertex by index from a polygon, get the number of vertices from a polygon
// #TODO: maybe optional optionally allow to directly use point operator on vertices
// #TODO: how to handle succ, pred, getting the vertex order?

std::vector<size_t> PolygonNesting(const PolygonSet& polygonSet)
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

    // find subchains for all polygons and also construct the list of endpoints
    std::vector<Subchain> subchains;
    std::vector<Endpoint> endpoints;
    for (size_t i = 0; i < polygonSet.size(); ++i)
    {
        const Polygon& currentPolygon = *(polygonSet[i]);

        assert(currentPolygon.size() > 2);

        // depending on the winding order, we need the left turn or right turn test to check for reflex vertices
        auto convexityTest = (vertexOrder[i] == VertexOrder::CCW) ? &LeftTurnCollinear : &RightTurnCollinear;

        // we look for the leftmost vertex as it is a guaranteed starting point of a subchain
        // if two vertices have the same x-coordinate, we choose the one with the higher y-coordinate
        size_t leftMostVertex = 0;   
        for (size_t currentVertex = 1; currentVertex < currentPolygon.size(); ++currentVertex)
        {
            if ((currentPolygon[currentVertex].x < currentPolygon[leftMostVertex].x) 
                || ((currentPolygon[currentVertex].x == currentPolygon[leftMostVertex].x) 
                    && (currentPolygon[currentVertex].y > currentPolygon[leftMostVertex].y)))
            {
                leftMostVertex = currentVertex;
            }
        }

        // traverse vertices
        // subchain ends: 
        // - when the convex polygonal line ends (i.e., at a reflex vertex)
        // - when the next vertex would break the convex chain (i.e., connecting the last to the first vertex would lead to a reflex angle)
        // - when the next vertex would break the subchain (i.e., strict (!) x-monotony of the chain)
        // note: if there are multiple consecutive vertices with the same x-coordinate, we create a degenerated subchain
        size_t subChainEndVertex = currentPolygon.size();   // end vertex of the current subchain (in the beginning undefined)
        bool subChainEnded = false;                         // flag that is se when one of the terminating conditions for the current subchain is met
        bool increaseX = true;                              // for checking x-monotony (in the beginning, x increases as we start with the leftmost vertex)

        // current subchain
        Subchain currentSubchain = { };
        currentSubchain.polygon = i;
        currentSubchain.currentEdge = INVALID_INDEX;
        // current endpoint
        Endpoint currentEndpoint = { };
        currentEndpoint.polygon = i;
        currentEndpoint.polygonVertexIndex = leftMostVertex;

        size_t firstEndpointIndex = endpoints.size();   // index to the first endpoint created for this polygon, as we will need it later on

        // we terminate once a subchain ends where we started
        size_t currentVertex = leftMostVertex;
        while(subChainEndVertex != leftMostVertex)
        {
            size_t nextVertex = succ(currentVertex, currentPolygon);

            // #TODO: add handling of degenerated vertices: put all consecutive degenerated edges in a degenerated subchain
            
            // check conditions for ending subchain
            if (    // 1. we are back at the beginning
                    (currentVertex == leftMostVertex && subChainEndVertex != currentPolygon.size())
                    // 2. current vertex is a reflex vertex (or collinear) - convex polygonal line ends, therefore also the subchain
                    ||  (!convexityTest(currentPolygon[pred(currentVertex, currentPolygon)], currentPolygon[currentVertex], currentPolygon[nextVertex]))
                    // 3. next vertex breaks the convex chain - we need to terminate the polygonal chain
                    ||  ((currentSubchain.vertices.size() > 1) && !convexityTest(currentPolygon[nextVertex], currentPolygon[currentSubchain.vertices[0]], currentPolygon[succ(currentSubchain.vertices[0], currentPolygon)]))
                    // 4. next vertex breaks the strict monotony
                    ||  ((increaseX && currentPolygon[nextVertex].x <= currentPolygon[currentVertex].x) || (!increaseX && currentPolygon[nextVertex].x >= currentPolygon[currentVertex].x))
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

                // the current endpoint is the end vertex of the last chain and thus our start vertex
                size_t endpointIndex = increaseX ? 0 : (currentSubchain.vertices.size() - 1);
                assert(currentEndpoint.polygonVertexIndex == currentSubchain.vertices[endpointIndex]);
                currentEndpoint.subchains[1] = subchains.size() - 1;
                currentEndpoint.subchainVertexIndices[1] = endpointIndex;
                // append the endpoint to the list
                endpoints.push_back(currentEndpoint);

                size_t nextEndpointIndex = increaseX ? (currentSubchain.vertices.size() - 1) : 0;

                if (subChainEndVertex == leftMostVertex)
                {
                    // this is the last subchain - we need to put our information into the first endpoint
                    endpoints[firstEndpointIndex].subchains[0] = subchains.size() - 1;
                    endpoints[firstEndpointIndex].subchainVertexIndices[0]= nextEndpointIndex;
                    
                    assert(currentSubchain.vertices[nextEndpointIndex] == leftMostVertex);
                    assert(endpoints[firstEndpointIndex].polygonVertexIndex == currentSubchain.vertices[nextEndpointIndex]);
                }
                else
                {
                    // the new endpoint is our end vertex
                    currentEndpoint.subchains[0] = subchains.size() - 1;
                    currentEndpoint.subchainVertexIndices[0] = nextEndpointIndex;
                    currentEndpoint.polygonVertexIndex = currentVertex; 

                    // set new subchain start vertex and continue
                    currentSubchain = { };
                    currentSubchain.polygon = i;
                    currentSubchain.vertices.push_back(currentVertex);
                    currentSubchain.currentEdge = INVALID_INDEX;
                   
                    increaseX = (currentPolygon[nextVertex].x > currentPolygon[currentVertex].x);                
                    subChainEnded = false;
                }
            }

            currentVertex = nextVertex;
        }
    }

    assert(subchains.size() > 1);

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

    std::sort(endpoints.begin(), endpoints.end(), 
        // comparison functor: 1. by x coordinate, 2. by y coordinate       
        [&](Endpoint a, Endpoint b) 
        {
            size_t polyA = a.polygon;
            size_t polyB = b.polygon;
            size_t vertA = a.polygonVertexIndex;
            size_t vertB = b.polygonVertexIndex;

            assert((*(polygonSet[polyA]))[vertA].x != (*(polygonSet[polyB]))[vertB].x
                    || (*(polygonSet[polyA]))[vertA].y != (*(polygonSet[polyB]))[vertB].y);

            if ((*(polygonSet[polyA]))[vertA].x != (*(polygonSet[polyB]))[vertB].x)
            {
                return ((*(polygonSet[polyA]))[vertA].x < (*(polygonSet[polyB]))[vertB].x);
            }
            else
            {
                return ((*(polygonSet[polyA]))[vertA].y > (*(polygonSet[polyB]))[vertB].y);
            }
        });

#ifdef POLYGONNESTING_PRINT_DEBUG
    std::cout << std::endl << "Sorting result: " << std::endl;
    for (auto& e : endpoints)
    {
        size_t poly = e.polygon;
        size_t vert = e.polygonVertexIndex;
        float x = (*(polygonSet[poly]))[vert].x;
        float y = (*(polygonSet[poly]))[vert].y;
        std::cout << "[" << x << ", " << y << "] (" << e.subchains[0] << ", " << e.subchains[1] << ") ";
    }
    std::cout << std::endl;
#endif

    // step 3: setup the initial situation for the plane sweep

    std::vector<size_t> orderedSubchains;
    std::vector<std::vector<size_t> > orderedSubchainsForPolygons(polygonSet.size());
    std::vector<size_t> parents(polygonSet.size());
    std::fill(parents.begin(), parents.end(), INVALID_INDEX);

    // insert the two first subchains - make sure to take into account the "above" relationship when inserting them
    // #TODO: account for degenerated subchains
    size_t s1 = endpoints[0].subchains[0];
    size_t s2 = endpoints[0].subchains[1];

    assert(endpoints[0].polygon == subchains[s1].polygon);
    assert(subchains[s1].polygon == subchains[s2].polygon);
    assert(subchains[s1].vertices.size() > 1);
    assert(subchains[s2].vertices.size() > 1);

    size_t p = endpoints[0].polygon;

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

    orderedSubchainsForPolygons[p].push_back(orderedSubchains[0]);
    orderedSubchainsForPolygons[p].push_back(orderedSubchains[1]);

    // we remember the current edge where the sweep line stands by its start vertex
    // (the initial INVALID_INDEX means that a subchains has not been visited yet)
    subchains[s1].currentEdge = 0;
    subchains[s2].currentEdge = 0;

    // step 4: perform plane sweep from left to right

    SubchainCompareFunctor compare(polygonSet, subchains, 0.f);

    for (size_t sweepLineIndex = 1; sweepLineIndex < endpoints.size(); ++sweepLineIndex)
    {
        // get the next subchain endpoint v_i
        Endpoint& currentEndpoint = endpoints[sweepLineIndex];

        // curent sweep line x-coordinate
        float sweepLineCoord = (*(polygonSet[currentEndpoint.polygon]))[currentEndpoint.polygonVertexIndex].x;
        compare.SetSweepLineCoord(sweepLineCoord);

        // check if both subchains of this endpoint have already been visited - if so, remove them from the list of ordered subchains and continue
        if ((subchains[currentEndpoint.subchains[0]].currentEdge != INVALID_INDEX) && subchains[currentEndpoint.subchains[1]].currentEdge != INVALID_INDEX)
        {
                // find the two subchains and remove them from the list of ordered subchains
                auto firstIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), currentEndpoint.subchains[0], compare);
                assert(firstIterator != orderedSubchains.end());
                orderedSubchains.erase(firstIterator);
                
                auto secondIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), currentEndpoint.subchains[1], compare);
                assert(secondIterator != orderedSubchains.end());
                orderedSubchains.erase(secondIterator);

                // do the same for the list of ordered subchains for the specific polygon of these subchains
                auto firstPolygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), orderedSubchainsForPolygons[currentEndpoint.polygon].end(), currentEndpoint.subchains[0], compare);
                assert(firstPolygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());
                orderedSubchainsForPolygons[currentEndpoint.polygon].erase(firstPolygonIterator);

                auto secondPolygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), orderedSubchainsForPolygons[currentEndpoint.polygon].end(), currentEndpoint.subchains[1], compare);
                assert(secondPolygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());
                orderedSubchainsForPolygons[currentEndpoint.polygon].erase(secondPolygonIterator);

                continue;
        }

        // step 4(a)
        
        // we can either have one subchain inserted or none of them - check which is the case
        size_t s1 = currentEndpoint.subchains[0];
        size_t s2 = currentEndpoint.subchains[1];
        bool s1Inserted = subchains[s1].currentEdge != INVALID_INDEX;
        bool s2Inserted = subchains[s2].currentEdge != INVALID_INDEX;

        // we will insert the new subchains later and in order to use the binary search, we set their current edge to 0
        // #TODO: account for degenerated subchains - do not set their edge index!
        if (!s1Inserted)
        {
            subchains[s1].currentEdge = 0;
        }
        if (!s2Inserted)
        {
            subchains[s2].currentEdge = 0;
        }

        // find position in subchain order
        // this will either give us the position of the inserted subchain or the position where we need to insert if none is inserted
        // #TODO: account for degenerated subchain - if one of the subchains for this endpoint is degenerated, use the other one
        auto positionIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s2Inserted ? s2 : s1, compare);

        // step 4(b)
        auto neighbor = orderedSubchains.end();
        if (!orderedSubchains.empty() && (positionIterator != orderedSubchains.begin()))
        {
            // there is at least one subchain above the current subchain - get the direct neighbor
            neighbor = positionIterator - 1;

            size_t polygon = subchains[*neighbor].polygon;
            
            // only proceed if the neighbor is a different polygon
            if (polygon != currentEndpoint.polygon)
            {
                // check in the list of subchains for this specific polygon to to know how many subchains are above the current subchain
                auto polygonIterator = std::lower_bound(orderedSubchainsForPolygons[polygon].begin(), orderedSubchainsForPolygons[polygon].end(), *neighbor, compare);
                assert(polygonIterator != orderedSubchainsForPolygons[polygon].end());
                size_t numAbove = std::distance(orderedSubchainsForPolygons[polygon].begin(), polygonIterator) + 1;
                // if odd: set the parent to the polygon, else set it to the parent of the polygon
                if (numAbove % 2 == 1)
                {
                    parents[currentEndpoint.polygon] = polygon;
                }
                else
                {
                    parents[currentEndpoint.polygon] = parents[polygon];
                }
            }
        }

        // step 4(c)
        // find the position in the polygon ordering as well
        auto polygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), 
                orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s2Inserted ? s2 : s1, compare);

        // in case one of the subchains has already been visited we replace it by the other one in both orderings
        // #TODO: account for degenerated subchains - only insert the non-degenerated subchain!
        if (s1Inserted)
        {
            assert(positionIterator != orderedSubchains.end());
            assert(polygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());

            *positionIterator = s2;
            *polygonIterator = s2;
        }
        else if (s2Inserted)
        {
            assert(positionIterator != orderedSubchains.end());
            assert(polygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());

            *positionIterator = s1;
            *polygonIterator = s1;
        }
        else
        {
            // otherwise, we insert both in the ordering
            orderedSubchains.insert(positionIterator, s1);
            orderedSubchainsForPolygons[currentEndpoint.polygon].insert(polygonIterator, s1);

            positionIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s2, compare);
            polygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), 
                orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s2, compare);

            orderedSubchains.insert(positionIterator, s2);
            orderedSubchainsForPolygons[currentEndpoint.polygon].insert(polygonIterator, s2);
        }
    }

    // #TODO: use better representation for result...
    return parents;
}

 
