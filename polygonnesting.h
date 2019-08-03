#include <vector>

#include <algorithm>
#include <assert.h>
#include <functional>
#include <type_traits>
#include <unordered_map>

#ifdef POLYGONNESTING_PRINT_DEBUG
#include <iostream>
#endif

/**
 * Class for computing the nesting structure of multiple polygons.
 *
 * In order to use this class, you must provide several functors to process the template types.
 * Note that VALUE_TYPE must be a floating point type.
 *
 * The functors are:
 *      PolygonNesting::VertexOrder GetVertexOrderFunctor(const POLYGON_TYPE*) - used once initially to get the vertex order of a polygon
 *      size_t GetNumVerticesFunctor(const POLYGON_TYPE*) - used to retrieve the number of vertices in a polygon
 *      const VERTEX_TYPE& GetVertexFunctor(const POLYGON_TYPE*,size_t) - used to get a const reference to the vertex with a specific index in the polygon
 *      VALUE_TYPE GetXFunctor(const VERTEX_TYPE&) - used to get the x-coordinate of a vertex
 *      VALUE_TYPE GetYFunctor(const VERTEX_TYPE&) - used to get the y-coordinate of a vertex
 *
 *      Note that the GetVertexFunctor, GetXFunctor, and GetYFunctor may be called many times (so they should be implemented efficiently).
 *
 * In order to use the class, add polygons with the AddPolygon() method and then compute the nesting structure with ComputePolygonNesting().
 * Use GetParent() to retrieve the parent of a specific polygon after computing the nesting structure.
 * In order to clear the data (e.g., to process a different set of polygons), call the Clear()-method.
 */
template<typename POLYGON_TYPE, typename VERTEX_TYPE, typename VALUE_TYPE = float>
class PolygonNesting
{
    static_assert(std::is_floating_point<VALUE_TYPE>::value, "VALUE_TYPE must be floating point type");

public:
    
    // enum class for vertex order of polygons
    enum class VertexOrder
    {
        CCW,    // counter-clockwise orientation
        CW      // clock-wise orientation
    };

    // alias names for functors that need to be provided for the template types
    using GetVertexOrderFunctor = std::function<VertexOrder(const POLYGON_TYPE*)>;
    using GetNumVerticesFunctor = std::function<size_t(const POLYGON_TYPE*)>;
    using GetVertexFunctor = std::function<const VERTEX_TYPE&(const POLYGON_TYPE*,size_t)>;
    using GetXFunctor = std::function<VALUE_TYPE(const VERTEX_TYPE&)>;
    using GetYFunctor = std::function<VALUE_TYPE(const VERTEX_TYPE&)>;

    PolygonNesting(GetVertexOrderFunctor getVertexOrder, GetNumVerticesFunctor getNumVertices, GetVertexFunctor getVertex, GetXFunctor getX, GetYFunctor getY)
    : mf_GetVertexOrder(getVertexOrder)
    , mf_GetNumVertices(getNumVertices)
    , mf_GetVertex(getVertex)
    , mf_GetX(getX)
    , mf_GetY(getY)
    {
    } 
      
    PolygonNesting() = delete;

    /// Adds a polygon to the polygon set
    void AddPolygon(const POLYGON_TYPE* polygon)
    {
        m_polygonSet.push_back(polygon);
    } 

    /// Computes the polygon nesting of the added polygons
    void ComputePolygonNesting();

    /// returns the parent polygon of a polygon if any, else nullptr
    const POLYGON_TYPE* GetParent(const POLYGON_TYPE* polygon) const
    {
        auto it = m_parents.find(polygon);
        if (it != m_parents.end())
        {
            return it->second;
        }

        return nullptr;
    }

    /// Clears the added polygon list and the polygon nesting results
    void Clear()
    {
        m_polygonSet.clear();
        m_parents.clear();
    }

private:

    //  ************************
    //  internal data structures
    //  ************************

    constexpr static size_t INVALID_INDEX = static_cast<size_t>(-1);

    // struct for the list of subchains
    struct Subchain
    {
        size_t polygon;                 ///< index into the polygon set for the polygon containing this subchain
        std::vector<size_t> vertices;   ///< indices into the polygon's vertex list

        size_t currentEdge;             ///< index into the polygon's vertex list for the left vertex of the edge in the subchain currently intersected by the vertical sweep line

        bool degenerate;                ///< is set to true if the subchain contains only vertical edges
    };

    // endpoint struct for the ordered list of endpoints
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

    //  ****************
    //  helper functions
    //  ****************

    /// returns the cyclic successor index for a vertex index
    size_t succ(size_t vertexIndex, const POLYGON_TYPE* polygon)
    {
        return (vertexIndex + 1) % mf_GetNumVertices(polygon);
    }

    /// returns the cyclic predecessor index for a vertex index
    size_t pred(size_t vertexIndex, const POLYGON_TYPE* polygon)
    {
        assert(mf_GetNumVertices(polygon) > 2);
        
        return (vertexIndex == 0) ? mf_GetNumVertices(polygon) - 1 : vertexIndex - 1;
    }

    /// Computes the orientation of three consecutive vertices 
    VALUE_TYPE OrientationTest(const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc) {
        VALUE_TYPE acx = mf_GetX(pa) - mf_GetX(pc);
        VALUE_TYPE bcx = mf_GetX(pb) - mf_GetX(pc);
        VALUE_TYPE acy = mf_GetY(pa) - mf_GetY(pc);
        VALUE_TYPE bcy = mf_GetY(pb) - mf_GetY(pc);

        return acx * bcy - acy * bcx;
    }

    /// Checks if vertex pc lies to the left of the directed line through vertices pa and pb or on the line
    bool LeftTurnCollinear(const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc) {
        return (OrientationTest(pa, pb, pc) >= 0);
    }

    /// Checks if vertex pc lies to the right of the directed line through vertices pa and pb or on the line
    bool RightTurnCollinear(const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc) {
        return (OrientationTest(pa, pb, pc) <= 0);
    }

    //  ************
    //  data members
    //  ************

    // set of polygons that have been added
    std::vector<const POLYGON_TYPE*> m_polygonSet;

    // map containing the parent of each polygon (if it has one)
    std::unordered_map<const POLYGON_TYPE*, const POLYGON_TYPE*> m_parents;

    //  ***************
    //  functor members
    //  ***************

    GetVertexOrderFunctor   mf_GetVertexOrder;
    GetNumVerticesFunctor   mf_GetNumVertices;
    GetVertexFunctor        mf_GetVertex;
    GetXFunctor             mf_GetX;
    GetYFunctor             mf_GetY;

    // nested class that is used as a comparator
    class SubchainCompareFunctor
    {
        public:

            SubchainCompareFunctor(const std::vector<const POLYGON_TYPE*>& polygonSet, std::vector<Subchain>& subchains, VALUE_TYPE sweepLineCoord,
                    GetVertexFunctor getVertex, GetXFunctor getX, GetYFunctor getY)
            : m_polygonSet(polygonSet)
            , m_subchains(subchains)
            , m_sweepLineCoord(sweepLineCoord)
            , mf_GetVertex(getVertex)
            , mf_GetX(getX)
            , mf_GetY(getY)
            {
            }

            SubchainCompareFunctor() = delete;

            void SetSweepLineCoord(VALUE_TYPE s) { m_sweepLineCoord = s; }

            /// compares two subchains given by their index in the set of subchains with respect to an x-coordinate of the sweep line
            /// NOTE: this will increment the current edge of the subchains in order to progress along the subchain according to the sweep line x-coordinate
            /// Make sure that consecutive calls to this function with the same set of subchains have mononously increasing values of sweepLineCoord
            bool operator() (const size_t& a, const size_t& b) const
            {
                if (a == b)
                {
                    return false;
                }

                auto& GetVertex = mf_GetVertex;
                auto& GetX = mf_GetX;
                auto& GetY = mf_GetY;

                const POLYGON_TYPE* polygonA = m_polygonSet[m_subchains[a].polygon];
                const POLYGON_TYPE* polygonB = m_polygonSet[m_subchains[b].polygon];

                // we need to compare the subchains regarding their y-coordinate of the intersection with the sweep line
                VALUE_TYPE yCoordA = GetSweepLineIntersection(polygonA, const_cast<Subchain&>(m_subchains[a]), m_sweepLineCoord);
                VALUE_TYPE yCoordB = GetSweepLineIntersection(polygonB, const_cast<Subchain&>(m_subchains[b]), m_sweepLineCoord);

                if (yCoordA == yCoordB)
                {
                    // if this happens, a shared endpoint of two subchains has been met
                    assert(GetX(GetVertex(polygonA, m_subchains[a].vertices[m_subchains[a].currentEdge])) == m_sweepLineCoord);
                    assert(GetX(GetVertex(polygonB, m_subchains[b].vertices[m_subchains[b].currentEdge])) == m_sweepLineCoord);
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

                        return (GetY(GetVertex(polygonA, m_subchains[a].vertices[1])) > GetY(GetVertex(polygonB, m_subchains[b].vertices[1])));
                    }
                    else 
                    {
                        // case 3: we put the subchains with the higher y-coordinate in the preceding vertex first
                        return (GetY(GetVertex(polygonA, m_subchains[a].vertices[m_subchains[a].currentEdge - 1])) > GetY(GetVertex(polygonB, m_subchains[b].vertices[m_subchains[b].currentEdge - 1])));
                    }
                }
                else
                {
                    return (yCoordA > yCoordB);
                }
            }

        private:

            // helper function to get the intersection between the subchain and the sweepline
                   
            VALUE_TYPE GetSweepLineIntersection(const POLYGON_TYPE* polygon, Subchain& subchain, VALUE_TYPE sweepLineCoord) const
            {         
                auto& GetVertex = mf_GetVertex;
                auto& GetX = mf_GetX;
                auto& GetY = mf_GetY;

                assert(GetX(GetVertex(polygon, subchain.vertices.front())) <= sweepLineCoord);
                assert(GetX(GetVertex(polygon, subchain.vertices.back())) >= sweepLineCoord);

                // advance current edge until it intersects the sweep line
                while(GetX(GetVertex(polygon, subchain.vertices[subchain.currentEdge])) < sweepLineCoord)
                {
                    subchain.currentEdge++;
                }

                // note: vertices [currentEdge - 1] and [currentEdge] correspond to the edge intersecting the sweep line
                size_t edge = subchain.currentEdge;
                assert(edge < subchain.vertices.size());

                VALUE_TYPE yCoordA = GetY(GetVertex(polygon, subchain.vertices[edge]));

                if (GetX(GetVertex(polygon, subchain.vertices[edge])) > sweepLineCoord)
                {
                    assert(edge > 0);
                    // sweep line intersects the edge and not the vertex
                    // we can compute the intersection by linear interpolation
                    VALUE_TYPE x1 = GetX(GetVertex(polygon, subchain.vertices[edge - 1]));
                    VALUE_TYPE y1 = GetY(GetVertex(polygon, subchain.vertices[edge - 1]));

                    VALUE_TYPE x2 = GetX(GetVertex(polygon, subchain.vertices[edge]));
                    VALUE_TYPE y2 = GetY(GetVertex(polygon, subchain.vertices[edge]));

                    assert(x1 != x2);
                    assert(sweepLineCoord >= x1);
                    assert(sweepLineCoord <= x2);

                    VALUE_TYPE ratio = (sweepLineCoord - x1) / (x2 - x1);
                    yCoordA = (static_cast<VALUE_TYPE>(1) - ratio) * y1 + ratio * y2;
                }

                return yCoordA;
            }

            const std::vector<const POLYGON_TYPE*>& m_polygonSet;
            std::vector<Subchain>& m_subchains;
            VALUE_TYPE m_sweepLineCoord;

            // functors
            PolygonNesting::GetVertexFunctor        mf_GetVertex;
            PolygonNesting::GetXFunctor             mf_GetX;
            PolygonNesting::GetYFunctor             mf_GetY;
    };
};

// polygon nesting algorithm implementation
template<typename POLYGON_TYPE, typename VERTEX_TYPE, typename VALUE_TYPE>
void PolygonNesting<POLYGON_TYPE, VERTEX_TYPE, VALUE_TYPE>::ComputePolygonNesting()
{
    // some alias names to make the code more readable
    auto& GetVertexOrder = mf_GetVertexOrder;
    auto& GetNumVertices = mf_GetNumVertices;
    auto& GetVertex = mf_GetVertex;
    auto& GetX = mf_GetX;
    auto& GetY = mf_GetY;

    // all parents are initialized as nullptr
    m_parents.clear();
    m_parents.reserve(m_polygonSet.size());
    for (auto p : m_polygonSet)
    {
        m_parents[p] = nullptr;
    }

    if (m_polygonSet.empty())
    {
        return;
    }

    // find subchains for all polygons and also construct the list of endpoints
    std::vector<Subchain> subchains;
    std::vector<Endpoint> endpoints;

    // this will be set depending on the vertex order of a polygon
    std::function<bool(const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc)> convexityTest;
    
    for (size_t i = 0; i < m_polygonSet.size(); ++i)
    {
        const POLYGON_TYPE* currentPolygon = m_polygonSet[i];
        size_t numVertices = GetNumVertices(currentPolygon);

        assert(GetNumVertices(currentPolygon) > 2);

        // depending on the winding order, we need the left turn or right turn test to check for reflex vertices
        if (GetVertexOrder(currentPolygon) == VertexOrder::CCW)
        {
            convexityTest = [this](const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc) -> bool { return LeftTurnCollinear(pa, pb, pc); };
        }
        else
        { 
            convexityTest = [this](const VERTEX_TYPE& pa, const VERTEX_TYPE& pb, const VERTEX_TYPE& pc) -> bool { return RightTurnCollinear(pa, pb, pc); };
        }

        // we look for the leftmost vertex as it is a guaranteed starting point of a subchain
        // if two vertices have the same x-coordinate, we choose the one with the higher y-coordinate
        size_t leftMostVertex = 0;   
        for (size_t currentVertex = 1; currentVertex < numVertices; ++currentVertex)
        {
            if (    (GetX(GetVertex(currentPolygon, currentVertex)) < GetX(GetVertex(currentPolygon, leftMostVertex))) 
                 || (   (GetX(GetVertex(currentPolygon, currentVertex)) == GetX(GetVertex(currentPolygon, leftMostVertex))) 
                     && (GetY(GetVertex(currentPolygon, currentVertex)) > GetY(GetVertex(currentPolygon, leftMostVertex)))))
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
        enum class SubchainDirection
        {
            LEFT,       // subchain's x-coordinate is strictly increasing
            RIGHT       // subchain's x-coordinate is strictly decreasing
        };

        size_t subChainEndVertex = numVertices;             // end vertex of the current subchain (in the beginning undefined)
        bool subChainEnded = false;                         // flag that is se when one of the terminating conditions for the current subchain is met
        SubchainDirection currentDirection = SubchainDirection::LEFT;

        // current subchain
        Subchain currentSubchain = { };
        currentSubchain.polygon = i;
        currentSubchain.currentEdge = INVALID_INDEX;
        currentSubchain.degenerate = false;
        // current endpoint
        Endpoint currentEndpoint = { };
        currentEndpoint.polygon = i;
        currentEndpoint.polygonVertexIndex = leftMostVertex;

        size_t firstEndpointIndex = endpoints.size();   // index to the first endpoint created for this polygon, as we will need it later on

        currentSubchain.vertices.push_back(leftMostVertex);
        size_t currentVertex = succ(leftMostVertex, currentPolygon);
        // we terminate once a subchain ends where we started
        while(subChainEndVertex != leftMostVertex)
        {
            assert(currentSubchain.vertices.size() > 0);

            size_t nextVertex = succ(currentVertex, currentPolygon);

            // if we just started a subchain, we check its direction
            if (currentSubchain.vertices.size() == 1)
            {
                currentSubchain.degenerate = false;
                VALUE_TYPE firstX = GetX(GetVertex(currentPolygon, currentSubchain.vertices[0]));
                VALUE_TYPE secondX = GetX(GetVertex(currentPolygon, currentVertex));
                if (firstX < secondX)
                {
                    currentDirection = SubchainDirection::LEFT;
                }
                else if (firstX > secondX)
                {
                    currentDirection = SubchainDirection::RIGHT;
                }
                else
                {
                    currentSubchain.degenerate = true;
                    // we still set the direction to reverse the vertex order within the subchain if necessary
                    VALUE_TYPE firstY = GetY(GetVertex(currentPolygon, currentSubchain.vertices[0]));
                    VALUE_TYPE secondY = GetY(GetVertex(currentPolygon, currentVertex));
                    if (firstY > secondY)
                    {
                        currentDirection = SubchainDirection::LEFT;
                    }
                    else
                    {
                        assert(firstY != secondY);
                        currentDirection = SubchainDirection::RIGHT;
                    }
                }
            }

            // check conditions for ending subchain
            if (    // 1. we are back at the beginning
                    (currentVertex == leftMostVertex && subChainEndVertex != numVertices)
                    // 2. we have a degenerate subchain and the next vertex would not have the same x-coordinate
                    || (currentSubchain.degenerate && (GetX(GetVertex(currentPolygon, nextVertex)) != GetX(GetVertex(currentPolygon, currentVertex))))
                    // 3. not degenerate, current vertex is a reflex vertex (or collinear) - convex polygonal line ends, therefore also the subchain
                    ||  (!currentSubchain.degenerate && 
                        !convexityTest(GetVertex(currentPolygon, pred(currentVertex, currentPolygon)), GetVertex(currentPolygon, currentVertex), GetVertex(currentPolygon, nextVertex)))
                    // 4. next vertex breaks the convex chain - we need to terminate the polygonal chain
                    ||  (!currentSubchain.degenerate && 
                        (currentSubchain.vertices.size() > 1) 
                        && !convexityTest(GetVertex(currentPolygon, nextVertex), GetVertex(currentPolygon, currentSubchain.vertices[0]), GetVertex(currentPolygon, succ(currentSubchain.vertices[0], currentPolygon))))
                    // 5. next vertex breaks the strict monotony
                    ||  (!currentSubchain.degenerate &&
                        (((currentDirection == SubchainDirection::LEFT) && GetX(GetVertex(currentPolygon, nextVertex)) <= GetX(GetVertex(currentPolygon, currentVertex))) 
                        || ((currentDirection == SubchainDirection::RIGHT) && GetX(GetVertex(currentPolygon, nextVertex)) >= GetX(GetVertex(currentPolygon, currentVertex)))))
                )
            {
                subChainEnded = true;
                subChainEndVertex = currentVertex;
            }

            currentSubchain.vertices.push_back(currentVertex);

            if (subChainEnded)
            {
                // subchain ended - insert into list of subchains
                if (currentDirection == SubchainDirection::RIGHT)
                {
                    std::reverse(currentSubchain.vertices.begin(), currentSubchain.vertices.end());
                }
                subchains.push_back(currentSubchain);

                // the current endpoint is the end vertex of the last chain and thus our start vertex
                size_t endpointIndex = (currentDirection == SubchainDirection::LEFT) ? 0 : (currentSubchain.vertices.size() - 1);
                assert(currentEndpoint.polygonVertexIndex == currentSubchain.vertices[endpointIndex]);
                currentEndpoint.subchains[1] = subchains.size() - 1;
                currentEndpoint.subchainVertexIndices[1] = endpointIndex;
                // append the endpoint to the list
                endpoints.push_back(currentEndpoint);

                size_t nextEndpointIndex = (currentDirection == SubchainDirection::LEFT) ? (currentSubchain.vertices.size() - 1) : 0;

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
        if (s.degenerate)
        {
            std::cout << "(degenerate)";
        }
        std::cout << std::endl;
    }
#endif

    // step 2: sort the endpoints of all subchains

    std::sort(endpoints.begin(), endpoints.end(), 
        // comparison functor: 1. by x coordinate, 2. by y coordinate       
        [&](Endpoint a, Endpoint b) 
        {
            const POLYGON_TYPE* polyA = m_polygonSet[a.polygon];
            const POLYGON_TYPE* polyB = m_polygonSet[b.polygon];
            size_t vertA = a.polygonVertexIndex;
            size_t vertB = b.polygonVertexIndex;

            assert(GetX(GetVertex(polyA, vertA)) != GetX(GetVertex(polyB, vertB))
                   || GetY(GetVertex(polyA, vertA)) != GetY(GetVertex(polyB, vertB)));

            if (GetX(GetVertex(polyA, vertA)) != GetX(GetVertex(polyB, vertB)))
            {
                return (GetX(GetVertex(polyA, vertA)) < GetX(GetVertex(polyB, vertB)));
            }
            else
            {
                return (GetX(GetVertex(polyA, vertA)) > GetX(GetVertex(polyB, vertB)));
            }
        });

#ifdef POLYGONNESTING_PRINT_DEBUG
    std::cout << std::endl << "Sorting result: " << std::endl;
    for (auto& e : endpoints)
    {
        const POLYGON_TYPE* poly = m_polygonSet[e.polygon];
        size_t vert = e.polygonVertexIndex;
        VALUE_TYPE x = GetX(GetVertex(poly, vert));
        VALUE_TYPE y = GetY(GetVertex(poly, vert));
        std::cout << "[" << x << ", " << y << "] (" << e.subchains[0] << ", " << e.subchains[1] << ") ";
    }
    std::cout << std::endl;
#endif

    // step 3: setup the initial situation for the plane sweep

    std::vector<size_t> orderedSubchains;
    std::vector<std::vector<size_t> > orderedSubchainsForPolygons(m_polygonSet.size());

    // insert the two first subchains - make sure to take into account the "above" relationship when inserting them
    size_t s1 = endpoints[0].subchains[0];
    size_t s2 = endpoints[0].subchains[1];

    assert(endpoints[0].polygon == subchains[s1].polygon);
    assert(subchains[s1].polygon == subchains[s2].polygon);
    assert(subchains[s1].vertices.size() > 1);
    assert(subchains[s2].vertices.size() > 1);

    size_t p = endpoints[0].polygon;

    // only insert non-degenerate subchains
    if (subchains[s1].degenerate)
    {
        orderedSubchains.push_back(s2);
    }
    else if (subchains[s2].degenerate)
    {
        orderedSubchains.push_back(s1);
    }
    else
    {
        // insert the subchain above the other one first
        if (GetY(GetVertex(m_polygonSet[p], subchains[s1].vertices[1])) > GetY(GetVertex(m_polygonSet[p], subchains[s2].vertices[1])))
        {
            orderedSubchains.push_back(s1);
            orderedSubchains.push_back(s2);
        }
        else
        {
            orderedSubchains.push_back(s2);
            orderedSubchains.push_back(s1);
        }
    }

    assert(orderedSubchains.size() > 0); // there is at least one non-degenerate subchain sine we start at the top left vertex

    orderedSubchainsForPolygons[p].push_back(orderedSubchains[0]);
    if (orderedSubchains.size() > 1)
    {
        orderedSubchainsForPolygons[p].push_back(orderedSubchains[1]);
    }

    // we remember the current edge where the sweep line stands by its start vertex
    // (the initial INVALID_INDEX means that a subchains has not been visited yet)
    if (!subchains[s1].degenerate)
    {
        subchains[s1].currentEdge = 0;
    }
    
    if (!subchains[s2].degenerate)
    {
        subchains[s2].currentEdge = 0;
    }

    // step 4: perform plane sweep from left to right

    SubchainCompareFunctor compare(m_polygonSet, subchains, 0.f, mf_GetVertex, mf_GetX, mf_GetY);

    for (size_t sweepLineIndex = 1; sweepLineIndex < endpoints.size(); ++sweepLineIndex)
    {
        // get the next subchain endpoint v_i
        Endpoint& currentEndpoint = endpoints[sweepLineIndex];
        // get the subchains connected to this endpoint
        size_t s1 = currentEndpoint.subchains[0];
        size_t s2 = currentEndpoint.subchains[1];
        bool s1Inserted = (subchains[s1].currentEdge != INVALID_INDEX);
        bool s2Inserted = (subchains[s2].currentEdge != INVALID_INDEX);
        bool s1Degenerate = subchains[s1].degenerate;
        bool s2Degenerate = subchains[s2].degenerate;

        // curent sweep line x-coordinate
        VALUE_TYPE sweepLineCoord = GetX(GetVertex(m_polygonSet[currentEndpoint.polygon], currentEndpoint.polygonVertexIndex));
            
        compare.SetSweepLineCoord(sweepLineCoord);

        // check if both subchains of this endpoint have already been visited (or are degenerate) - if so, remove them from the list of ordered subchains and continue
        if ((s1Degenerate || s1Inserted) && (s2Degenerate || s2Inserted))
        {
            if (!s1Degenerate)
            {
                auto firstIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s1, compare);
                assert(firstIterator != orderedSubchains.end());
                orderedSubchains.erase(firstIterator);

                auto firstPolygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s1, compare);
                assert(firstPolygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());
                orderedSubchainsForPolygons[currentEndpoint.polygon].erase(firstPolygonIterator);
            }
            
            if (!s2Degenerate)
            {
                auto secondIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s2, compare);
                assert(secondIterator != orderedSubchains.end());
                orderedSubchains.erase(secondIterator);

                auto secondPolygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s2, compare);
                assert(secondPolygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());
                orderedSubchainsForPolygons[currentEndpoint.polygon].erase(secondPolygonIterator);
            }

            continue;
        }

        // step 4(a)
        
        // we can either have one subchain inserted or none of them - check which is the case
        // we need to use the non-degenerate subchain if one is degenerate
        size_t subchainToSearch = s1;
        if (s1Degenerate || s2Inserted)
        {
            subchainToSearch = s2;
        }

        // we will insert the new subchains later and in order to use the binary search, we set their current edge to 0
        if (!s1Degenerate && !s1Inserted)
        {
            subchains[s1].currentEdge = 0;
        }

        if (!s2Degenerate && !s2Inserted)
        {
            subchains[s2].currentEdge = 0;
        }

        // find position in subchain order
        // this will either give us the position of the inserted subchain or the position where we need to insert if none is inserted
        auto positionIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), subchainToSearch, compare);

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
                    m_parents[m_polygonSet[currentEndpoint.polygon]] = m_polygonSet[polygon];
                }
                else
                {
                    m_parents[m_polygonSet[currentEndpoint.polygon]] = m_parents[m_polygonSet[polygon]];
                }
            }
        }

        // step 4(c)
        // find the position in the polygon ordering as well
        auto polygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), 
                orderedSubchainsForPolygons[currentEndpoint.polygon].end(), subchainToSearch, compare);

        // in case one of the subchains has already been visited we replace it by the other one in both orderings
        if (s1Inserted && !s2Degenerate)
        {
            assert(positionIterator != orderedSubchains.end());
            assert(polygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());

            *positionIterator = s2;
            *polygonIterator = s2;
        }
        else if (s2Inserted && !s1Degenerate)
        {
            assert(positionIterator != orderedSubchains.end());
            assert(polygonIterator != orderedSubchainsForPolygons[currentEndpoint.polygon].end());

            *positionIterator = s1;
            *polygonIterator = s1;
        }
        else
        {
            // otherwise, we insert the subchains in the ordering
            if (!s1Degenerate && !s1Inserted)
            {
                positionIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s1, compare);
                polygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), 
                    orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s1, compare);

                orderedSubchains.insert(positionIterator, s1);
                orderedSubchainsForPolygons[currentEndpoint.polygon].insert(polygonIterator, s1);
            }

            if (!s2Degenerate && !s2Inserted)
            {
                positionIterator = std::lower_bound(orderedSubchains.begin(), orderedSubchains.end(), s2, compare);
                polygonIterator = std::lower_bound(orderedSubchainsForPolygons[currentEndpoint.polygon].begin(), 
                    orderedSubchainsForPolygons[currentEndpoint.polygon].end(), s2, compare);

                orderedSubchains.insert(positionIterator, s2);
                orderedSubchainsForPolygons[currentEndpoint.polygon].insert(polygonIterator, s2);
            }
        }
    }
}
