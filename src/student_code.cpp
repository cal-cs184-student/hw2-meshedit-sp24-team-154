#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    if (points.size() <= 1) {
        return points;
    }
    vector<Vector2D> ans;
    ans.reserve(points.size()-1);

    for (int i = 0; i < points.size()-1; ++i) {
        ans.push_back((1-t)*points[i]+t*points[i+1]);
    }
    return ans;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    if (points.size() <= 1) {
        return points;
    }

    vector<Vector3D> ans;
    ans.reserve(points.size()-1);

    for (int i = 0; i < points.size()-1; ++i) {
        ans.push_back((1-t)*points[i]+t*points[i+1]);
    }

    return ans;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    vector<Vector3D> ans = points;
    while (ans.size() > 1) {
        ans = evaluateStep(ans, t);
    }
    return ans[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    vector<Vector3D> rows;
    unsigned long n = controlPoints.size();
    rows.reserve(n);

    for (int i = 0; i < n; ++i) {
        rows.push_back(evaluate1D(controlPoints[i], u));
    }

    Vector3D ans = evaluate1D(rows, v);
    return ans;
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    Vector3D norms(0,0,0);

    HalfedgeCIter h = halfedge();

    if (h->face()->isBoundary()) {
        h = h->twin()->next();
    }
    HalfedgeCIter h_orig = h;

    do {
        Vector3D p0 = h->vertex()->position;
        Vector3D p1 = h->next()->vertex()->position;
        Vector3D p2 = h->next()->next()->vertex()->position;

        Vector3D normal = cross(p1 - p0, p2 - p0);

        norms += normal;

        h = h->twin()->next();
    } while (h != h_orig);

    if (norms.norm() > 0) {
        norms.normalize();
    }

    return norms;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (e0->isBoundary()) {
        return e0;
    }

    HalfedgeIter h = e0->halfedge(), a = h->next(), a_twin = a->twin(), b = a->next(), b_twin = b->twin(), h_twin = h->twin(),
    c = h_twin->next(), c_twin = c->twin(), d = c->next(), d_twin = d->twin();

    VertexIter a_vertex = h->vertex(), h_twin_vertex = h_twin->vertex(), b_vertex = b->vertex(), d_vertex = d->vertex();
    EdgeIter a_edge = a->edge(), b_edge = b->edge(), c_edge = c->edge(), d_edge = d->edge();
    FaceIter h_face = h->face(), h_twin_face = h_twin->face();

    h->vertex() = b_vertex;
    h->next() = a; h->twin() = h_twin;  h->edge() = e0; h->face() = h_face;
    a->next() = b; a->twin() = d_twin; a->vertex() = d_vertex; a->edge() = d_edge; a->face() = h_face;
    b->next() = h; b->twin() = a_twin; b->vertex() = h_twin_vertex; b->edge() = a_edge; b->face() = h_face;
    h_twin->next() = c; h_twin->twin() = h; h_twin->vertex() = d_vertex; h_twin->edge() = e0; h_twin->face() = h_twin_face;
    c->next() = d; c->twin() = b_twin; c->vertex() = b_vertex; c->edge() = b_edge; c->face() = h_twin_face;
    d->next() = h_twin; d->twin() = c_twin; d->vertex() = a_vertex; d->edge() = c_edge; d->face() = h_twin_face;

    a_twin->twin() = b;
    a_twin->vertex() = b_vertex; a_twin->edge() = a_edge;
    b_twin->twin() = c;
    b_twin->vertex() = a_vertex; b_twin->edge() = b_edge;
    c_twin->twin() = d;
    c_twin->vertex() = d_vertex; c_twin->edge() = c_edge;
    d_twin->twin() = a;
    d_twin->vertex() = h_twin_vertex; d_twin->edge() = d_edge;

    a_vertex->halfedge() = d; h_twin_vertex->halfedge() = b; b_vertex->halfedge() = c; d_vertex->halfedge() = a;

    e0->halfedge() = h; a_edge->halfedge() = b; b_edge->halfedge() = c; c_edge->halfedge() = d; d_edge->halfedge() = a;

    h_face->halfedge() = h; h_twin_face->halfedge() = h_twin;

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if (e0->isBoundary()) {
        return VertexIter();
    }

    HalfedgeIter bc = e0->halfedge(), cb = bc->twin(), ca = bc->next(), ab = ca->next(), bd = cb->next(), dc = bd->next(),
    ca_outer = ca->twin(), ab_outer = ab->twin(), bd_outer = bd->twin(), dc_outer = dc->twin();

    EdgeIter e1 = ca->edge(), e2 = ab->edge(), e3 = bd->edge(), e4 = dc->edge();

    VertexIter b = bc->vertex(), c = cb->vertex(), d = dc->vertex(), a = ab->vertex();

    FaceIter bca = bc->face();
    FaceIter cbd = cb->face();

    VertexIter m = newVertex();
    EdgeIter e5 = newEdge();
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();
    FaceIter mca = newFace();
    FaceIter mbd = newFace();

    HalfedgeIter am = newHalfedge(), md = newHalfedge(), mc = newHalfedge(), ma = newHalfedge(), dm = newHalfedge(), cm = newHalfedge();

    e0->halfedge() = bc;
    e1->halfedge() = ca;
    e2->halfedge() = mc;
    e3->halfedge() = dm;
    e4->halfedge() = dc;
    e5->halfedge() = ab;
    e6->halfedge() = am;
    e7->halfedge() = bd;

    m->position = (c->position + b->position)/2;
    m->isNew = true;
    e5->isNew = true;
    e7->isNew = true;

    bc->setNeighbors(ca,cb,m,e0,bca);

    ab->setNeighbors(bc,md,a,e5,bca);
    cb->setNeighbors(bd,bc,c,e0,cbd);
    bd->setNeighbors(dc,cm,m,e7,cbd);

    am->setNeighbors(md,ma,b,e6,mca);
    md->setNeighbors(mc,ab,m,e5,mca);
    mc->setNeighbors(am,ab_outer,a,e2,mca);
    ma->setNeighbors(dm,am,m,e6,mbd);
    dm->setNeighbors(cm,bd_outer,b,e3,mbd);
    cm->setNeighbors(ma,bd,d,e7,mbd);
    ca_outer->setNeighbors(ca_outer->next(),ca,a,e1,ca_outer->face());
    ab_outer->setNeighbors(ab_outer->next(),mc,b,e2,ab_outer->face());
    bd_outer->setNeighbors(bd_outer->next(),dm,d,e3,bd_outer->face());
    dc_outer->setNeighbors(dc_outer->next(),dc,c,e4,dc_outer->face());

    b->halfedge() = am;

    a->halfedge() = mc;

    m->halfedge() = bc;

    mca->halfedge() = am;
    mbd->halfedge() = ma;

    return m;
  }


  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
        Vector3D sum(0.0, 0.0, 0.0);
        HalfedgeIter h = v->halfedge();
        float n = 0;
        do {
            sum += h->twin()->vertex()->position;
            n++;
            h = h->twin()->next();
        } while (h != v->halfedge());
        float u;
        if (n == 3) {
            u = 3.0;
        } else {
            u = 3.0/(8.0*n);
        }
        v->isNew = false;
        v->newPosition = (1.0 - n * u) * v->position + u * sum;
    }

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {

        HalfedgeIter h = e->halfedge();
        HalfedgeIter b = h->next()->next();
        HalfedgeIter a_twin = h->twin();
        HalfedgeIter d = a_twin->next()->next();

        e->isNew = false;
        e->newPosition = 3.0 / 8.0 * (h->vertex()->position + a_twin->vertex()->position) + 1.0 / 8.0 * (b->vertex()->position + d->vertex()->position);
    }

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {
        if (!e->halfedge()->vertex()->isNew && !e->halfedge()->twin()->vertex()->isNew) {
            VertexIter newV = mesh.splitEdge(e);
            newV->isNew = true;
            newV->position = e->newPosition;
        }

        if (e->isNew) {
            if (e->halfedge()->vertex()->isNew != e->halfedge()->twin()->vertex()->isNew) {
                mesh.flipEdge(e);
            }
        }
    }

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
        if (!v->isNew) {
            v->position = v->newPosition;
        }
        v->isNew = false;
    }
  }
}
