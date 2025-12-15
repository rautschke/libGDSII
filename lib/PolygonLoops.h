#ifndef LIBGDSII_POLYGON_LOOPS_H
#define LIBGDSII_POLYGON_LOOPS_H

/*
 * PolygonLoops.h
 *
 * Geometry helpers to extend libGDSII polygon handling to support
 * "flattened" polygons that encode holes as interior closed sub-loops.
 *
 * The input format matches libGDSII's Entity::XY representation:
 *   flatXY = [x0,y0,x1,y1,...] (no requirement that last==first)
 *
 * Output rings are OPEN (no duplicated closing point).
 *
 * Portability notes:
 *  - Uses only the C/C++ standard library.
 *  - Written in a style that is easy to port to C++03/11.
 */

#include <vector>

namespace libGDSII {

struct XY
{
  double x;
  double y;

  XY() : x(0.0), y(0.0) {}
  XY(double xx, double yy) : x(xx), y(yy) {}
};

typedef std::vector<XY>      Loop;
typedef std::vector<Loop>    LoopList;

struct PolygonLoopsResult
{
  Loop     outer;   // open ring (>=3 points)
  LoopList holes;   // open rings (>=3 points each)
};

/** Convert flat (x0,y0,x1,y1,...) to a Loop (same order). */
Loop FlatXYToLoop(const std::vector<double>& flatXY);

/** Convert a Loop to flat (x0,y0,x1,y1,...) */
std::vector<double> LoopToFlatXY(const Loop& loop);

/**
 * Separate outer loop and interior loops (holes) from an arbitrary
 * "vertex walk" that may revisit earlier vertices to close sub-loops.
 *
 * tol:
 *   - tol == 0.0 => exact vertex matching (bitwise double equality).
 *   - tol  > 0.0 => tolerance matching by quantization (round(x/tol), round(y/tol)).
 *
 * Outer ring is chosen as the loop with maximum |area|.
 * Orientation is normalized: outer CCW, holes CW.
 */
PolygonLoopsResult SeparatePolygonLoops(const Loop& path, double tol /*=0.0*/);
PolygonLoopsResult SeparatePolygonLoops(const std::vector<double>& flatXY, double tol /*=0.0*/);

/**
 * Remove "shared-but-collinear" vertices from holes:
 * if a hole vertex is used by >=2 rings (outer or other holes) AND lies on
 * a straight segment within that hole, it is removed.
 *
 * This fixes the exact failure mode you showed:
 *   holes include a shared junction point on a straight edge
 *   (e.g. (-28.337, 19.829)) that should not be part of the hole ring.
 *
 * eps: collinearity tolerance (scale-aware).
 */
void RemoveSharedCollinearVerticesInHoles(const Loop& outer,
                                          LoopList& holes,
                                          double tol /*=0.0*/,
                                          double eps /*=1e-12*/);

/** Signed area (shoelace); positive => CCW. */
double SignedArea(const Loop& ring);

/** Ensure ring orientation. */
void EnsureCCW(Loop& ring);
void EnsureCW(Loop& ring);

} // namespace libGDSII

#endif // LIBGDSII_POLYGON_LOOPS_H
