/*
 * PolygonLoops.cc
 *
 * See PolygonLoops.h for API and behavior.
 */
#include <cstdint>
#include "PolygonLoops.h"

#include <math.h>
#include <string.h>   // memcpy
#include <map>
#include <set>
#include <algorithm>

namespace libGDSII {

namespace {

/*---------------------------------------------*/
/* bitwise (exact) keying for doubles          */
/*---------------------------------------------*/
static unsigned long long DoubleBits(double x)
{
  unsigned long long u = 0ULL;
  /* memcpy is strictly-aliasing safe */
  memcpy(&u, &x, sizeof(double));
  return u;
}

static long long RoundToLL(double v)
{
  /* portable rounding (works for negative values too) */
  return (v >= 0.0) ? (long long)(v + 0.5) : (long long)(v - 0.5);
}

/*---------------------------------------------*/
/* Key type usable in std::map / std::set      */
/*---------------------------------------------*/
#include <cmath>

struct Key {
  std::int64_t a;
  std::int64_t b;
  Key() : a(0), b(0) {}
  Key(std::int64_t aa, std::int64_t bb) : a(aa), b(bb) {}
};

static bool operator<(const Key& k1, const Key& k2) {
  if (k1.a < k2.a) return true;
  if (k1.a > k2.a) return false;
  return k1.b < k2.b;
}

static bool operator==(const Key& k1, const Key& k2) {
  return (k1.a == k2.a) && (k1.b == k2.b);
}

static Key MakeKey(const XY& p, double tol) {
  if (tol > 0.0) {
    const std::int64_t ix = (std::int64_t) llround(p.x / tol);
    const std::int64_t iy = (std::int64_t) llround(p.y / tol);
    return Key(ix, iy);
  }
  // exact mode (bitwise)
  return Key((std::int64_t)DoubleBits(p.x), (std::int64_t)DoubleBits(p.y));
}

static bool SamePoint(const XY& p1, const XY& p2, double tol)
{
  return MakeKey(p1, tol) == MakeKey(p2, tol);
}

static void RemoveConsecutiveDuplicates(Loop& pts, double tol)
{
  if (pts.empty()) return;
  Loop out;
  out.reserve(pts.size());
  out.push_back(pts[0]);
  for (size_t i = 1; i < pts.size(); ++i)
    if (!SamePoint(out.back(), pts[i], tol))
      out.push_back(pts[i]);
  pts.swap(out);
}

static void StripClosingDuplicate(Loop& ring, double tol)
{
  if (ring.size() >= 2 && SamePoint(ring.front(), ring.back(), tol))
    ring.pop_back();
}

static bool Collinear(const XY& a, const XY& b, const XY& c, double eps)
{
  const double abx = b.x - a.x;
  const double aby = b.y - a.y;
  const double acx = c.x - a.x;
  const double acy = c.y - a.y;

  const double cross = abx * acy - aby * acx;

  /* scale-aware tolerance */
  const double scale = std::max(1.0,
                        std::max(fabs(abx), std::max(fabs(aby),
                        std::max(fabs(acx), fabs(acy)))));

  return fabs(cross) <= eps * scale * scale;
}

} // anonymous namespace

/*---------------------------------------------*/
/* Public helpers                              */
/*---------------------------------------------*/
Loop FlatXYToLoop(const std::vector<double>& flatXY)
{
  Loop out;
  out.reserve(flatXY.size() / 2);
  for (size_t i = 0; i + 1 < flatXY.size(); i += 2)
    out.push_back(XY(flatXY[i], flatXY[i + 1]));
  return out;
}

std::vector<double> LoopToFlatXY(const Loop& loop)
{
  std::vector<double> out;
  out.reserve(2 * loop.size());
  for (size_t i = 0; i < loop.size(); ++i)
  {
    out.push_back(loop[i].x);
    out.push_back(loop[i].y);
  }
  return out;
}

double SignedArea(const Loop& ring)
{
  const size_t n = ring.size();
  if (n < 3) return 0.0;

  double a = 0.0;
  for (size_t i = 0; i < n; ++i)
  {
    const size_t j = (i + 1) % n;
    a += ring[i].x * ring[j].y - ring[j].x * ring[i].y;
  }
  return 0.5 * a;
}

void EnsureCCW(Loop& ring)
{
  if (SignedArea(ring) < 0.0)
    std::reverse(ring.begin(), ring.end());
}

void EnsureCW(Loop& ring)
{
  if (SignedArea(ring) > 0.0)
    std::reverse(ring.begin(), ring.end());
}

/*---------------------------------------------*/
/* SeparatePolygonLoops                         */
/*---------------------------------------------*/
PolygonLoopsResult SeparatePolygonLoops(const Loop& pathIn, double tol)
{
  PolygonLoopsResult res;

  Loop path = pathIn;
  RemoveConsecutiveDuplicates(path, tol);
  StripClosingDuplicate(path, tol);

  if (path.size() < 3)
    return res;

  Loop active;
  std::vector<Key> activeKeys;
  std::map<Key, size_t> indexByKey; /* vertex key -> index in active */

  LoopList extracted;

  bool hasPrev = false;
  Key prevK;

  for (size_t ip = 0; ip < path.size(); ++ip)
  {
    const XY& p = path[ip];
    const Key k = MakeKey(p, tol);

    if (hasPrev && k == prevK)
      continue;
    prevK = k;
    hasPrev = true;

    std::map<Key, size_t>::iterator it = indexByKey.find(k);
    if (it != indexByKey.end())
    {
      const size_t j = it->second;

      Loop loop;
      for (size_t i = j; i < active.size(); ++i)
        loop.push_back(active[i]);
      loop.push_back(p); /* closes at p */

      RemoveConsecutiveDuplicates(loop, tol);
      StripClosingDuplicate(loop, tol);

      if (loop.size() >= 3)
        extracted.push_back(loop);

      /* Collapse active to the revisited vertex */
      active.resize(j + 1);
      activeKeys.resize(j + 1);

      /* rebuild map */
      indexByKey.clear();
      for (size_t i = 0; i < activeKeys.size(); ++i)
        indexByKey[activeKeys[i]] = i;
    }
    else
    {
      indexByKey[k] = active.size();
      active.push_back(p);
      activeKeys.push_back(k);
    }
  }

  LoopList candidates = extracted;

  RemoveConsecutiveDuplicates(active, tol);
  StripClosingDuplicate(active, tol);
  if (active.size() >= 3)
    candidates.push_back(active);

  if (candidates.empty())
    return res;

  /* choose outer as max |area| */
  size_t outerIdx = 0;
  double best = fabs(SignedArea(candidates[0]));
  for (size_t i = 1; i < candidates.size(); ++i)
  {
    const double a = fabs(SignedArea(candidates[i]));
    if (a > best)
    {
      best = a;
      outerIdx = i;
    }
  }

  res.outer = candidates[outerIdx];
  EnsureCCW(res.outer);

  for (size_t i = 0; i < candidates.size(); ++i)
  {
    if (i == outerIdx) continue;
    Loop h = candidates[i];
    EnsureCW(h);
    res.holes.push_back(h);
  }

  return res;
}

PolygonLoopsResult SeparatePolygonLoops(const std::vector<double>& flatXY, double tol)
{
  return SeparatePolygonLoops(FlatXYToLoop(flatXY), tol);
}

/*---------------------------------------------*/
/* RemoveSharedCollinearVerticesInHoles         */
/*---------------------------------------------*/
void RemoveSharedCollinearVerticesInHoles(const Loop& outer, LoopList& holes, double tol, double eps)
{
  /* Count which vertices are present in >=2 DIFFERENT loops */
  std::map<Key, int> counts;

  {
    std::set<Key> uniq;
    for (size_t i = 0; i < outer.size(); ++i)
      uniq.insert(MakeKey(outer[i], tol));
    for (std::set<Key>::iterator it = uniq.begin(); it != uniq.end(); ++it)
      counts[*it] += 1;
  }

  for (size_t h = 0; h < holes.size(); ++h)
  {
    std::set<Key> uniq;
    for (size_t i = 0; i < holes[h].size(); ++i)
      uniq.insert(MakeKey(holes[h][i], tol));
    for (std::set<Key>::iterator it = uniq.begin(); it != uniq.end(); ++it)
      counts[*it] += 1;
  }

  std::set<Key> shared;
  for (std::map<Key, int>::iterator it = counts.begin(); it != counts.end(); ++it)
    if (it->second >= 2)
      shared.insert(it->first);

  /* Clean each hole */
  for (size_t h = 0; h < holes.size(); ++h)
  {
    Loop v = holes[h];
    RemoveConsecutiveDuplicates(v, tol);
    StripClosingDuplicate(v, tol);

    if (v.size() < 3)
    {
      holes[h].clear();
      continue;
    }

    bool changed = true;
    while (changed && v.size() >= 3)
    {
      changed = false;

      const size_t n = v.size();
      std::vector<char> removeMask(n, 0);

      for (size_t i = 0; i < n; ++i)
      {
        const Key ki = MakeKey(v[i], tol);
        if (shared.find(ki) == shared.end())
          continue;

        const XY& prev = v[(i + n - 1) % n];
        const XY& cur  = v[i];
        const XY& next = v[(i + 1) % n];

        if (Collinear(prev, cur, next, eps))
          removeMask[i] = 1;
      }

      Loop newV;
      newV.reserve(n);
      for (size_t i = 0; i < n; ++i)
        if (!removeMask[i])
          newV.push_back(v[i]);

      if (newV.size() >= 3 && newV.size() != v.size())
      {
        v.swap(newV);
        changed = true;
      }
    }

    EnsureCW(v);
    holes[h] = v;
  }

  /* Remove any holes that collapsed to <3 vertices */
  LoopList filtered;
  for (size_t h = 0; h < holes.size(); ++h)
    if (holes[h].size() >= 3)
      filtered.push_back(holes[h]);
  holes.swap(filtered);
}

} // namespace libGDSII
