import unittest
import random
import sys

sys.path.insert(0, ".")
from voronoi import (
    bisectorIntersection,
    findL1Bisector,
    samePoint,
    generateVoronoiPoints,
    cleanData,
    generateL1Voronoi,
    distance,
)


class TestBisectorIntersection(unittest.TestCase):
    def _s(self, x, y):
        return {"site": [x, y], "bisectors": []}

    def test_normal_shared_site(self):
        A = self._s(4, 6)
        B = self._s(3, 10)
        C1 = self._s(10, 6)
        C2 = self._s(10, 6)
        bAC = findL1Bisector(A, C1, 30, 30)
        bBC = findL1Bisector(B, C2, 30, 30)
        P = bisectorIntersection(bAC, bBC)
        self.assertIsInstance(P, list)
        self.assertEqual(len(P), 2)

    def test_non_shared_bisectors(self):
        A = self._s(0, 0)
        B = self._s(10, 0)
        C = self._s(0, 20)
        D = self._s(10, 20)
        bAB = findL1Bisector(A, B, 30, 30)
        bCD = findL1Bisector(C, D, 30, 30)
        P = bisectorIntersection(bAB, bCD)
        self.assertTrue(P is False or isinstance(P, list))

    def test_no_sites_key(self):
        zline = {"points": [[-4, 4], [200, 4]]}
        C = self._s(-8, -2)
        bBC = findL1Bisector(self._s(-1, 1), C, 200, 200)
        P = bisectorIntersection(zline, bBC)
        self.assertTrue(P is False or isinstance(P, list))


class TestGenerateVoronoiPoints(unittest.TestCase):
    def test_basic(self):
        pts = [[4, 6], [3, 10], [10, 6], [1, 2]]
        result = generateVoronoiPoints(pts, 30, 30, distance)
        self.assertEqual(len(result), 900)


class TestCleanData(unittest.TestCase):
    def test_nudge(self):
        data = [[4, 6], [6, 4]]
        cleanData(data)
        self.assertNotEqual(data[1], [6, 4])


class TestGenerateL1Voronoi(unittest.TestCase):
    def test_four_sites(self):
        sites = [[4, 6], [3, 10], [10, 6], [1, 2]]
        result = generateL1Voronoi(sites, 30, 30, nudgeData=True)
        self.assertEqual(len(result), 4)
        for site in result:
            self.assertIn("d", site)
            self.assertIn("neighbors", site)
            self.assertIn("polygonPoints", site)
            self.assertGreaterEqual(len(site["polygonPoints"]), 3)

    def test_16_random_sites(self):
        random.seed(42)
        sites = [
            [int(random.random() * 400), int(random.random() * 400)] for _ in range(16)
        ]
        result = generateL1Voronoi(sites, 400, 400, nudgeData=True)
        self.assertEqual(len(result), 16)
        for site in result:
            self.assertGreaterEqual(len(site["neighbors"]), 1)


class TestRandomStress(unittest.TestCase):
    @staticmethod
    def _random_normal(sharpness):
        return sum(random.random() for _ in range(sharpness)) / sharpness

    def _gen_sites(self, seed, n):
        random.seed(seed)
        return [
            [int(self._random_normal(2) * 400), int(self._random_normal(2) * 400)]
            for _ in range(n)
        ]

    def test_89_sites(self):
        sites = self._gen_sites(1, 89)
        result = generateL1Voronoi(sites, 400, 400)
        self.assertGreaterEqual(len(result), 1)
        for site in result:
            self.assertGreaterEqual(len(site["neighbors"]), 1)
            self.assertGreaterEqual(len(site["polygonPoints"]), 3)
            self.assertIn("d", site)
            self.assertTrue(site["d"].startswith("M "))
            self.assertTrue(site["d"].endswith(" Z"))

    def test_multi_seed_stress(self):
        for seed in range(20):
            sites = self._gen_sites(seed, 67)
            result = generateL1Voronoi(sites, 400, 400)
            for site in result:
                self.assertGreaterEqual(len(site["neighbors"]), 1)
                self.assertGreaterEqual(len(site["polygonPoints"]), 3)
                for bisector in site["bisectors"]:
                    pts = bisector["points"]
                    self.assertGreaterEqual(
                        len(pts), 2, f"seed {seed}: bisector has <2 points"
                    )
                    for k in range(len(pts) - 1):
                        self.assertFalse(
                            samePoint(pts[k], pts[k + 1]),
                            f"seed {seed}: zero-length segment in bisector",
                        )

    def test_intersection_correctness(self):
        for seed in range(10):
            sites_raw = self._gen_sites(seed, 30)
            site_objs = [{"site": p, "bisectors": []} for p in sites_raw]
            for i in range(len(site_objs)):
                for j in range(i + 1, len(site_objs)):
                    for k in range(len(site_objs)):
                        if k == i or k == j:
                            continue
                        A = site_objs[i]
                        B = site_objs[j]
                        C = site_objs[k]
                        try:
                            bAC = findL1Bisector(A, C, 400, 400)
                        except ValueError:
                            continue
                        try:
                            bBC = findL1Bisector(B, C, 400, 400)
                        except ValueError:
                            continue
                        P = bisectorIntersection(bAC, bBC)
                        if P is False:
                            continue
                        self.assertTrue(
                            self._point_on_bisector(P, bAC),
                            f"seed {seed}: P={P} not on bAC {A['site']}-{C['site']}",
                        )
                        self.assertTrue(
                            self._point_on_bisector(P, bBC),
                            f"seed {seed}: P={P} not on bBC {B['site']}-{C['site']}",
                        )

    def test_polygon_simple(self):
        for seed in range(10):
            sites = self._gen_sites(seed, 20)
            result = generateL1Voronoi(sites, 400, 400)
            for site in result:
                pts = site["polygonPoints"]
                n = len(pts)
                for a in range(n):
                    for c in range(a + 2, n):
                        if a == 0 and c == n - 1:
                            continue
                        b = (a + c) // 2
                        if self._segments_cross(
                            pts[a], pts[b], pts[b if b + 1 < n else 0], pts[c]
                        ):
                            self.fail(
                                f"seed {seed} self-intersecting polygon at site {site['site']}"
                            )

    @staticmethod
    def _segments_cross(p0, p1, p2, p3):
        d1 = (p1[0] - p0[0], p1[1] - p0[1])
        d2 = (p3[0] - p2[0], p3[1] - p2[1])
        denom = d2[1] * d1[0] - d2[0] * d1[1]
        if denom == 0:
            return False
        ua = (d2[0] * (p0[1] - p2[1]) - d2[1] * (p0[0] - p2[0])) / denom
        ub = (d1[0] * (p0[1] - p2[1]) - d1[1] * (p0[0] - p2[0])) / denom
        return 0 < ua < 1 and 0 < ub < 1

    @staticmethod
    def _point_on_bisector(pt, bisector):
        def point_on_segment(p, s0, s1):
            dx = s1[0] - s0[0]
            dy = s1[1] - s0[1]
            if dx == 0 and dy == 0:
                return p[0] == s0[0] and p[1] == s0[1]
            if abs(dx) > abs(dy):
                t = (p[0] - s0[0]) / dx
            elif abs(dx) < abs(dy):
                t = (p[1] - s0[1]) / dy
            else:
                t = (p[0] - s0[0]) / dx
            if not (0 <= t <= 1):
                return False
            return (
                abs(p[0] - (s0[0] + t * dx)) < 1e-9
                and abs(p[1] - (s0[1] + t * dy)) < 1e-9
            )

        for k in range(len(bisector["points"]) - 1):
            if point_on_segment(pt, bisector["points"][k], bisector["points"][k + 1]):
                return True
        return False


if __name__ == "__main__":
    unittest.main()
