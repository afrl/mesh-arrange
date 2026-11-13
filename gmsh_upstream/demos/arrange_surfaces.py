import gmsh
import math
import sys


def add_geo_box(xmin, ymin, zmin, xmax, ymax, zmax, lc):
    pts = [
        gmsh.model.geo.addPoint(xmin, ymin, zmin, lc),
        gmsh.model.geo.addPoint(xmax, ymin, zmin, lc),
        gmsh.model.geo.addPoint(xmax, ymax, zmin, lc),
        gmsh.model.geo.addPoint(xmin, ymax, zmin, lc),
        gmsh.model.geo.addPoint(xmin, ymin, zmax, lc),
        gmsh.model.geo.addPoint(xmax, ymin, zmax, lc),
        gmsh.model.geo.addPoint(xmax, ymax, zmax, lc),
        gmsh.model.geo.addPoint(xmin, ymax, zmax, lc),
    ]
    lns = [
        gmsh.model.geo.addLine(pts[0], pts[1]),
        gmsh.model.geo.addLine(pts[1], pts[2]),
        gmsh.model.geo.addLine(pts[2], pts[3]),
        gmsh.model.geo.addLine(pts[3], pts[0]),
        gmsh.model.geo.addLine(pts[4], pts[5]),
        gmsh.model.geo.addLine(pts[5], pts[6]),
        gmsh.model.geo.addLine(pts[6], pts[7]),
        gmsh.model.geo.addLine(pts[7], pts[4]),
        gmsh.model.geo.addLine(pts[0], pts[4]),
        gmsh.model.geo.addLine(pts[1], pts[5]),
        gmsh.model.geo.addLine(pts[2], pts[6]),
        gmsh.model.geo.addLine(pts[3], pts[7]),
    ]
    faces = []
    loops = [
        [-lns[0], -lns[3], -lns[2], -lns[1]],
        [lns[4], lns[5], lns[6], lns[7]],
        [lns[0], lns[9], -lns[4], -lns[8]],
        [lns[1], lns[10], -lns[5], -lns[9]],
        [lns[2], lns[11], -lns[6], -lns[10]],
        [lns[3], lns[8], -lns[7], -lns[11]],
    ]
    for lp in loops:
        cl = gmsh.model.geo.addCurveLoop(lp)
        faces.append(gmsh.model.geo.addPlaneSurface([cl]))
    sl = gmsh.model.geo.addSurfaceLoop(faces)
    return gmsh.model.geo.addVolume([sl])


def compute_tet_volume_and_min_sicn(vol_tags):
    total_vol = 0.0
    min_sicn = 1.0
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = {}
    for i, tag in enumerate(node_tags):
        coords[int(tag)] = (
            node_coords[3 * i],
            node_coords[3 * i + 1],
            node_coords[3 * i + 2],
        )
    for vtag in vol_tags:
        elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(3, vtag)
        for etype, etags, enodes in zip(elem_types, elem_tags, elem_nodes):
            if etype != 4:
                continue
            quals = gmsh.model.mesh.getElementQualities(etags, "minSICN")
            if len(quals):
                min_sicn = min(min_sicn, min(quals))
            for k in range(0, len(enodes), 4):
                p0 = coords[int(enodes[k])]
                p1 = coords[int(enodes[k + 1])]
                p2 = coords[int(enodes[k + 2])]
                p3 = coords[int(enodes[k + 3])]
                ux, uy, uz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
                vx, vy, vz = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
                wx, wy, wz = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]
                det = (
                    ux * (vy * wz - vz * wy)
                    - uy * (vx * wz - vz * wx)
                    + uz * (vx * wy - vy * wx)
                )
                total_vol += abs(det) / 6.0
    return total_vol, min_sicn


gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Verbosity", 2)

# ---------------------------------------------------------------------------
# 1. Priority boolean arrangement + void cut (toolType >= 0) + RVE clip (< 0)
# ---------------------------------------------------------------------------
gmsh.model.add("arrange_boolean_clip")
v_mat = add_geo_box(-1.2, -1.2, -1.2, 1.2, 1.2, 1.2, 0.45)
v_inc = add_geo_box(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.35)
v_void = add_geo_box(-0.2, -0.2, -0.2, 0.2, 0.2, 0.2, 0.25)
v_rve = add_geo_box(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.40)
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

out_vols, out_map = gmsh.model.mesh.arrangeSurfaces(
    [(3, v_mat), (3, v_inc)],
    [(3, v_void), (3, v_rve)],
    objectPriorities=[0, 10],
    toolTypes=[1, -1],
    toolPriorities=[20, 30],
    removeOverlaps=False,
    optimizeSurfaces=False,
)

mat_children = [t for (d, t) in out_map[0] if d == 3]
inc_children = [t for (d, t) in out_map[1] if d == 3]
gmsh.model.addPhysicalGroup(3, mat_children, 1, "Matrix")
gmsh.model.addPhysicalGroup(3, inc_children, 2, "Inclusion")

gmsh.model.mesh.generate(3)
gmsh.model.mesh.optimize("", True)

vol_mat, sicn_mat = compute_tet_volume_and_min_sicn(mat_children)
vol_inc, sicn_inc = compute_tet_volume_and_min_sicn(inc_children)
assert math.isclose(vol_mat, 7.0, rel_tol=1e-6), f"Expected 7.0, got {vol_mat}"
assert math.isclose(vol_inc, 0.936, rel_tol=1e-6), f"Expected 0.936, got {vol_inc}"
assert min(sicn_mat, sicn_inc) > 0.0, f"Non-positive SICN: {min(sicn_mat, sicn_inc)}"

# ---------------------------------------------------------------------------
# 2. Equal-priority K=3 multi-material harmonic medial overlap removal
# ---------------------------------------------------------------------------
gmsh.model.add("arrange_harmonic_overlap")
v0 = add_geo_box(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.45)
v1 = add_geo_box(-0.55, -0.35, -0.40, 0.35, 0.35, 0.40, 0.32)
v2 = add_geo_box(-0.20, -0.55, -0.35, 0.55, 0.25, 0.45, 0.32)
v3 = add_geo_box(-0.30, -0.25, -0.55, 0.30, 0.55, 0.55, 0.32)
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

out_vols2, out_map2 = gmsh.model.mesh.arrangeSurfaces(
    [(3, v0), (3, v1), (3, v2), (3, v3)],
    [],
    objectPriorities=[0, 10, 10, 10],
    removeOverlaps=True,
    optimizeSurfaces=True,
)

all_child_tags = []
for idx, name in enumerate(["Matrix", "PhaseA", "PhaseB", "PhaseC"]):
    ctags = [t for (d, t) in out_map2[idx] if d == 3]
    assert len(ctags) > 0, f"Missing child volumes for {name}"
    gmsh.model.addPhysicalGroup(3, ctags, idx + 1, name)
    all_child_tags.extend(ctags)

gmsh.model.mesh.generate(3)
gmsh.model.mesh.optimize("", True)

vol_total, min_sicn = compute_tet_volume_and_min_sicn(all_child_tags)
assert math.isclose(vol_total, 8.0, rel_tol=1e-6), f"Expected 8.0, got {vol_total}"
assert min_sicn > 0.0, f"Non-positive SICN: {min_sicn}"

if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
