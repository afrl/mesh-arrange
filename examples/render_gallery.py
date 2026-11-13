#!/usr/bin/env python3
"""Render the MeshArrange README demonstration gallery from exported C++ demo bundles.

Populates 3D PyVista off-screen studio scenes (and runs Gmsh 3D tetrahedral meshing
on the discrete B-Rep from Demo 5) to generate publication-quality multi-panel figures
in `docs/images/`.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Auto-relaunch inside xvfb-run when running headless without an X11 display
if not os.environ.get("DISPLAY") and not os.environ.get("MESHARRANGE_UNDER_XVFB"):
    xvfb = shutil.which("xvfb-run")
    if xvfb:
        env = os.environ.copy()
        env["MESHARRANGE_UNDER_XVFB"] = "1"
        sys.exit(
            subprocess.call(
                [xvfb, "-a", "-s", "-screen 0 1600x1300x24", sys.executable, *sys.argv],
                env=env,
            )
        )

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import pyvista as pv
import vtk

pv.OFF_SCREEN = True

# Curated high-clarity scientific palette
PALETTE = {
    "coral": "#e76f51",
    "cobalt": "#3b82f6",
    "teal": "#2a9d8f",
    "amber": "#f4a261",
    "gold": "#eab308",
    "violet": "#8b5cf6",
    "crimson": "#dc2626",
    "emerald": "#10b981",
    "slate_light": "#cbd5e1",
    "resin_blue": "#93c5fd",
    "seam_dark": "#0f172a",
    "edge": "#334155",
}


def make_polydata(
    vertices: list | np.ndarray,
    triangles: list | np.ndarray,
    orient_outward_from: np.ndarray | None = None,
) -> pv.PolyData:
    pts = np.asarray(vertices, dtype=np.float64)
    tris = np.asarray(triangles, dtype=np.int64).copy()
    if len(tris) == 0:
        return pv.PolyData()
    if orient_outward_from is not None:
        v0 = pts[tris[:, 0]]
        v1 = pts[tris[:, 1]]
        v2 = pts[tris[:, 2]]
        nrm = np.cross(v1 - v0, v2 - v0)
        cent = (v0 + v1 + v2) / 3.0
        flip = np.einsum("ij,ij->i", nrm, cent - orient_outward_from) < 0.0
        tris[flip, 1], tris[flip, 2] = tris[flip, 2].copy(), tris[flip, 1].copy()

    faces = np.hstack([np.full((len(tris), 1), 3, dtype=np.int64), tris]).ravel()
    poly = pv.PolyData(pts, faces).clean()
    poly.compute_normals(
        cell_normals=True,
        point_normals=False,
        consistent_normals=True,
        auto_orient_normals=False,
        inplace=True,
    )
    return poly


def make_volume_shell_polydata(result_dict: dict, bvol: dict) -> pv.PolyData:
    verts = np.asarray(result_dict["surface_mesh"]["vertices"], dtype=np.float64)
    all_tris = []
    for fid, sign in bvol["signed_bounding_faces"]:
        ftris = np.asarray(result_dict["brep_faces"][fid]["triangles"], dtype=np.int64).copy()
        if sign < 0 and len(ftris) > 0:
            ftris[:, [1, 2]] = ftris[:, [2, 1]]
        if len(ftris) > 0:
            all_tris.append(ftris)
    if not all_tris:
        return pv.PolyData()
    tris = np.vstack(all_tris)
    faces = np.hstack([np.full((len(tris), 1), 3, dtype=np.int64), tris]).ravel()
    poly = pv.PolyData(verts, faces).clean()
    poly.compute_normals(
        cell_normals=True,
        point_normals=False,
        consistent_normals=True,
        auto_orient_normals=True,
        inplace=True,
    )
    return poly


def make_owner_shells(result_dict: dict) -> dict[int, pv.PolyData]:
    """Assemble smoothed closed shells per surviving solid owner_object from brep_volumes."""
    by_owner: dict[int, list[pv.PolyData]] = {}
    for bvol in result_dict["brep_volumes"]:
        if len(bvol["owner_objects"]) == 1:
            owner = int(bvol["owner_objects"][0])
            by_owner.setdefault(owner, []).append(make_volume_shell_polydata(result_dict, bvol))
    return {owner: pv.merge(parts) for owner, parts in by_owner.items()}


def make_studio_plotter(window_size: tuple[int, int] = (1450, 1250)) -> pv.Plotter:
    pl = pv.Plotter(off_screen=True, window_size=window_size, lighting="none")
    pl.set_background("white")
    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    try:
        pl.enable_depth_peeling(number_of_peels=8)
    except Exception:
        pass

    pl.add_light(pv.Light(position=(3.8, -3.6, 4.6), focal_point=(0, 0, 0), color="white", intensity=0.78))
    pl.add_light(pv.Light(position=(-4.0, -2.2, 2.8), focal_point=(0, 0, 0), color="#f1f5f9", intensity=0.45))
    pl.add_light(pv.Light(position=(0.2, 4.5, 3.4), focal_point=(0, 0, 0), color="#fffbeb", intensity=0.30))
    pl.add_light(pv.Light(light_type="headlight", intensity=0.18))
    return pl


def add_styled_mesh(
    pl: pv.Plotter,
    mesh: pv.DataSet,
    color: str,
    opacity: float = 1.0,
    show_edges: bool = True,
    edge_color: str = PALETTE["edge"],
    line_width: float = 0.80,
    specular: float = 0.16,
) -> None:
    if mesh.n_cells == 0:
        return
    pl.add_mesh(
        mesh,
        color=color,
        opacity=opacity,
        show_edges=show_edges,
        edge_color=edge_color,
        line_width=line_width,
        smooth_shading=False,
        ambient=0.25,
        diffuse=0.77,
        specular=specular,
        specular_power=24,
    )


def add_brep_seams_and_vertices(
    pl: pv.Plotter,
    result_dict: dict,
    tube_radius: float = 0.011,
    vertex_radius: float = 0.024,
    curve_color: str = PALETTE["seam_dark"],
    vertex_color: str = PALETTE["crimson"],
) -> None:
    verts = np.asarray(result_dict["surface_mesh"]["vertices"], dtype=np.float64)
    seg_list = []
    for bc in result_dict["brep_curves"]:
        poly = bc["polyline"]
        if len(poly) < 2:
            continue
        pts = verts[np.asarray(poly, dtype=np.int64)]
        n_pts = len(pts)
        lines = np.hstack([[n_pts], np.arange(n_pts, dtype=np.int64)])
        pdata = pv.PolyData(pts, lines=lines)
        seg_list.append(pdata.tube(radius=tube_radius, n_sides=14))
    if seg_list:
        pl.add_mesh(
            pv.merge(seg_list),
            color=curve_color,
            smooth_shading=True,
            ambient=0.25,
            diffuse=0.75,
            specular=0.40,
            specular_power=30,
        )

    if vertex_radius > 0 and result_dict["brep_vertices"]:
        spheres = [
            pv.Sphere(radius=vertex_radius, center=bv["point"], theta_resolution=18, phi_resolution=18)
            for bv in result_dict["brep_vertices"]
        ]
        pl.add_mesh(
            pv.merge(spheres),
            color=vertex_color,
            smooth_shading=True,
            ambient=0.25,
            diffuse=0.75,
            specular=0.50,
            specular_power=35,
        )


def set_iso_camera(pl: pv.Plotter, azimuth: float = 38.0, elevation: float = 28.0, zoom: float = 1.08) -> None:
    pl.camera_position = "yz"
    pl.camera.azimuth = azimuth
    pl.camera.elevation = elevation
    pl.reset_camera()
    pl.camera.zoom(zoom)


def capture_cropped(pl: pv.Plotter, pad_frac: float = 0.055) -> np.ndarray:
    img = pl.screenshot(return_img=True)
    pl.close()
    mask = np.any(img < 252, axis=2)
    if not np.any(mask):
        return img
    ys, xs = np.where(mask)
    y0, y1 = int(ys.min()), int(ys.max())
    x0, x1 = int(xs.min()), int(xs.max())
    h, w = y1 - y0 + 1, x1 - x0 + 1
    pad = int(max(h, w) * pad_frac)
    y0 = max(0, y0 - pad)
    y1 = min(img.shape[0], y1 + pad)
    x0 = max(0, x0 - pad)
    x1 = min(img.shape[1], x1 + pad)
    cropped = img[y0:y1, x0:x1]

    ch, cw = cropped.shape[:2]
    target_ratio = 1.15
    if cw / ch < target_ratio:
        new_w = int(ch * target_ratio)
        new_h = ch
    else:
        new_w = cw
        new_h = int(cw / target_ratio)
    canvas = np.full((new_h, new_w, 3), 255, dtype=np.uint8)
    oy = (new_h - ch) // 2
    ox = (new_w - cw) // 2
    canvas[oy : oy + ch, ox : ox + cw] = cropped
    return canvas


def save_three_panel_card(
    out_path: Path,
    main_title: str,
    subtitle: str,
    panels: list[tuple[np.ndarray, str, str]],
    legend_items: list[tuple[str, str]],
) -> None:
    fig = plt.figure(figsize=(16.8, 6.8), dpi=200, facecolor="#ffffff")
    gs = fig.add_gridspec(
        nrows=3,
        ncols=3,
        height_ratios=[0.17, 0.71, 0.12],
        wspace=0.045,
        hspace=0.14,
        left=0.022,
        right=0.978,
        top=0.965,
        bottom=0.025,
    )

    ax_hdr = fig.add_subplot(gs[0, :])
    ax_hdr.axis("off")
    ax_hdr.text(
        0.0,
        0.76,
        main_title,
        fontsize=16.2,
        fontweight="bold",
        color="#0f172a",
        va="center",
        ha="left",
        family="DejaVu Sans",
    )
    ax_hdr.text(
        0.0,
        0.24,
        subtitle,
        fontsize=10.8,
        color="#475569",
        va="center",
        ha="left",
        family="DejaVu Sans",
    )

    for idx, (img, ptitle, pcaption) in enumerate(panels):
        ax = fig.add_subplot(gs[1, idx])
        ax.imshow(img, interpolation="lanczos")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor("#cbd5e1")
            spine.set_linewidth(1.1)
        ax.set_title(ptitle, fontsize=11.8, fontweight="bold", color="#1e293b", pad=7, family="DejaVu Sans")
        ax.set_xlabel(pcaption, fontsize=9.8, color="#475569", labelpad=6, family="DejaVu Sans")

    ax_leg = fig.add_subplot(gs[2, :])
    ax_leg.axis("off")
    handles = [
        mpatches.Patch(facecolor=col, edgecolor="#334155", linewidth=0.8, label=lbl)
        for lbl, col in legend_items
    ]
    ax_leg.legend(
        handles=handles,
        loc="center",
        ncol=min(len(handles), 6),
        frameon=True,
        facecolor="#f8fafc",
        edgecolor="#e2e8f0",
        fontsize=10.0,
        handlelength=1.4,
        handleheight=0.9,
        columnspacing=1.5,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(".tmp.png")
    fig.savefig(tmp_path, dpi=200, facecolor="#ffffff")
    plt.close(fig)

    im = Image.open(tmp_path)
    if im.width > 2160:
        new_h = int(round(im.height * (2160 / im.width)))
        im = im.resize((2160, new_h), Image.Resampling.LANCZOS)
    im.save(out_path, optimize=True)
    tmp_path.unlink(missing_ok=True)
    print(f"Saved {out_path} ({im.size[0]}x{im.size[1]})")


# ==============================================================================
# Demo 1: Exact Non-Manifold Arrangement & Complete B-Rep Extraction
# ==============================================================================
def render_demo1(bundle_path: Path, out_dir: Path) -> np.ndarray:
    data = json.loads(bundle_path.read_text())
    inputs = data["inputs"]
    res = data["result"]

    # Panel A: Intersecting inputs
    pl = make_studio_plotter()
    add_styled_mesh(pl, make_polydata(inputs[0]["vertices"], inputs[0]["triangles"]), PALETTE["cobalt"], opacity=0.86)
    add_styled_mesh(pl, make_polydata(inputs[1]["vertices"], inputs[1]["triangles"]), PALETTE["coral"], opacity=0.82)
    set_iso_camera(pl, azimuth=-48, elevation=26, zoom=1.08)
    img_a = capture_cropped(pl)

    # Panel B: Conforming BRepFaces (extracted from auto-oriented closed BRepVolume shells
    # so every triangle normal is 100% consistent) + 1D BRepCurves + 0D BRepVertices
    pl = make_studio_plotter()
    face_colors = [PALETTE["cobalt"], PALETTE["teal"], PALETTE["amber"], PALETTE["violet"]]
    drawn_fids: set[int] = set()
    # First orient each closed BRepVolume shell outward via VTK, tracking face IDs per cell
    for bvol in res["brep_volumes"]:
        verts = np.asarray(res["surface_mesh"]["vertices"], dtype=np.float64)
        tris_list = []
        fid_list = []
        for fid, sign in bvol["signed_bounding_faces"]:
            ftris = np.asarray(res["brep_faces"][fid]["triangles"], dtype=np.int64).copy()
            if sign < 0 and len(ftris) > 0:
                ftris[:, [1, 2]] = ftris[:, [2, 1]]
            if len(ftris) > 0:
                tris_list.append(ftris)
                fid_list.append(np.full(len(ftris), fid, dtype=np.int32))
        tris_all = np.vstack(tris_list)
        fids_all = np.concatenate(fid_list)
        faces_vtk = np.hstack([np.full((len(tris_all), 1), 3, dtype=np.int64), tris_all]).ravel()
        shell = pv.PolyData(verts, faces_vtk).clean()
        shell.cell_data["fid"] = fids_all
        shell.compute_normals(
            cell_normals=True,
            point_normals=False,
            consistent_normals=True,
            auto_orient_normals=True,
            inplace=True,
        )
        for fid, _ in bvol["signed_bounding_faces"]:
            if fid in drawn_fids:
                continue
            drawn_fids.add(fid)
            sub = shell.extract_cells(shell.cell_data["fid"] == fid).extract_surface(algorithm="dataset_surface")
            bf = res["brep_faces"][fid]
            is_internal = (bf["bounded_regions"][0] >= 0) and (bf["bounded_regions"][1] >= 0)
            # Keep torus opaque and sphere semi-translucent so interior lens BRepFaces show through!
            op = 1.0 if (not is_internal and bf["origin_objects"] == [0]) else (0.58 if not is_internal else 1.0)
            add_styled_mesh(pl, sub, face_colors[fid % len(face_colors)], opacity=op)

    add_brep_seams_and_vertices(pl, res, tube_radius=0.015, vertex_radius=0.034)
    set_iso_camera(pl, azimuth=-48, elevation=26, zoom=1.08)
    img_b = capture_cropped(pl)

    # Panel C: Exploded 3D BRepVolume shells (Torus\Sphere, Overlap A∩B, Sphere\Torus)
    pl = make_studio_plotter()
    for bvol in res["brep_volumes"]:
        shell = make_volume_shell_polydata(res, bvol)
        owners = bvol["owner_objects"]
        if bvol["is_overlap"]:
            shift = np.array([0.0, 0.0, 0.0])
            col = PALETTE["amber"]
        elif owners == [0]:
            shift = np.array([-0.36, 0.0, 0.0])
            col = PALETTE["cobalt"]
        else:
            shift = np.array([+0.38, 0.0, 0.0])
            col = PALETTE["coral"]
        shell.points = shell.points + shift
        add_styled_mesh(pl, shell, col, opacity=1.0)
    set_iso_camera(pl, azimuth=-52, elevation=25, zoom=1.06)
    img_c = capture_cropped(pl)

    nbv = len(res["brep_volumes"])
    nbf = len(res["brep_faces"])
    nbc = len(res["brep_curves"])
    nbp = len(res["brep_vertices"])
    save_three_panel_card(
        out_dir / "demo1_brep_arrangement.png",
        "1. Exact Surface Arrangement & Non-Manifold B-Rep Topology",
        f"Exact CGAL::Epeck intersection of closed triangular meshes without overlap removal — extracting {nbv} BRepVolumes, {nbf} BRepFaces, {nbc} BRepCurves, and {nbp} BRepVertices.",
        [
            (img_a, "(a) Intersecting Input Surfaces", "Closed Torus & Icosphere with non-conforming overlap"),
            (img_b, "(b) Conforming BRepFaces & Seams", f"{nbf} manifold BRepFaces + {nbc} BRepCurve loops + {nbp} BRepVertices"),
            (img_c, "(c) Exploded 3D BRepVolume Shells", "Watertight shells: Torus \\ Sphere, Overlap (A ∩ B), Sphere \\ Torus"),
        ],
        [
            ("Torus \\ Sphere (BRepVolume)", PALETTE["cobalt"]),
            ("Overlap Lens A ∩ B (is_overlap=true)", PALETTE["amber"]),
            ("Sphere \\ Torus (BRepVolume)", PALETTE["coral"]),
            ("Exterior Sphere BRepFace", PALETTE["teal"]),
            ("1D BRepCurve Seam", PALETTE["seam_dark"]),
            ("0D BRepVertex Anchor", PALETTE["crimson"]),
        ],
    )
    return img_b


# ==============================================================================
# Demo 2: Priority Hierarchy & Boolean Domain Clipping
# ==============================================================================
def render_demo2(bundle_path: Path, out_dir: Path) -> None:
    data = json.loads(bundle_path.read_text())
    inputs = data["inputs"]
    res = data["result"]
    shells = make_owner_shells(res)

    # Panel A: Input scene showing priorities & clipper roles
    pl = make_studio_plotter()
    add_styled_mesh(pl, make_polydata(inputs[0]["vertices"], inputs[0]["triangles"]), PALETTE["coral"], opacity=0.95)
    add_styled_mesh(pl, make_polydata(inputs[1]["vertices"], inputs[1]["triangles"]), PALETTE["amber"], opacity=0.95)
    add_styled_mesh(pl, make_polydata(inputs[2]["vertices"], inputs[2]["triangles"]), PALETTE["resin_blue"], opacity=0.20, show_edges=False)
    add_styled_mesh(pl, make_polydata(inputs[3]["vertices"], inputs[3]["triangles"]), PALETTE["crimson"], opacity=0.45, edge_color="#991b1b")
    dom_poly = make_polydata(inputs[4]["vertices"], inputs[4]["triangles"])
    feat_edges = dom_poly.extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False)
    pl.add_mesh(feat_edges.tube(radius=0.012), color=PALETTE["seam_dark"], smooth_shading=True)
    set_iso_camera(pl, azimuth=38, elevation=28, zoom=1.06)
    img_a = capture_cropped(pl)

    # Panel B: Clipped multi-material RVE assembly with translucent matrix & seams
    pl = make_studio_plotter()
    for cid, poly in shells.items():
        if cid == 0:
            add_styled_mesh(pl, poly.copy(), PALETTE["coral"], opacity=1.0)
        elif cid == 1:
            add_styled_mesh(pl, poly.copy(), PALETTE["amber"], opacity=1.0)
        else:
            add_styled_mesh(pl, poly.copy(), PALETTE["resin_blue"], opacity=0.28, edge_color="#64748b", line_width=0.6)
    add_brep_seams_and_vertices(pl, res, tube_radius=0.009, vertex_radius=0.019)
    set_iso_camera(pl, azimuth=38, elevation=28, zoom=1.08)
    img_b = capture_cropped(pl)

    # Panel C: Exploded watertight solids showing carved cavities + central through-bore
    pl = make_studio_plotter()
    for cid, poly in shells.items():
        p = poly.copy()
        if cid == 0:
            p.points = p.points + np.array([0.0, 0.0, +0.46])
            add_styled_mesh(pl, p, PALETTE["coral"], opacity=1.0)
        elif cid == 1:
            p.points = p.points + np.array([0.0, 0.0, -0.46])
            add_styled_mesh(pl, p, PALETTE["amber"], opacity=1.0)
        else:
            add_styled_mesh(pl, p, PALETTE["slate_light"], opacity=1.0)
    set_iso_camera(pl, azimuth=38, elevation=26, zoom=1.06)
    img_c = capture_cropped(pl)

    save_three_panel_card(
        out_dir / "demo2_priority_clipping.png",
        "2. Priority Hierarchy & Boolean Clippers (SOLID, CLIP_INSIDE, CLIP_OUTSIDE)",
        "Single-pass RVE boolean assembly: high-priority inclusions carve out of a background matrix while a cylindrical bore and unit-cell box clip all solids.",
        [
            (img_a, "(a) Raw Inputs & Clipper Roles", "Oversized rods (p=+2), matrix (p=0), bore & unit-cell clippers (p=+5)"),
            (img_b, "(b) Clipped Conforming RVE Assembly", "Trimmed unit cell with translucent matrix, drilled rods & BRepSeams"),
            (img_c, "(c) Exploded Watertight Solids", "Drilled inclusions separated from sculpted matrix cavities"),
        ],
        [
            ("Inclusion X (SOLID, p=+2)", PALETTE["coral"]),
            ("Inclusion Y (SOLID, p=+2)", PALETTE["amber"]),
            ("Sculpted Matrix (SOLID, p=0)", PALETTE["slate_light"]),
            ("Through-Bore (CLIP_INSIDE, p=+5)", PALETTE["crimson"]),
            ("RVE Cutter / BRepSeams", PALETTE["seam_dark"]),
        ],
    )


# ==============================================================================
# Demo 3: Pairwise Harmonic Medial Interpenetration Removal (K = 2 Woven Tows)
# ==============================================================================
def render_demo3(bundle_path: Path, out_dir: Path) -> np.ndarray:
    data = json.loads(bundle_path.read_text())
    res = data["result"]
    before = data["before_result"]
    shells = make_owner_shells(res)

    # Panel A: Before — translucent crossing tows highlighting 3D overlap lens in crimson
    pl = make_studio_plotter()
    for bvol in before["brep_volumes"]:
        shell = make_volume_shell_polydata(before, bvol)
        if bvol["is_overlap"]:
            add_styled_mesh(pl, shell, PALETTE["crimson"], opacity=1.0, edge_color="#7f1d1d", line_width=0.95)
        elif bvol["owner_objects"] == [0]:
            add_styled_mesh(pl, shell, PALETTE["cobalt"], opacity=0.34, show_edges=True, edge_color="#64748b", line_width=0.55)
        else:
            add_styled_mesh(pl, shell, PALETTE["teal"], opacity=0.34, show_edges=True, edge_color="#64748b", line_width=0.55)
    set_iso_camera(pl, azimuth=36, elevation=27, zoom=1.10)
    img_a = capture_cropped(pl)

    # Panel B: After — disjoint smoothed harmonic partition with contact seam loop
    pl = make_studio_plotter()
    for cid, poly in shells.items():
        col = PALETTE["cobalt"] if cid == 0 else PALETTE["teal"]
        add_styled_mesh(pl, poly, col, opacity=1.0)
    add_brep_seams_and_vertices(pl, res, tube_radius=0.011, vertex_radius=0.024)
    set_iso_camera(pl, azimuth=36, elevation=27, zoom=1.10)
    img_b = capture_cropped(pl)

    # Panel C: Exploded view highlighting the smooth saddle harmonic medial sheet in amber!
    pl = make_studio_plotter()
    verts = np.asarray(res["surface_mesh"]["vertices"], dtype=np.float64)
    for bvol in res["brep_volumes"]:
        owner = bvol["owner_objects"][0]
        dz = +0.34 if owner == 0 else -0.34
        base_col = PALETTE["cobalt"] if owner == 0 else PALETTE["teal"]
        for fid, sign in bvol["signed_bounding_faces"]:
            bf = res["brep_faces"][fid]
            ftris = np.asarray(bf["triangles"], dtype=np.int64).copy()
            if sign < 0 and len(ftris) > 0:
                ftris[:, [1, 2]] = ftris[:, [2, 1]]
            poly = make_polydata(verts + np.array([0.0, 0.0, dz]), ftris)
            is_medial = (bf["bounded_regions"][0] >= 0) and (bf["bounded_regions"][1] >= 0)
            if is_medial and poly.n_cells > 0:
                poly = poly.smooth_taubin(n_iter=20, pass_band=0.08, boundary_smoothing=False)
                poly.compute_normals(cell_normals=True, point_normals=False, consistent_normals=True, inplace=True)
            add_styled_mesh(pl, poly, PALETTE["amber"] if is_medial else base_col, opacity=1.0, line_width=0.60)
    set_iso_camera(pl, azimuth=34, elevation=25, zoom=1.08)
    img_c = capture_cropped(pl)

    save_three_panel_card(
        out_dir / "demo3_harmonic_medial.png",
        "3. Harmonic Medial Interpenetration Removal on Curved Solids (K = 2)",
        "Equal-priority overlapping composite tows are partitioned along a clamped-cotangent harmonic zero-isosurface — producing a smooth, conforming saddle contact sheet.",
        [
            (img_a, "(a) Interpenetrating Woven Tows", "3D overlap lens (A ∩ B) highlighted inside crossover"),
            (img_b, "(b) Disjoint Harmonic Arrangement", "Zero overlap • 100% union volume conserved • smoothed seams"),
            (img_c, "(c) Exploded Saddle Medial Interface", "Matching doubly-curved harmonic contact patch on both solids"),
        ],
        [
            ("Warp Tow (SOLID, p=1)", PALETTE["cobalt"]),
            ("Weft Tow (SOLID, p=1)", PALETTE["teal"]),
            ("Initial Overlap Region (A ∩ B)", PALETTE["crimson"]),
            ("Harmonic Medial Contact Sheet", PALETTE["amber"]),
            ("1D BRepCurve Seam", PALETTE["seam_dark"]),
        ],
    )
    return img_c


# ==============================================================================
# Demo 4: Simultaneous Multi-Material Overlap Partitioning (K = 4 Cluster)
# ==============================================================================
def render_demo4(bundle_path: Path, out_dir: Path) -> np.ndarray:
    data = json.loads(bundle_path.read_text())
    inputs = data["inputs"]
    res = data["result"]
    before = data["before_result"]
    shells = make_owner_shells(res)

    sph_colors = [PALETTE["coral"], PALETTE["cobalt"], PALETTE["teal"], PALETTE["violet"]]

    # Panel A: 4 mutually overlapping spheres (translucent + K=4 core highlighted in gold)
    pl = make_studio_plotter()
    for idx, inp in enumerate(inputs):
        add_styled_mesh(pl, make_polydata(inp["vertices"], inp["triangles"]), sph_colors[idx], opacity=0.35, edge_color="#64748b", line_width=0.55)
    for bvol in before["brep_volumes"]:
        if len(bvol["owner_objects"]) == 4:
            core = make_volume_shell_polydata(before, bvol)
            add_styled_mesh(pl, core, PALETTE["gold"], opacity=1.0, edge_color="#854d0e", line_width=1.0)
    set_iso_camera(pl, azimuth=66, elevation=20, zoom=1.10)
    img_a = capture_cropped(pl)

    # Panel B: Conforming 4-way partitioned arrangement + triple-junction BRepCurves
    pl = make_studio_plotter()
    for cid, poly in shells.items():
        add_styled_mesh(pl, poly, sph_colors[cid % 4], opacity=1.0)
    add_brep_seams_and_vertices(pl, res, tube_radius=0.012, vertex_radius=0.026)
    set_iso_camera(pl, azimuth=66, elevation=20, zoom=1.10)
    img_b = capture_cropped(pl)

    # Panel C: Radially exploded 4 watertight cells viewed straight into the open tetrahedral gap
    # so all 4 wedges, their 6 smooth medial sheets, and the central gold quadruple vertex are visible!
    pl = make_studio_plotter()
    tet_dirs = [
        np.array([+1.0, +1.0, +1.0]),
        np.array([+1.0, -1.0, -1.0]),
        np.array([-1.0, +1.0, -1.0]),
        np.array([-1.0, -1.0, +1.0]),
    ]
    tet_dirs = [d / np.linalg.norm(d) for d in tet_dirs]
    verts = np.asarray(res["surface_mesh"]["vertices"], dtype=np.float64)
    for bvol in res["brep_volumes"]:
        cid = int(bvol["owner_objects"][0])
        disp = 0.34 * tet_dirs[cid]
        for fid, sign in bvol["signed_bounding_faces"]:
            bf = res["brep_faces"][fid]
            ftris = np.asarray(bf["triangles"], dtype=np.int64).copy()
            if sign < 0 and len(ftris) > 0:
                ftris[:, [1, 2]] = ftris[:, [2, 1]]
            poly = make_polydata(verts + disp, ftris)
            is_medial = (bf["bounded_regions"][0] >= 0) and (bf["bounded_regions"][1] >= 0)
            if is_medial and poly.n_cells > 0:
                poly = poly.smooth_taubin(n_iter=20, pass_band=0.08, boundary_smoothing=False)
                poly.compute_normals(cell_normals=True, point_normals=False, consistent_normals=True, inplace=True)
            add_styled_mesh(pl, poly, sph_colors[cid % 4], opacity=1.0, line_width=0.55)

    # Highlight the central quadruple BRepVertex at the origin!
    for bv in res["brep_vertices"]:
        pt = np.asarray(bv["point"])
        if np.linalg.norm(pt) < 0.15:
            pl.add_mesh(
                pv.Sphere(radius=0.055, center=pt, theta_resolution=24, phi_resolution=24),
                color=PALETTE["gold"],
                smooth_shading=True,
                specular=0.6,
            )
    # Viewing near the -Y axis projects the 4 tetrahedral directions (+,+,+), (+,-,-), (-,+,-), (-,-,+)
    # cleanly into the 4 screen quadrants around the central quadruple vertex!
    pl.camera_position = "xz"
    pl.camera.azimuth = 12
    pl.camera.elevation = 10
    pl.reset_camera()
    pl.camera.zoom(1.08)
    img_c = capture_cropped(pl)

    n_raw = len(before["brep_volumes"])
    save_three_panel_card(
        out_dir / "demo4_multiway_junction.png",
        "4. Simultaneous Multi-Material Junction Partitioning (K = 4 Sphere Cluster)",
        f"Coupled K-field simplex harmonic solver resolves {n_raw} overlapping subregions (including a central K=4 core) into 4 order-independent watertight cells.",
        [
            (img_a, "(a) 4 Overlapping Spheres (K = 4 Core)", f"{n_raw} raw subregions; central 4-way overlap core in gold"),
            (img_b, "(b) Coupled 4-Way Harmonic Partition", "Conforming triple-sheet BRepCurves & zero residual overlap"),
            (img_c, "(c) Exploded Spherical-Voronoi Cells", "6 smooth medial sheets meeting at central quadruple BRepVertex"),
        ],
        [
            ("Sphere 0 Cell", PALETTE["coral"]),
            ("Sphere 1 Cell", PALETTE["cobalt"]),
            ("Sphere 2 Cell", PALETTE["teal"]),
            ("Sphere 3 Cell", PALETTE["violet"]),
            ("K=4 Core / Quadruple Vertex", PALETTE["gold"]),
            ("Triple-Sheet BRepCurve Seam", PALETTE["seam_dark"]),
        ],
    )
    return img_c


# ==============================================================================
# Demo 5: End-to-End Gmsh Discrete GModel Population & 3D Conforming Tet Meshing
# ==============================================================================
def populate_gmsh_and_mesh_3d(
    res: dict, out_vtu_path: Path | None = None
) -> dict[int, pv.UnstructuredGrid]:
    """Populate a Gmsh discrete GModel from `ArrangeResult` and generate a conforming 3D tet mesh."""
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("mesh_arrange_rve")

    coords = np.asarray(res["surface_mesh"]["vertices"], dtype=np.float64)
    assigned_vids: set[int] = set()

    # 1. 0D discreteVertex entities
    for i, bv in enumerate(res["brep_vertices"]):
        vtag = i + 1
        vid = int(bv["global_vid"])
        gmsh.model.addDiscreteEntity(0, vtag)
        gmsh.model.mesh.addNodes(0, vtag, [vid + 1], coords[vid].tolist())
        gmsh.model.mesh.addElementsByType(vtag, 15, [], [vid + 1])
        assigned_vids.add(vid)

    # 2. 1D discreteEdge entities
    for i, bc in enumerate(res["brep_curves"]):
        etag = i + 1
        bnd = []
        if bc["start_bv"] >= 0:
            bnd.append(bc["start_bv"] + 1)
        if bc["end_bv"] >= 0:
            bnd.append(bc["end_bv"] + 1)
        gmsh.model.addDiscreteEntity(1, etag, bnd)

        poly = [int(v) for v in bc["polyline"]]
        new_nodes = [v for v in poly[1:-1] if v not in assigned_vids]
        if new_nodes:
            gmsh.model.mesh.addNodes(1, etag, [v + 1 for v in new_nodes], coords[new_nodes].ravel().tolist())
            assigned_vids.update(new_nodes)

        if len(poly) >= 2:
            line_nodes = []
            for k in range(len(poly) - 1):
                line_nodes.extend([poly[k] + 1, poly[k + 1] + 1])
            gmsh.model.mesh.addElementsByType(etag, 1, [], line_nodes)

    # 3. 2D discreteFace entities
    for i, bf in enumerate(res["brep_faces"]):
        ftag = i + 1
        bnd = [cid + 1 for cid, _ in bf["signed_bounding_curves"]]
        gmsh.model.addDiscreteEntity(2, ftag, bnd)

        ftris = np.asarray(bf["triangles"], dtype=np.int64)
        unique_v = [int(v) for v in np.unique(ftris) if int(v) not in assigned_vids]
        if unique_v:
            gmsh.model.mesh.addNodes(2, ftag, [v + 1 for v in unique_v], coords[unique_v].ravel().tolist())
            assigned_vids.update(unique_v)

        if len(ftris) > 0:
            gmsh.model.mesh.addElementsByType(ftag, 2, [], (ftris + 1).ravel().tolist())

        if bf["embedded_curves"]:
            gmsh.model.mesh.embed(1, [cid + 1 for cid in bf["embedded_curves"]], 2, ftag)
        if bf["embedded_vertices"]:
            gmsh.model.mesh.embed(0, [vid + 1 for vid in bf["embedded_vertices"]], 2, ftag)

    # 4. 3D volumes bounded by the discrete BRepFace shells
    for i, bvol in enumerate(res["brep_volumes"]):
        sl = gmsh.model.geo.addSurfaceLoop([fid + 1 for fid, _ in bvol["signed_bounding_faces"]])
        gmsh.model.geo.addVolume([sl], i + 1)
    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.11)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # 3D Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(3)

    ntags, all_coords, _ = gmsh.model.mesh.getNodes()
    coords_arr = np.asarray(all_coords, dtype=np.float64).reshape(-1, 3)
    nid_to_idx = {int(nt): idx for idx, nt in enumerate(ntags)}

    owner_grids: dict[int, pv.UnstructuredGrid] = {}
    owner_blocks: dict[int, list[np.ndarray]] = {}
    all_vol_blocks: list[np.ndarray] = []
    all_vol_ids: list[np.ndarray] = []
    all_owner_ids: list[np.ndarray] = []

    for i, bvol in enumerate(res["brep_volumes"]):
        rtag = i + 1
        owner = int(bvol["owner_objects"][0])
        etypes, etags, enodes = gmsh.model.mesh.getElements(3, rtag)
        for et, tags, nodes in zip(etypes, etags, enodes):
            if et == 4 and len(tags) > 0:
                blk = np.asarray(nodes, dtype=np.int64).reshape(-1, 4)
                owner_blocks.setdefault(owner, []).append(blk)
                all_vol_blocks.append(blk)
                all_vol_ids.append(np.full(len(blk), i, dtype=np.int32))
                all_owner_ids.append(np.full(len(blk), owner, dtype=np.int32))

    if out_vtu_path is not None and all_vol_blocks:
        all_conn = np.vstack(all_vol_blocks)
        unids = np.unique(all_conn)
        pts = np.array([coords_arr[nid_to_idx[int(n)]] for n in unids], dtype=np.float64)
        lut = np.zeros(int(np.max(unids)) + 1, dtype=np.int64)
        lut[unids] = np.arange(len(unids), dtype=np.int64)
        mconn = lut[all_conn].copy()

        # Ensure positive VTK_TETRA Jacobian across all cells
        p0 = pts[mconn[:, 0]]
        u = pts[mconn[:, 1]] - p0
        v = pts[mconn[:, 2]] - p0
        w = pts[mconn[:, 3]] - p0
        det = np.einsum("ij,ij->i", np.cross(u, v), w)
        neg = det < 0.0
        if np.any(neg):
            mconn[neg, 1], mconn[neg, 2] = mconn[neg, 2].copy(), mconn[neg, 1].copy()

        cells = np.hstack([np.full((len(mconn), 1), 4, dtype=np.int64), mconn]).ravel()
        ctypes = np.full(len(mconn), vtk.VTK_TETRA, dtype=np.uint8)
        full_grid = pv.UnstructuredGrid(cells, ctypes, pts)
        vol_arr = np.concatenate(all_vol_ids)
        own_arr = np.concatenate(all_owner_ids)
        full_grid.cell_data["BRepVolumeID"] = vol_arr
        full_grid.cell_data["OwnerObjectID"] = own_arr
        full_grid.cell_data["IsResin"] = (own_arr == 2).astype(np.int32)
        out_vtu_path.parent.mkdir(parents=True, exist_ok=True)
        full_grid.save(str(out_vtu_path), binary=True)
        print(f"Saved {out_vtu_path} ({full_grid.n_points:,} nodes, {full_grid.n_cells:,} tets)")

    for owner, blocks in owner_blocks.items():
        conn = np.vstack(blocks)
        unids = np.unique(conn)
        pts = np.array([coords_arr[nid_to_idx[int(n)]] for n in unids], dtype=np.float64)
        lut = np.zeros(int(np.max(unids)) + 1, dtype=np.int64)
        lut[unids] = np.arange(len(unids), dtype=np.int64)
        mconn = lut[conn]
        cells = np.hstack([np.full((len(mconn), 1), 4, dtype=np.int64), mconn]).ravel()
        ctypes = np.full(len(mconn), vtk.VTK_TETRA, dtype=np.uint8)
        owner_grids[owner] = pv.UnstructuredGrid(cells, ctypes, pts)

    gmsh.finalize()
    return owner_grids


def render_demo5(
    bundle_path: Path, out_dir: Path, meshes_dir: Path | None = None
) -> tuple[np.ndarray, np.ndarray]:
    data = json.loads(bundle_path.read_text())
    inputs = data["inputs"]
    res = data["result"]
    shells = make_owner_shells(res)
    gmsh_vtu = (meshes_dir / "demo5_gmsh_rve_gmsh3d.vtu") if meshes_dir is not None else None
    owner_grids = populate_gmsh_and_mesh_3d(res, out_vtu_path=gmsh_vtu)

    total_tets = sum(g.n_cells for g in owner_grids.values())
    print(f"Demo 5 Gmsh 3D total tetrahedra: {total_tets:,}")

    # Panel A: Input overlapping tows + RVE unit-cell frame
    pl = make_studio_plotter()
    add_styled_mesh(pl, make_polydata(inputs[0]["vertices"], inputs[0]["triangles"]), PALETTE["coral"], opacity=0.92)
    add_styled_mesh(pl, make_polydata(inputs[1]["vertices"], inputs[1]["triangles"]), PALETTE["cobalt"], opacity=0.92)
    add_styled_mesh(pl, make_polydata(inputs[2]["vertices"], inputs[2]["triangles"]), PALETTE["resin_blue"], opacity=0.16, show_edges=False)
    rve_poly = make_polydata(inputs[3]["vertices"], inputs[3]["triangles"])
    rve_edges = rve_poly.extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False)
    pl.add_mesh(rve_edges.tube(radius=0.011), color=PALETTE["seam_dark"], smooth_shading=True)
    set_iso_camera(pl, azimuth=38, elevation=28, zoom=1.08)
    img_a = capture_cropped(pl)

    # Panel B: Populated Gmsh discrete B-Rep (GFace patches + GEdge seams + GVertex corners)
    pl = make_studio_plotter()
    for cid, poly in shells.items():
        if cid == 0:
            add_styled_mesh(pl, poly, PALETTE["coral"], opacity=1.0)
        elif cid == 1:
            add_styled_mesh(pl, poly, PALETTE["cobalt"], opacity=1.0)
        else:
            add_styled_mesh(pl, poly, PALETTE["resin_blue"], opacity=0.26, edge_color="#64748b", line_width=0.55)
    add_brep_seams_and_vertices(pl, res, tube_radius=0.010, vertex_radius=0.021)
    set_iso_camera(pl, azimuth=38, elevation=28, zoom=1.10)
    img_b = capture_cropped(pl)

    # Panel C: Gmsh 3D Conforming Tetrahedral Cutaway (front-diagonal crinkle clip exposing crossover & resin tets!)
    pl = make_studio_plotter()
    # Normal pointing toward camera (+X, -Y) with invert=True keeps the back half facing the camera!
    clip_normal = (0.707, -0.707, 0.0)
    clip_origin = (0.02, -0.02, 0.0)
    for owner in sorted(owner_grids.keys()):
        clipped = owner_grids[owner].clip(normal=clip_normal, origin=clip_origin, invert=True, crinkle=True)
        if clipped.n_cells == 0:
            continue
        surf = clipped.extract_surface(algorithm="dataset_surface")
        surf.compute_normals(cell_normals=True, point_normals=False, consistent_normals=True, auto_orient_normals=True, inplace=True)
        col = PALETTE["coral"] if owner == 0 else (PALETTE["cobalt"] if owner == 1 else PALETTE["slate_light"])
        add_styled_mesh(pl, surf, col, opacity=1.0, edge_color="#1e293b", line_width=0.65)
    set_iso_camera(pl, azimuth=-42, elevation=26, zoom=1.10)
    img_c = capture_cropped(pl)

    nbv = len(res["brep_volumes"])
    nbf = len(res["brep_faces"])
    nbc = len(res["brep_curves"])
    save_three_panel_card(
        out_dir / "demo5_gmsh_volume_mesh.png",
        "5. End-to-End Gmsh Discrete GModel & Conforming 3D Tetrahedral Meshing",
        f"Discrete B-Rep ({nbv} GRegions, {nbf} GFaces, {nbc} GEdges) populated from ArrangeResult and meshed directly in Gmsh 3D ({total_tets:,} conforming tetrahedra).",
        [
            (img_a, "(a) Overlapping Tows & RVE Cutter", "Interpenetrating warp/weft tows + matrix + CLIP_OUTSIDE box"),
            (img_b, "(b) Discrete B-Rep (GModel Entities)", f"{nbf} oriented GFaces + {nbc} GEdge seams + {len(res['brep_vertices'])} GVertices"),
            (img_c, "(c) Gmsh 3D Tetrahedral Cutaway", f"Crinkle cutaway of {total_tets:,} conforming multi-material tetrahedra"),
        ],
        [
            ("Warp Tow 3D Tets (GRegion)", PALETTE["coral"]),
            ("Weft Tow 3D Tets (GRegion)", PALETTE["cobalt"]),
            ("Resin Matrix 3D Tets (GRegion)", PALETTE["slate_light"]),
            ("1D Discrete GEdge Seams", PALETTE["seam_dark"]),
            ("0D Discrete GVertex Corners", PALETTE["crimson"]),
        ],
    )
    return img_a, img_c


# ==============================================================================
# Hero Banner: 4-Stage Visual Overview at the Top of README.md
# ==============================================================================
def render_hero_banner(
    out_dir: Path,
    img_input_rve: np.ndarray,
    img_brep: np.ndarray,
    img_multiway: np.ndarray,
    img_gmsh_tet: np.ndarray,
) -> None:
    fig = plt.figure(figsize=(17.0, 5.15), dpi=200, facecolor="#ffffff")
    gs = fig.add_gridspec(
        nrows=2,
        ncols=4,
        height_ratios=[0.16, 0.84],
        wspace=0.04,
        hspace=0.14,
        left=0.018,
        right=0.982,
        top=0.955,
        bottom=0.075,
    )

    ax_hdr = fig.add_subplot(gs[0, :])
    ax_hdr.axis("off")
    ax_hdr.text(
        0.5,
        0.55,
        "MeshArrange — Exact 3D Surface Arrangements, Harmonic Overlap Removal & Gmsh B-Rep Extraction",
        fontsize=14.8,
        fontweight="bold",
        color="#0f172a",
        va="center",
        ha="center",
        family="DejaVu Sans",
    )

    tiles = [
        (img_input_rve, "1. Intersecting & Clipper Inputs", "Arbitrary solids + boolean priority rules"),
        (img_brep, "2. Non-Manifold B-Rep Topology", "Exact GVertex / GEdge / GFace / GRegion"),
        (img_multiway, "3. Harmonic Medial Partitioning", "Smooth K >= 2 overlap removal"),
        (img_gmsh_tet, "4. Conforming 3D Volume Meshing", "Zero-gap multi-material Gmsh tets"),
    ]
    for idx, (img, title, sub) in enumerate(tiles):
        ax = fig.add_subplot(gs[1, idx])
        ax.imshow(img, interpolation="lanczos")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor("#cbd5e1")
            spine.set_linewidth(1.1)
        ax.set_title(title, fontsize=11.2, fontweight="bold", color="#1e293b", pad=7, family="DejaVu Sans")
        ax.set_xlabel(sub, fontsize=9.4, color="#475569", labelpad=5, family="DejaVu Sans")

    out_path = out_dir / "hero_banner.png"
    tmp_path = out_path.with_suffix(".tmp.png")
    fig.savefig(tmp_path, dpi=200, facecolor="#ffffff")
    plt.close(fig)
    im = Image.open(tmp_path)
    if im.width > 2160:
        new_h = int(round(im.height * (2160 / im.width)))
        im = im.resize((2160, new_h), Image.Resampling.LANCZOS)
    im.save(out_path, optimize=True)
    tmp_path.unlink(missing_ok=True)
    print(f"Saved {out_path} ({im.size[0]}x{im.size[1]})")


def main() -> None:
    if len(sys.argv) < 3:
        print(
            "Usage: render_gallery.py <json_bundle_dir> <output_images_dir> [output_meshes_dir]",
            file=sys.stderr,
        )
        sys.exit(1)

    json_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    meshes_dir = Path(sys.argv[3]) if len(sys.argv) >= 4 else None
    out_dir.mkdir(parents=True, exist_ok=True)

    img_brep = render_demo1(json_dir / "demo1.json", out_dir)
    render_demo2(json_dir / "demo2.json", out_dir)
    render_demo3(json_dir / "demo3.json", out_dir)
    img_multiway = render_demo4(json_dir / "demo4.json", out_dir)
    img_rve_in, img_gmsh_tet = render_demo5(json_dir / "demo5.json", out_dir, meshes_dir=meshes_dir)
    render_hero_banner(out_dir, img_rve_in, img_brep, img_multiway, img_gmsh_tet)

    if meshes_dir is not None and meshes_dir.is_dir():
        for vtu_file in sorted(meshes_dir.glob("*.vtu")):
            grid = pv.read(str(vtu_file))
            grid.save(str(vtu_file), binary=True)


if __name__ == "__main__":
    main()
