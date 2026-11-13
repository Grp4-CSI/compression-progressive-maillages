"""Visualiseur OpenGL d'un flux .pmc, niveau par niveau.

    python tools/view_pmc.py out/fandisk.pmc

Touches : flèche droite / gauche = niveau suivant / précédent, Début / Fin = premier /
dernier niveau, W = arêtes, C = surbrillance des faces du niveau, R = vue initiale,
Échap = quitter. Souris : clic gauche pour tourner, clic droit pour déplacer, molette
pour zoomer. Nécessite glfw et PyOpenGL (pip install glfw PyOpenGL).
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import bitstream  # noqa: E402
from decoder import Decoder  # noqa: E402

BASE_COLOR = (0.55, 0.72, 0.85)
NEW_COLOR = (1.0, 0.75, 0.15)


@dataclass
class Snapshot:
    level: int
    vertices: np.ndarray        # (n, 3) positions déquantifiées
    faces: np.ndarray           # (m, 3)
    new: np.ndarray             # (m,) vrai pour les faces créées ou modifiées à ce niveau
    bits: int                   # bits du flux reçus pour afficher ce niveau


def build_snapshots(data: bytes) -> list[Snapshot]:
    """Décode le flux et garde le maillage obtenu après M0 puis après chaque niveau."""
    pm = bitstream.decode(data)
    _, stats = bitstream.encode(pm)             # encodage déterministe : mêmes tailles de blocs
    bits = stats.cumulative_bits()
    decoder = Decoder(pm)
    touched: set = set()
    snaps = []
    while True:
        fids = list(decoder.mesh.faces)
        faces = np.array([decoder.mesh.faces[f] for f in fids], dtype=np.int64).reshape(-1, 3)
        new = np.array([f in touched for f in fids], dtype=bool)
        if len(pm.fixed_faces):
            faces = np.vstack([faces, pm.fixed_faces])
            new = np.concatenate([new, np.zeros(len(pm.fixed_faces), dtype=bool)])
        pos = pm.quantizer.dequantize(decoder.positions[: decoder.count])
        snaps.append(Snapshot(decoder.level, pos, faces, new, bits[decoder.level]))
        if decoder.done:
            return snaps
        touched = {fid for ref in decoder.step() for fid, _, _ in ref.faces}


def triangle_arrays(snap: Snapshot, highlight: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tableaux « une ligne par sommet de triangle » (positions, normales, couleurs) pour OpenGL."""
    tri = snap.vertices[snap.faces]                                   # (m, 3, 3)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    colors = np.tile(BASE_COLOR, (len(tri), 1))
    if highlight:
        colors[snap.new] = NEW_COLOR
    rep = lambda a: np.repeat(a, 3, axis=0).astype(np.float32)
    return tri.reshape(-1, 3).astype(np.float32), rep(n), rep(colors)


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    try:
        import glfw
        from OpenGL import GL as gl, GLU as glu
    except ImportError:
        sys.exit("glfw et PyOpenGL sont nécessaires : pip install glfw PyOpenGL")

    with open(sys.argv[1], "rb") as f:
        snaps = build_snapshots(f.read())
    n_total = len(snaps[-1].vertices)
    all_pts = snaps[-1].vertices
    center = (all_pts.min(axis=0) + all_pts.max(axis=0)) / 2
    radius = float(np.linalg.norm(all_pts - center, axis=1).max()) or 1.0
    state = {"i": 0, "wire": True, "hl": True, "ax": 20.0, "ay": -30.0, "zoom": 3.0 * radius,
             "tx": 0.0, "ty": 0.0, "left": False, "right": False, "x": 0.0, "y": 0.0, "arrays": None}

    if not glfw.init():
        sys.exit("impossible d'initialiser GLFW")
    window = glfw.create_window(1200, 850, "pmc", None, None)
    glfw.make_context_current(window)
    glfw.swap_interval(1)

    def refresh():
        s = snaps[state["i"]]
        state["arrays"] = triangle_arrays(s, state["hl"])
        glfw.set_window_title(window, f"{os.path.basename(sys.argv[1])} — niveau {s.level}/{len(snaps) - 1} — "
                                      f"{len(s.vertices)} sommets, {len(s.faces)} faces — "
                                      f"{s.bits / 8:.0f} octets reçus ({s.bits / n_total:.1f} bits/sommet)")

    def on_key(win, key, scancode, action, mods):
        if action not in (glfw.PRESS, glfw.REPEAT):
            return
        i = state["i"]
        if key == glfw.KEY_ESCAPE:
            glfw.set_window_should_close(win, True)
        elif key == glfw.KEY_RIGHT:
            state["i"] = min(i + 1, len(snaps) - 1)
        elif key == glfw.KEY_LEFT:
            state["i"] = max(i - 1, 0)
        elif key == glfw.KEY_HOME:
            state["i"] = 0
        elif key == glfw.KEY_END:
            state["i"] = len(snaps) - 1
        elif key == glfw.KEY_W:
            state["wire"] = not state["wire"]
        elif key == glfw.KEY_C:
            state["hl"] = not state["hl"]
        elif key == glfw.KEY_R:
            state.update(ax=20.0, ay=-30.0, zoom=3.0 * radius, tx=0.0, ty=0.0)
        refresh()

    def on_button(win, button, action, mods):
        pressed = action != glfw.RELEASE
        if button == glfw.MOUSE_BUTTON_LEFT:
            state["left"] = pressed
        elif button == glfw.MOUSE_BUTTON_RIGHT:
            state["right"] = pressed
        state["x"], state["y"] = glfw.get_cursor_pos(win)

    def on_cursor(win, x, y):
        dx, dy = x - state["x"], y - state["y"]
        state["x"], state["y"] = x, y
        if state["left"]:
            state["ay"] += dx * 0.4
            state["ax"] += dy * 0.4
        elif state["right"]:
            f = state["zoom"] * 0.0015
            state["tx"] += dx * f
            state["ty"] -= dy * f

    def on_scroll(win, dx, dy):
        state["zoom"] = min(max(state["zoom"] * 0.9 ** dy, radius * 0.2), radius * 50)

    glfw.set_key_callback(window, on_key)
    glfw.set_mouse_button_callback(window, on_button)
    glfw.set_cursor_pos_callback(window, on_cursor)
    glfw.set_scroll_callback(window, on_scroll)
    refresh()

    gl.glEnable(gl.GL_DEPTH_TEST)
    gl.glEnable(gl.GL_LIGHTING)
    gl.glEnable(gl.GL_LIGHT0)
    gl.glEnable(gl.GL_COLOR_MATERIAL)
    gl.glLightModeli(gl.GL_LIGHT_MODEL_TWO_SIDE, gl.GL_TRUE)
    gl.glEnable(gl.GL_POLYGON_OFFSET_FILL)
    gl.glPolygonOffset(1.0, 1.0)

    while not glfw.window_should_close(window):
        glfw.poll_events()
        w, h = glfw.get_framebuffer_size(window)
        if w == 0 or h == 0:
            continue
        gl.glViewport(0, 0, w, h)
        gl.glClearColor(0.12, 0.12, 0.15, 1)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        glu.gluPerspective(45, w / h, radius * 0.01, radius * 200)
        gl.glMatrixMode(gl.GL_MODELVIEW)
        gl.glLoadIdentity()
        gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, (0.0, 0.0, 1.0, 0.0))
        gl.glTranslatef(state["tx"], state["ty"], -state["zoom"])
        gl.glRotatef(state["ax"], 1, 0, 0)
        gl.glRotatef(state["ay"], 0, 1, 0)
        gl.glTranslatef(*(-center))

        pts, normals, colors = state["arrays"]
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
        gl.glEnableClientState(gl.GL_COLOR_ARRAY)
        gl.glVertexPointer(3, gl.GL_FLOAT, 0, pts)
        gl.glNormalPointer(gl.GL_FLOAT, 0, normals)
        gl.glColorPointer(3, gl.GL_FLOAT, 0, colors)
        gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
        gl.glDrawArrays(gl.GL_TRIANGLES, 0, len(pts))
        gl.glDisableClientState(gl.GL_COLOR_ARRAY)
        gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
        if state["wire"]:
            gl.glDisable(gl.GL_LIGHTING)
            gl.glColor3f(0.08, 0.08, 0.08)
            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
            gl.glDrawArrays(gl.GL_TRIANGLES, 0, len(pts))
            gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
            gl.glEnable(gl.GL_LIGHTING)
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        glfw.swap_buffers(window)
    glfw.terminate()


if __name__ == "__main__":
    main()
