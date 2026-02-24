#Args:
shading = True
VR = False
debug = True
edges = True
edges_color = "#000000"

from scene import *
from math import sin, cos, tan, radians, sqrt, fabs
import motion, time
import numpy as np
import ui, console

# Shortcuts
_time = time.time
_get_att = motion.get_attitude
_cos = cos; _sin = sin; _sqrt = sqrt
_fill = fill; _abs = abs
_float = float
_int = int
_range = range
_min = min; _max = max
_fabs = fabs
_enumerate = enumerate
_triangle_strip = triangle_strip
eps = 1e-6

CUBE_EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
]


class Gizmo(ShapeNode):
    def __init__(self, position, dister):
        super().__init__()
        self.sizer = 10
        self.dister = dister
        self.z_max = self.z_min = 10
        self.thick = 4
        self.alpha = 0.5
        self.color = "#000"
        self.size = (5, 5)
        self.position = position
        h = self.sizer * 0.5
        verts = [(-h, -h, -h), ( h, -h, -h), ( h,  h, -h), (-h,  h, -h),
                 (-h, -h,  h), ( h, -h,  h), ( h,  h,  h), (-h,  h,  h)]
        path = ui.Path()
        for s, e in CUBE_EDGES:
            path.move_to(*verts[s][:2])
            path.line_to(*verts[e][:2])
        path.line_width = self.thick
        self.x_cube = ShapeNode(path, stroke_color="#FF0000", parent=self)
        self.y_cube = ShapeNode(path, stroke_color="#0000FF", parent=self)
        self.z_cube = ShapeNode(path, stroke_color="#008000", parent=self)

    def update_gizmo(self, roll, pitch, yaw):
        cr, sr = _cos(roll), _sin(roll)
        cp, sp = _cos(pitch), _sin(pitch)
        cy, sy_ = _cos(yaw), _sin(yaw)
        m00 = cy*cr + sy_*sp*sr;  m01 = -cy*sr + sy_*sp*cr;  m02 = sy_*cp
        m10 = cp*sr;              m11 = cp*cr;                m12 = -sp
        m20 = -sy_*cr + cy*sp*sr; m21 = sy_*sr + cy*sp*cr;   m22 = cy*cp
        cx, cy_ = self.size[0]*0.5, self.size[1]*0.5
        d = self.dister
        z_max = self.z_max
        z_min = self.z_min
        for cube, vx, vy, vz in (
            (self.x_cube, d, 0, 0),
            (self.y_cube, 0, d, 0),
            (self.z_cube, 0, 0, d),
        ):
            tx = m00*vx + m01*vy + m02*vz
            ty = m10*vx + m11*vy + m12*vz
            tz = m20*vx + m21*vy + m22*vz
            s = max(1.0, z_max - (tz / z_min)) * 0.1
            cube.position = (cx + tx, cy_ + ty)
            cube.size = (s*5, s*5)


def _parse_hex(col):
    """Pre-parse '#rrggbb' into (r, g, b) ints."""
    return (int(col[1:3], 16), int(col[3:5], 16), int(col[5:7], 16))


class Game(Scene):
    def setup(self):
        motion.start_updates()
        self.background_color = '#ffffff'
        self.pos = [50.0, 50.0, 50.0]
        self.vel = [0.0, 0.0, 0.0]
        self._t = _time()
        self.pitch = self.yaw = self.roll = 0.0
        self.fov_tan = tan(radians(70) * 0.5)
        self.Plus = radians(90)
        self.res = 10

        if debug:
            if VR:
                w, h = self.size.w, self.size.h
                self.gizmos = [
                    Gizmo(position=(w / 4, h / 2), dister=50),
                    Gizmo(position=(w * 3 / 4, h / 2), dister=50)
                ]
            else:
                self.gizmos = [Gizmo(position=(self.size.w / 2, self.size.h / 2), dister=50)]
            for g in self.gizmos:
                self.add_child(g)

        B = []
        raw_blocks = [
            (0, 0, 0, 10, '#ff0000'), (0, 15, 0, 10, '#00ff00'),
            (0, 30, 0, 10, '#0000ff'), (0, 45, 0, 10, '#ffff00'),
            (0, 60, 0, 10, '#ff00ff'), (15, 0, 0, 10, '#ff0000'),
            (15, 15, 0, 10, '#00ff00'), (15, 30, 0, 10, '#0000ff'),
            (15, 45, 0, 10, '#ffff00'), (15, 60, 0, 10, '#ff00ff'),
            (30, 0, 0, 10, '#ff0000'), (30, 15, 0, 10, '#00ff00'),
            (30, 30, 0, 10, '#0000ff'), (30, 45, 0, 10, '#ffff00'),
            (30, 60, 0, 10, '#ff00ff'), (45, 0, 0, 10, '#ff0000'),
            (45, 15, 0, 10, '#00ff00'), (45, 30, 0, 10, '#0000ff'),
            (45, 45, 0, 10, '#ffff00'), (45, 60, 0, 10, '#ff00ff'),
            (60, 0, 0, 10, '#ff0000'), (60, 15, 0, 10, '#00ff00'),
            (60, 30, 0, 10, '#0000ff'), (60, 45, 0, 10, '#ffff00'),
            (60, 60, 0, 10, '#ff00ff'), (0, 0, 15, 10, '#ff0000'),
            (0, 15, 15, 10, '#00ff00'), (0, 30, 15, 10, '#0000ff'),
            (0, 45, 15, 10, '#ffff00'), (0, 60, 15, 10, '#ff00ff'),
            (15, 0, 15, 10, '#ff0000'), (15, 15, 15, 10, '#00ff00'),
            (15, 30, 15, 10, '#0000ff'), (15, 45, 15, 10, '#ffff00'),
            (15, 60, 15, 10, '#ff00ff'), (30, 0, 15, 10, '#ff0000'),
            (30, 15, 15, 10, '#00ff00'), (30, 30, 15, 10, '#0000ff'),
            (30, 45, 15, 10, '#ffff00'), (30, 60, 15, 10, '#ff00ff'),
            (45, 0, 15, 10, '#ff0000'), (45, 15, 15, 10, '#00ff00'),
            (45, 30, 15, 10, '#0000ff'), (45, 45, 15, 10, '#ffff00'),
            (45, 60, 15, 10, '#ff00ff'), (60, 0, 15, 10, '#ff0000'),
            (60, 15, 15, 10, '#00ff00'), (60, 30, 15, 10, '#0000ff'),
            (60, 45, 15, 10, '#ffff00'), (60, 60, 15, 10, '#ff00ff'),
        ]

        # --- Optimization: pre-parse colors and pre-compute centers ---
        for x, y, z, s, c in raw_blocks:
            xm = x + s; ym = y + s; zm = z + s
            rgb = _parse_hex(c)
            B.append((
                _float(x), _float(y), _float(z),
                _float(xm), _float(ym), _float(zm),
                _float(x + s * 0.5), _float(y + s * 0.5), _float(z + s * 0.5),
                c, rgb[0], rgb[1], rgb[2]
            ))
        self.blocks = B

        # --- Optimization: pre-compute pixel NDC values once (resize recalcs) ---
        self._precompute_screen(self.size.w, self.size.h)

    def _precompute_screen(self, w, h):
        """Cache screen-space ray multipliers so we don't recompute per-frame."""
        res = self.res
        fov = self.fov_tan
        aspect = w / h
        inv_sx2 = 2.0 * res / w
        inv_sy2 = 2.0 * res / h

        if VR:
            sx = _int(w // res) // 2
        else:
            sx = _int(w // res)
        sy = _int(h // res)

        # Store as tuples for fast iteration (avoid repeated attr lookup)
        self._px_mul = tuple(((px + 0.5) * inv_sx2 - 1.0) * fov * aspect for px in _range(sx))
        self._py_mul = tuple(((py + 0.5) * inv_sy2 - 1.0) * fov for py in _range(sy))
        self._sx = sx
        self._sy = sy
        self._cached_w = w
        self._cached_h = h

    def motion_update(self):
        t = _time()
        self._t = t
        roll, pitch, yaw = _get_att()
        self.roll = -pitch
        self.pitch = roll + self.Plus
        self.yaw = -yaw

    def render(self):
        posx, posy, posz = self.pos
        w, h, res = self.size.w, self.size.h, self.res

        # Recompute screen cache only if size changed
        if w != self._cached_w or h != self._cached_h:
            self._precompute_screen(w, h)

        sx = self._sx
        sy = self._sy
        px_mul = self._px_mul
        py_mul = self._py_mul

        if VR:
            x_offsets = (0, w * 0.5)
        else:
            x_offsets = (0,)

        # --- Rotation matrix (computed once per frame) ---
        roll_ = self.roll; pitch_ = self.pitch; yaw_ = self.yaw
        cr, sr = _cos(roll_), _sin(roll_)
        cp, sp = _cos(pitch_), _sin(pitch_)
        cy, sy_ = _cos(yaw_), _sin(yaw_)
        r00 = cy*cr + sy_*sp*sr;  r01 = -cy*sr + sy_*sp*cr;  r02 = sy_*cp
        r10 = cp*sr;              r11 = cp*cr;                r12 = -sp
        r20 = -sy_*cr + cy*sp*sr; r21 = sy_*sr + cy*sp*cr;   r22 = cy*cp

        blocks = self.blocks
        nblocks = len(blocks)

        # --- Optimization: pre-sort blocks by distance to camera (front-to-back) ---
        # Allows early-out when a closer block is already found.
        block_dists = []
        for i in _range(nblocks):
            b = blocks[i]
            dx2 = b[6] - posx; dy2 = b[7] - posy; dz2 = b[8] - posz
            block_dists.append((dx2*dx2 + dy2*dy2 + dz2*dz2, i))
        block_dists.sort()
        sorted_indices = tuple(bd[1] for bd in block_dists)
        sorted_dists = tuple(bd[0] for bd in block_dists)

        # Light (normalized once)
        lx, ly, lz = 0.5773502691896258, 0.6928203230275509, -0.4618802153517006
        # ^ Pre-normalized (1, 1.2, -0.8) / ||(1, 1.2, -0.8)||
        MIN = 0.3
        do_shading = shading
        do_edges = edges

        for x_off in x_offsets:
            hit_buf   = [[False]*sx for _ in _range(sy)]
            color_buf = [[''] * sx for _ in _range(sy)]

            for pyi in _range(len(py_mul)):
                ny = py_mul[pyi]
                row_hit = hit_buf[pyi]
                row_col = color_buf[pyi]

                for pxi in _range(len(px_mul)):
                    nx = px_mul[pxi]

                    # --- Ray direction (rotated) ---
                    dx = r00*nx + r01*ny + r02
                    dy = r10*nx + r11*ny + r12
                    dz = r20*nx + r21*ny + r22
                    inv_len = 1.0 / _sqrt(dx*dx + dy*dy + dz*dz)
                    rx = dx*inv_len; ry = dy*inv_len; rz = dz*inv_len

                    # Ray origin (offset by one unit along ray)
                    x0 = posx + rx; y0 = posy + ry; z0 = posz + rz

                    # Inverse ray directions for slab test
                    arx = _fabs(rx); ary = _fabs(ry); arz = _fabs(rz)
                    inv_rx = 1.0/rx if arx > eps else 1e18
                    inv_ry = 1.0/ry if ary > eps else 1e18
                    inv_rz = 1.0/rz if arz > eps else 1e18

                    best_t = 1e18
                    best_idx = -1

                    # --- AABB ray intersection (front-to-back with early skip) ---
                    for si in _range(nblocks):
                        # If the closest possible distance for remaining blocks
                        # exceeds best hit, we can break early.
                        # (squared dist to center vs. squared ray-t is an approximation,
                        #  but sorted order still helps skip far blocks.)
                        idx = sorted_indices[si]
                        b = blocks[idx]
                        xmin = b[0]; ymin = b[1]; zmin = b[2]
                        xmax = b[3]; ymax = b[4]; zmax = b[5]

                        t1 = (xmin - x0) * inv_rx; t2 = (xmax - x0) * inv_rx
                        if t1 > t2: t1, t2 = t2, t1
                        u1 = (ymin - y0) * inv_ry; u2 = (ymax - y0) * inv_ry
                        if u1 > u2: u1, u2 = u2, u1
                        v1 = (zmin - z0) * inv_rz; v2 = (zmax - z0) * inv_rz
                        if v1 > v2: v1, v2 = v2, v1

                        t_near = t1 if t1 > u1 else u1
                        if v1 > t_near: t_near = v1
                        t_far = t2 if t2 < u2 else u2
                        if v2 < t_far: t_far = v2

                        if t_near <= t_far and t_far > 0.0 and t_near < best_t:
                            best_t = t_near
                            best_idx = idx

                    if best_idx < 0:
                        continue

                    b = blocks[best_idx]
                    c = b[9]  # hex color string

                    # --- Combined shading + edges (compute hit point once) ---
                    if do_shading or do_edges:
                        hx = x0 + rx * best_t
                        hy = y0 + ry * best_t
                        hz = z0 + rz * best_t
                        xmin = b[0]; ymin = b[1]; zmin = b[2]
                        xmax = b[3]; ymax = b[4]; zmax = b[5]
                        bs = _max(xmax-xmin, ymax-ymin, zmax-zmin)

                        if do_shading:
                            tol = bs * 1e-3
                            dhx_min = _fabs(hx - xmin); dhx_max = _fabs(hx - xmax)
                            dhy_min = _fabs(hy - ymin); dhy_max = _fabs(hy - ymax)
                            dhz_min = _fabs(hz - zmin); dhz_max = _fabs(hz - zmax)
                            if   dhx_min < tol: nl = -lx
                            elif dhx_max < tol: nl =  lx
                            elif dhy_min < tol: nl = -ly
                            elif dhy_max < tol: nl =  ly
                            elif dhz_min < tol: nl = -lz
                            elif dhz_max < tol: nl =  lz
                            else:               nl = 0.0
                            if nl < 0.0: nl = 0.0
                            intensity = MIN + (1.0 - MIN) * nl
                            # Pre-parsed RGB from blocks tuple
                            r0 = b[10]; g0 = b[11]; b0 = b[12]
                            ri = _int(r0 * intensity)
                            gi = _int(g0 * intensity)
                            bi = _int(b0 * intensity)
                            c = f"#{ri:02x}{gi:02x}{bi:02x}"

                        if do_edges:
                            tol_e = bs * 0.01
                            if not do_shading:
                                dhx_min = _fabs(hx - xmin); dhx_max = _fabs(hx - xmax)
                                dhy_min = _fabs(hy - ymin); dhy_max = _fabs(hy - ymax)
                                dhz_min = _fabs(hz - zmin); dhz_max = _fabs(hz - zmax)
                            cnt = (
                                (dhx_min < tol_e or dhx_max < tol_e) +
                                (dhy_min < tol_e or dhy_max < tol_e) +
                                (dhz_min < tol_e or dhz_max < tol_e)
                            )
                            if cnt >= 2:
                                c = edges_color

                    row_hit[pxi] = True
                    row_col[pxi] = c

            # --- Draw: run-length encoded triangle strips ---
            for pyi in _range(sy):
                row_hit = hit_buf[pyi]
                row_col = color_buf[pyi]
                y0_ = pyi * res
                y1_ = y0_ + res
                pxi = 0
                while pxi < sx:
                    if not row_hit[pxi]:
                        pxi += 1
                        continue
                    run_color = row_col[pxi]
                    start = pxi
                    pxi += 1
                    while pxi < sx and row_hit[pxi] and row_col[pxi] == run_color:
                        pxi += 1
                    x0_ = x_off + start * res
                    x1_ = x_off + pxi * res
                    _fill(run_color)
                    _triangle_strip([(x0_, y0_), (x1_, y0_), (x0_, y1_), (x1_, y1_)])

    def stop(self):
        motion.stop_updates()

    def update(self):
        self.motion_update()
        self.render()
        if debug:
            for g in self.gizmos:
                g.update_gizmo(roll=self.roll, pitch=self.pitch, yaw=self.yaw)

if __name__ == '__main__':
    run(Game(), show_fps=True)
