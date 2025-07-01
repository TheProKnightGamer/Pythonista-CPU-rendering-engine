# Toggle this to enable stereo-VR duplication
shading = True # ~0.1–.3 fps lost
VR = False   # ~0.1–0.3 fps lost
debug = True     # ~0.1–0.4 fps lost
edges = True #~ 0.1-0.5 fps lost
edges_color = "#000000"
from scene import *
from math import sin, cos, tan, radians, sqrt, fabs
import motion, time
import numpy as np
import ui, console

# Shortcuts
_time = time.time
_get_att = motion.get_attitude
_get_acc = motion.get_user_acceleration
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
        cy, sy_ = _cos(yaw),   _sin(yaw)
        m = [
            [ cy*cr + sy_*sp*sr,  -cy*sr + sy_*sp*cr,  sy_*cp ],
            [          cp*sr,           cp*cr,       -sp   ],
            [ -sy_*cr + cy*sp*sr,  sy_*sr + cy*sp*cr,  cy*cp ]
        ]
        cx, cy_ = self.size[0]*0.5, self.size[1]*0.5
        d = self.dister
        for cube, v in ((self.x_cube, (d,0,0)),
                    (self.y_cube, (0,d,0)),
                    (self.z_cube, (0,0,d))):
            tx = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2]
            ty = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2]
            tz = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]
            s = max(1.0, self.z_max - (tz/self.z_min)) * 0.1
            cube.position = (cx + tx, cy_ + ty)
            cube.size = (s*5, s*5)

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
                w = self.size.w
                h = self.size.h
                self.gizmos = [
                    Gizmo(position=(w / 4, h / 2), dister=50),
                    Gizmo(position=(w * 3 / 4, h / 2), dister=50)
                ]
            else:
                self.gizmos = [Gizmo(position=(self.size.w / 2, self.size.h / 2), dister=50)]
            for g in self.gizmos:
                self.add_child(g)

        B = []
        #X, Y, Z, Size, Color.
        for x, y, z, s, c in [
 (0, 0, 0, 10, '#ff0000'),
 (0, 15, 0, 10, '#00ff00'),
 (0, 30, 0, 10, '#0000ff'),
 (0, 45, 0, 10, '#ffff00'),
 (0, 60, 0, 10, '#ff00ff'),
 (15, 0, 0, 10, '#ff0000'),
 (15, 15, 0, 10, '#00ff00'),
 (15, 30, 0, 10, '#0000ff'),
 (15, 45, 0, 10, '#ffff00'),
 (15, 60, 0, 10, '#ff00ff'),
 (30, 0, 0, 10, '#ff0000'),
 (30, 15, 0, 10, '#00ff00'),
 (30, 30, 0, 10, '#0000ff'),
 (30, 45, 0, 10, '#ffff00'),
 (30, 60, 0, 10, '#ff00ff'),
 (45, 0, 0, 10, '#ff0000'),
 (45, 15, 0, 10, '#00ff00'),
 (45, 30, 0, 10, '#0000ff'),
 (45, 45, 0, 10, '#ffff00'),
 (45, 60, 0, 10, '#ff00ff'),
 (60, 0, 0, 10, '#ff0000'),
 (60, 15, 0, 10, '#00ff00'),
 (60, 30, 0, 10, '#0000ff'),
 (60, 45, 0, 10, '#ffff00'),
 (60, 60, 0, 10, '#ff00ff'),
 (0, 0, 15, 10, '#ff0000'),
 (0, 15, 15, 10, '#00ff00'),
 (0, 30, 15, 10, '#0000ff'),
 (0, 45, 15, 10, '#ffff00'),
 (0, 60, 15, 10, '#ff00ff'),
 (15, 0, 15, 10, '#ff0000'),
 (15, 15, 15, 10, '#00ff00'),
 (15, 30, 15, 10, '#0000ff'),
 (15, 45, 15, 10, '#ffff00'),
 (15, 60, 15, 10, '#ff00ff'),
 (30, 0, 15, 10, '#ff0000'),
 (30, 15, 15, 10, '#00ff00'),
 (30, 30, 15, 10, '#0000ff'),
 (30, 45, 15, 10, '#ffff00'),
 (30, 60, 15, 10, '#ff00ff'),
 (45, 0, 15, 10, '#ff0000'),
 (45, 15, 15, 10, '#00ff00'),
 (45, 30, 15, 10, '#0000ff'),
 (45, 45, 15, 10, '#ffff00'),
 (45, 60, 15, 10, '#ff00ff'),
 (60, 0, 15, 10, '#ff0000'),
 (60, 15, 15, 10, '#00ff00'),
 (60, 30, 15, 10, '#0000ff'),
 (60, 45, 15, 10, '#ffff00'),
 (60, 60, 15, 10, '#ff00ff'),
        ]:
            xm = x + s; ym = y + s; zm = z + s
            B.append((x, y, z, xm, ym, zm,
                      x + s*0.5, y + s*0.5, z + s*0.5,
                      c))
        self.blocks = B

    def motion_update(self):
        t = _time()
        dt = t - self._t
        self._t = t
        roll, pitch, yaw = _get_att()
        self.roll = -pitch
        self.pitch = roll + self.Plus
        self.yaw = -yaw

    def render(self):
        posx, posy, posz = self.pos
        w, h, res = self.size.w, self.size.h, self.res
        fov = self.fov_tan
        aspect = w / h
        full_sx = _int(w // res)
        sy = _int(h // res)
        if VR:
            sx = full_sx // 2
            x_offsets = (0, w // 2)
        else:
            sx = full_sx
            x_offsets = (0,)
        inv_sx2 = 2.0 * res / w
        inv_sy2 = 2.0 * res / h
        cr, sr = _cos(self.roll), _sin(self.roll)
        cp, sp = _cos(self.pitch), _sin(self.pitch)
        cy, sy_ = _cos(self.yaw),   _sin(self.yaw)
        r00 =  cy*cr + sy_*sp*sr;  r01 = -cy*sr + sy_*sp*cr;  r02 =  sy_*cp
        r10 =  cp*sr;              r11 =  cp*cr;              r12 = -sp
        r20 = -sy_*cr + cy*sp*sr;  r21 =  sy_*sr + cy*sp*cr;  r22 =  cy*cp
        blocks = self.blocks
        px_ndc = [(px + 0.5) * inv_sx2 - 1.0 for px in _range(_int(sx))]
        px_mul = [val * fov * aspect for val in px_ndc]
        py_ndc = [(py + 0.5) * inv_sy2 - 1.0 for py in _range(sy)]
        py_mul = [val * fov for val in py_ndc]
        #sun
        light = (1.0, 1.2, -0.8)
        MIN = 0.3
        inv_len_l = 1.0/_sqrt(light[0]**2 + light[1]**2 + light[2]**2)
        light = (light[0]*inv_len_l, light[1]*inv_len_l, light[2]*inv_len_l)
        for x_off in x_offsets:
            hit_buf   = [[False]*sx for _ in _range(sy)]
            color_buf = [[None ]*sx for _ in _range(sy)]
            for py, ny in _enumerate(py_mul):
                for px, nx in _enumerate(px_mul):
                    dx = r00*nx + r01*ny + r02
                    dy = r10*nx + r11*ny + r12
                    dz = r20*nx + r21*ny + r22
                    inv_len = 1.0 / _sqrt(dx*dx + dy*dy + dz*dz)
                    rx, ry, rz = dx*inv_len, dy*inv_len, dz*inv_len
                    inv_rx = 1.0/rx if _fabs(rx) > eps else _float('inf')
                    inv_ry = 1.0/ry if _fabs(ry) > eps else _float('inf')
                    inv_rz = 1.0/rz if _fabs(rz) > eps else _float('inf')
                    x0 = posx + rx;  y0 = posy + ry;  z0 = posz + rz
                    best_d = _float('inf')
                    best_col = None
                    best_t = 0.0
                    best_bd = None
                    for xmin, ymin, zmin, xmax, ymax, zmax, cx, cy2, cz, col in blocks:
                        t1 = (xmin - x0) * inv_rx;  t2 = (xmax - x0) * inv_rx
                        if t1 > t2: t1, t2 = t2, t1
                        u1 = (ymin - y0) * inv_ry;  u2 = (ymax - y0) * inv_ry
                        if u1 > u2: u1, u2 = u2, u1
                        v1 = (zmin - z0) * inv_rz;  v2 = (zmax - z0) * inv_rz
                        if v1 > v2: v1, v2 = v2, v1
                        t_near = t1 if t1>u1 else u1
                        if v1 > t_near: t_near = v1
                        t_far  = t2 if t2<u2 else u2
                        if v2 < t_far:  t_far  = v2
                        if t_near <= t_far and t_far > 0.0:
                            dx2 = cx - x0;  dy2 = cy2 - y0;  dz2 = cz - z0
                            d2 = dx2*dx2 + dy2*dy2 + dz2*dz2
                            if d2 < best_d:
                                best_d     = d2
                                best_col   = col
                                best_t     = t_near
                                best_bd    = (xmin, ymin, zmin, xmax, ymax, zmax)
                    if not best_col:
                        continue
                    c = best_col
                    if shading:
                        hx = x0 + rx * best_t
                        hy = y0 + ry * best_t
                        hz = z0 + rz * best_t
                        xmin, ymin, zmin, xmax, ymax, zmax = best_bd
                        bs  = _max(xmax-xmin, ymax-ymin, zmax-zmin)
                        tol = bs * 1e-3
                        if   _abs(hx - xmin) < tol:  normal = (-1,  0,  0)
                        elif _abs(hx - xmax) < tol:  normal = ( 1,  0,  0)
                        elif _abs(hy - ymin) < tol:  normal = ( 0, -1,  0)
                        elif _abs(hy - ymax) < tol:  normal = ( 0,  1,  0)
                        elif _abs(hz - zmin) < tol:  normal = ( 0,  0, -1)
                        elif _abs(hz - zmax) < tol:  normal = ( 0,  0,  1)
                        else:                       normal = ( 0,  0,  0)  
                        intensity = normal[0]*light[0] + normal[1]*light[1] + normal[2]*light[2]
                        if intensity < 0: intensity = 0.0
                        r0 = _int(best_col[1:3], 16)
                        g0 = _int(best_col[3:5], 16)
                        b0 = _int(best_col[5:7], 16)
                        intensity = MIN + (1.0 - MIN)*intensity
                        r = _int(r0 * intensity)
                        g = _int(g0 * intensity)
                        b = _int(b0 * intensity)
                        c = f"#{r:02x}{g:02x}{b:02x}"
                    hit_buf  [py][px] = True
                    color_buf[py][px] = c
                    if edges:
                        xmin, ymin, zmin, xmax, ymax, zmax = best_bd
                        hx = x0 + rx * best_t
                        hy = y0 + ry * best_t
                        hz = z0 + rz * best_t
                        bs  = _max(xmax-xmin, ymax-ymin, zmax-zmin)
                        tol = bs * 0.01
                        cnt = (
                            (_fabs(hx-xmin)<tol or _fabs(hx-xmax)<tol) +
                            (_fabs(hy-ymin)<tol or _fabs(hy-ymax)<tol) +
                            (_fabs(hz-zmin)<tol or _fabs(hz-zmax)<tol)
                        )
                        if cnt >= 2:
                            color_buf[py][px] = edges_color
            for py in _range(sy):
                row_hit = hit_buf[py]
                row_col = color_buf[py]
                y0 = py * res
                y1 = y0 + res
                px = 0
                coords = [None, None, None, None]
                while px < sx:
                    if not row_hit[px]:
                        px += 1
                        continue
                    run_color = row_col[px]
                    start = px
                    px += 1
                    while px < sx and row_hit[px] and row_col[px] == run_color:
                        px += 1
                    x0 = x_off + start * res
                    x1 = x_off + px    * res
                    coords[0] = (x0, y0)
                    coords[1] = (x1, y0)
                    coords[2] = (x0, y1)
                    coords[3] = (x1, y1)
                    _fill(run_color)
                    _triangle_strip(coords)
        
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
