#Args:
shading = True
VR = False
debug = True
edges = True
edges_color = "#000000"

from scene import *
from math import sin, cos, tan, radians, sqrt
import motion, time
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
_enumerate = enumerate
_triangle_strip = triangle_strip

CUBE_EDGES = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
]

# Cube face triangles for rasterization.
# Vertex layout for cube (xmin,ymin,zmin)→(xmax,ymax,zmax):
#   0:(xmin,ymin,zmin)  1:(xmax,ymin,zmin)  2:(xmax,ymax,zmin)  3:(xmin,ymax,zmin)
#   4:(xmin,ymin,zmax)  5:(xmax,ymin,zmax)  6:(xmax,ymax,zmax)  7:(xmin,ymax,zmax)
# Each entry: (v0_idx, v1_idx, v2_idx, normal_x, normal_y, normal_z)
CUBE_FACE_TRIS = (
    (4, 5, 7,  0,  0,  1),  # Front  z+
    (5, 6, 7,  0,  0,  1),
    (1, 0, 2,  0,  0, -1),  # Back   z-
    (0, 3, 2,  0,  0, -1),
    (5, 1, 6,  1,  0,  0),  # Right  x+
    (1, 2, 6,  1,  0,  0),
    (0, 4, 3, -1,  0,  0),  # Left   x-
    (4, 7, 3, -1,  0,  0),
    (7, 6, 3,  0,  1,  0),  # Top    y+
    (6, 2, 3,  0,  1,  0),
    (0, 1, 4,  0, -1,  0),  # Bottom y-
    (1, 5, 4,  0, -1,  0),
)

# Normalized sun direction: (1, 1.2, -0.8) / ||(1, 1.2, -0.8)||
_SUN_X = 0.5773502691896258
_SUN_Y = 0.6928203230275509
_SUN_Z = -0.4618802153517006
_MIN_INTENSITY = 0.3
_NEAR_CLIP = 0.1


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

        # Pre-compute flat-shaded triangle list (12 triangles per block).
        all_triangles = []
        for x, y, z, s, c in raw_blocks:
            xmin = _float(x); ymin = _float(y); zmin = _float(z)
            xmax = xmin + s;  ymax = ymin + s;  zmax = zmin + s
            r0, g0, b0 = _parse_hex(c)
            verts = (
                (xmin, ymin, zmin), (xmax, ymin, zmin),
                (xmax, ymax, zmin), (xmin, ymax, zmin),
                (xmin, ymin, zmax), (xmax, ymin, zmax),
                (xmax, ymax, zmax), (xmin, ymax, zmax),
            )
            for i0, i1, i2, nx, ny, nz in CUBE_FACE_TRIS:
                if shading:
                    dot = nx * _SUN_X + ny * _SUN_Y + nz * _SUN_Z
                    if dot < 0.0:
                        dot = 0.0
                    intensity = _MIN_INTENSITY + (1.0 - _MIN_INTENSITY) * dot
                    tri_color = '#{:02x}{:02x}{:02x}'.format(
                        _int(r0 * intensity),
                        _int(g0 * intensity),
                        _int(b0 * intensity),
                    )
                else:
                    tri_color = c
                all_triangles.append((verts[i0], verts[i1], verts[i2], tri_color))
        self.all_triangles = all_triangles

    def motion_update(self):
        t = _time()
        self._t = t
        roll, pitch, yaw = _get_att()
        self.roll = -pitch
        self.pitch = roll + self.Plus
        self.yaw = -yaw

    def render(self):
        posx, posy, posz = self.pos
        w, h = self.size.w, self.size.h
        half_h = h * 0.5
        fov_tan = self.fov_tan
        aspect = w / h

        if VR:
            # Each eye occupies half the screen width; project into w/4 half-viewport.
            half_w = w * 0.25
            x_offsets = (0.0, w * 0.5)
        else:
            half_w = w * 0.5
            x_offsets = (0.0,)

        fov_aspect = fov_tan * aspect

        # --- Rotation matrix (camera-to-world, same as original) ---
        roll_ = self.roll; pitch_ = self.pitch; yaw_ = self.yaw
        cr, sr = _cos(roll_), _sin(roll_)
        cp, sp = _cos(pitch_), _sin(pitch_)
        cy, sy_ = _cos(yaw_), _sin(yaw_)
        r00 = cy*cr + sy_*sp*sr;  r01 = -cy*sr + sy_*sp*cr;  r02 = sy_*cp
        r10 = cp*sr;              r11 = cp*cr;                r12 = -sp
        r20 = -sy_*cr + cy*sp*sr; r21 = sy_*sr + cy*sp*cr;   r22 = cy*cp

        # World-to-camera transform = transpose of camera-to-world (orthonormal).
        t00 = r00; t01 = r10; t02 = r20
        t10 = r01; t11 = r11; t12 = r21
        t20 = r02; t21 = r12; t22 = r22

        all_triangles = self.all_triangles
        do_edges = edges

        # --- Project triangle vertices into camera space, then screen space ---
        projected = []
        for v0, v1, v2, tri_color in all_triangles:
            dx = v0[0] - posx; dy = v0[1] - posy; dz = v0[2] - posz
            c0x = t00*dx + t01*dy + t02*dz
            c0y = t10*dx + t11*dy + t12*dz
            c0z = t20*dx + t21*dy + t22*dz
            if c0z < _NEAR_CLIP:
                continue

            dx = v1[0] - posx; dy = v1[1] - posy; dz = v1[2] - posz
            c1x = t00*dx + t01*dy + t02*dz
            c1y = t10*dx + t11*dy + t12*dz
            c1z = t20*dx + t21*dy + t22*dz
            if c1z < _NEAR_CLIP:
                continue

            dx = v2[0] - posx; dy = v2[1] - posy; dz = v2[2] - posz
            c2x = t00*dx + t01*dy + t02*dz
            c2y = t10*dx + t11*dy + t12*dz
            c2z = t20*dx + t21*dy + t22*dz
            if c2z < _NEAR_CLIP:
                continue

            inv_z0 = 1.0 / c0z
            s0x = half_w + (c0x * inv_z0 / fov_aspect) * half_w
            s0y = half_h + (c0y * inv_z0 / fov_tan) * half_h

            inv_z1 = 1.0 / c1z
            s1x = half_w + (c1x * inv_z1 / fov_aspect) * half_w
            s1y = half_h + (c1y * inv_z1 / fov_tan) * half_h

            inv_z2 = 1.0 / c2z
            s2x = half_w + (c2x * inv_z2 / fov_aspect) * half_w
            s2y = half_h + (c2y * inv_z2 / fov_tan) * half_h

            avg_depth = (c0z + c1z + c2z) / 3.0
            projected.append((avg_depth, s0x, s0y, s1x, s1y, s2x, s2y, tri_color))

        # --- Painter's algorithm: sort back-to-front (descending depth) ---
        projected.sort(reverse=True)

        # --- Draw triangles for each viewport (mono or VR) ---
        for x_off in x_offsets:
            for item in projected:
                _, s0x, s0y, s1x, s1y, s2x, s2y, tri_color = item
                _fill(tri_color)
                _triangle_strip([
                    (x_off + s0x, s0y),
                    (x_off + s1x, s1y),
                    (x_off + s2x, s2y),
                ])
                if do_edges:
                    stroke(edges_color)
                    stroke_weight(1)
                    line(x_off + s0x, s0y, x_off + s1x, s1y)
                    line(x_off + s1x, s1y, x_off + s2x, s2y)
                    line(x_off + s2x, s2y, x_off + s0x, s0y)

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
