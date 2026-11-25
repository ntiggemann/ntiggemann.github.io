###############################
# This script plots the tropical surface defined by a trivariate tropical polynomial
# The tropical numbers are defined with max as addition
# v1.4
###############################

''' Preamble '''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import re
import itertools

''' INPUT '''
# This is the only stuff you should change
the_polynomial = "x+y+z+-1z2+0"

''' Plotting range '''
# Only the part of the surface inside this box is computed and drawn
# If fit_range is True, the box is chosen automatically: a cube containing all vertices, hence all bounded edges
# and polygons, with some space around them for the unbounded parts. The boundaries below are then ignored.
# (Needs scipy.)

fit_range = False

a_x = -6  # Boundary x-axis left
b_x = 6   # Boundary x-axis right
a_y = -6  # Boundary y-axis left
b_y = 6   # Boundary y-axis right
a_z = -6
b_z = 6

'''END OF INPUT'''

''' If something weird happens, read this:'''
## Non-symbolic computing problems: One does not simply check floats for equality
# Two values count as equal if they differ by less than eps * (size of the numbers involved).
# The default should be fine. If a plot looks wrong for a polynomial of very high degree or with
# huge coefficients, this is the first thing to change and check again.

eps = 1e-9  # Relative tolerance

use_min = False # This is for later possible implementation of min-support. So far, it only affects the parser

''' Code - functions '''
# Don't change stuff after here - unless you really know what you are doing
def parse_polynomial_trop(poly_str, use_min):
    """
    Parses a polynomial string in the variables x, y, z and creates a dictionary.

    Args:
        poly_str (str): The polynomial string, e.g., "3.2x2y + yx + -1z3 + 9".
            Terms are separated by '+'. A term is an optional coefficient followed by variables with optional exponents.
            A missing coefficient means 0 (the tropical 1). A '-' in front of a coefficient makes it negative,
            so "x+-1" and "x-1" mean the same.
        use_min (bool): If a monomial occurs twice, keep the min of the coefficients (else the max)

    Returns:
        dict: (i,j,k) -> coefficient of x^i y^j z^k
    """
    s = poly_str
    for char in ' *^()':
        s = s.replace(char, '')
    s = s.replace('+-', '-')
    if not s:
        raise ValueError("The polynomial is empty.")

    # sign, coefficient, variables with exponents
    term_pattern = re.compile(r'([+-]?)(\d+\.?\d*|\.\d+)?((?:[xyz]\d*)*)')

    coefficients = {}
    pos = 0
    while pos < len(s):
        match = term_pattern.match(s, pos)
        sign, number, variables = match.groups()
        if (pos > 0 and not sign) or (not number and not variables):
            raise ValueError(f"Cannot read the polynomial at position {pos}: '{s[pos:]}'")
        if sign == '-' and not number:
            raise ValueError(f"A '-' has to be followed by a number (the coefficient), at position {pos}: '{s[pos:]}'")

        coefficient = float(sign + number) if number else 0.0

        exponents = {'x': 0, 'y': 0, 'z': 0}
        for var, exp in re.findall(r'([xyz])(\d*)', variables):
            exponents[var] += int(exp) if exp else 1
        key = (exponents['x'], exponents['y'], exponents['z'])

        # Store coefficient (tropical sum if the same monomial occurs more than once)
        if key in coefficients:
            if use_min:
                coefficients[key] = min(coefficient, coefficients[key])
            else:
                coefficients[key] = max(coefficient, coefficients[key])
        else:
            coefficients[key] = coefficient
        pos = match.end()

    return dict(sorted(coefficients.items(), key=lambda item: item[0], reverse=True))


def maxMonomsAt(point, exps, coeffs, tol):
    '''
    Args:
        point (np.array of length 3)
        exps (n x 3 np.array): exponents of the monomials
        coeffs (np.array of length n): coefficients of the monomials
        tol (float): two values count as equal if they differ by less than tol
    Returns:
        np.array of int: indices of the monomials attaining the maximum at point
    '''
    vals = coeffs + exps @ point
    return np.flatnonzero(vals >= vals.max() - tol)


def affineDim(points):
    '''
    Returns the dimension of the affine span of the rows of points
    '''
    if len(points) < 2:
        return 0
    return np.linalg.matrix_rank(points[1:] - points[0])


def planeBasis(normal):
    '''
    Returns an orthonormal basis u, v of the plane orthogonal to normal
    '''
    n = normal / np.linalg.norm(normal)
    e = np.zeros(3)
    e[np.argmin(np.abs(n))] = 1.0
    u = e - np.dot(e, n) * n
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v


def clipPolygon(poly, a, b, tol):
    '''
    One step of Sutherland-Hodgman: intersects a convex polygon with the half plane {y : a.y <= b}
    Also works for degenerate polygons (a segment given by two points, a single point)

    Args:
        poly (list of 2-dim np.arrays): vertices of a convex polygon, in order
        a (2-dim np.array of norm 1), b (float): the half plane
        tol (float): points with a.y <= b + tol count as inside
    Returns:
        list of 2-dim np.arrays: vertices of the intersection, in order
    '''
    result = []
    d = [np.dot(a, p) - b for p in poly]
    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        inside_i = d[i] <= tol
        inside_j = d[j] <= tol
        if inside_i:
            result.append(poly[i])
        if inside_i != inside_j:  # The side crosses the line a.y = b
            t = min(max(d[i] / (d[i] - d[j]), 0.0), 1.0)
            result.append(poly[i] + t * (poly[j] - poly[i]))
    return result


def cleanPolygon(poly, tol):
    '''
    Removes vertices that occur twice and vertices that lie on a side of the polygon (up to tol)
    '''
    pts = []
    for p in poly:
        if not pts or np.max(np.abs(p - pts[-1])) > tol:
            pts.append(p)
    while len(pts) > 1 and np.max(np.abs(pts[0] - pts[-1])) <= tol:
        pts.pop()
    changed = True
    while changed and len(pts) > 2:
        changed = False
        for i in range(len(pts)):
            p, q, r = pts[i - 1], pts[i], pts[(i + 1) % len(pts)]
            pr = r - p
            length = np.linalg.norm(pr)
            # distance of q to the line through p and r
            if length <= tol or abs(pr[0] * (q - p)[1] - pr[1] * (q - p)[0]) / length <= tol:
                pts.pop(i)
                changed = True
                break
    return pts


def polygonArea(poly):
    '''
    Shoelace formula for a polygon in the plane
    '''
    area = 0.0
    for i in range(len(poly)):
        p, q = poly[i], poly[(i + 1) % len(poly)]
        area += p[0] * q[1] - p[1] * q[0]
    return abs(area) / 2


def cellPiece(l1, l2, exps, coeffs, lower, upper, tol_val, tol_len):
    '''
    Computes the set {p in the box [lower, upper] : m1(p) = m2(p) >= m(p) for all monomials m}, where m1, m2 are
    the monomials with indices l1, l2. This set is convex: a polygon, a segment, a point or empty.

    Args:
        l1, l2 (int): indices of two monomials
        exps, coeffs: the polynomial, see maxMonomsAt
        lower, upper (np.arrays of length 3): corners of the box
        tol_val, tol_len (float): tolerances for values and lengths
    Returns:
        (int, list of 3-dim np.arrays): (2, vertices of the polygon in order) or (1, the two end points of the segment)
        or (0, [the point]) or (-1, []) if the set is empty
    '''
    normal = exps[l1] - exps[l2]
    offset = coeffs[l2] - coeffs[l1]  # The plane m1 = m2 is {p : normal.p = offset}
    normal_sq = np.dot(normal, normal)
    center = (lower + upper) / 2
    radius = np.linalg.norm(upper - lower) / 2  # Every point of the box has at most this distance to center

    # The plane does not even meet the ball around the box
    if abs(offset - np.dot(normal, center)) / np.sqrt(normal_sq) > radius + tol_len:
        return -1, []

    # Coordinates in the plane: p = p0 + s*u + t*v
    p0 = center + (offset - np.dot(normal, center)) / normal_sq * normal
    u, v = planeBasis(normal)
    uv = np.column_stack([u, v])

    # Half spaces g.p <= h: the box and m1 >= m for all other monomials m
    others = [l for l in range(len(coeffs)) if l != l1 and l != l2]
    box_g = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
    box_h = np.array([upper[0], -lower[0], upper[1], -lower[1], upper[2], -lower[2]])
    g = np.vstack([box_g, exps[others] - exps[l1]])
    h = np.concatenate([box_h, coeffs[l1] - coeffs[others]])

    # In plane coordinates
    g_2d = g @ uv
    h_2d = h - g @ p0

    # Start with a square containing the intersection of the plane with the box
    R = radius + 1
    poly = [np.array([-R, -R]), np.array([R, -R]), np.array([R, R]), np.array([-R, R])]

    for a_2d, b_2d in zip(g_2d, h_2d):
        norm_a = np.linalg.norm(a_2d)
        if norm_a < 1e-9:
            # The half space contains the whole plane or nothing of it
            if b_2d < -tol_val:
                return -1, []
            continue
        poly = clipPolygon(poly, a_2d / norm_a, b_2d / norm_a, tol_len)
        if not poly:
            return -1, []

    poly = cleanPolygon(poly, tol_len)
    if len(poly) >= 3 and polygonArea(poly) > 100 * tol_len * R:
        return 2, [p0 + uv @ q for q in poly]

    # Not 2-dimensional: a segment or a point. Take the two points with the largest distance.
    best = (0.0, poly[0], poly[0])
    for q1, q2 in itertools.combinations(poly, 2):
        dist = np.linalg.norm(q1 - q2)
        if dist > best[0]:
            best = (dist, q1, q2)
    if best[0] > 10 * tol_len:
        return 1, [p0 + uv @ best[1], p0 + uv @ best[2]]
    return 0, [p0 + uv @ ((best[1] + best[2]) / 2)]


def tropicalVertices(exps, coeffs, eps):
    '''
    Computes the vertices of the tropical surface. They correspond to the upper facets of the lifted Newton polytope
    conv{(a, c_a)}: an upper facet with outer normal (p, 1) gives the vertex p, as there the monomials of the facet are maximal.
    If the Newton polytope is not 3-dimensional, the surface has no vertices: the cells of minimal dimension are lines or
    planes. Then for each of them the point closest to the origin is returned.

    Args:
        exps, coeffs: the polynomial, see maxMonomsAt
        eps (float): relative tolerance
    Returns:
        (m x 3 np.array, int): the points and the dimension of the Newton polytope
    '''
    from scipy.spatial import ConvexHull, QhullError

    # Coordinates of the exponents in an orthonormal basis of the directions of the Newton polytope
    _, sing, vt = np.linalg.svd(exps[1:] - exps[0])
    dim = int(np.sum(sing > 1e-9))
    basis = vt[:dim].T
    y = (exps - exps[0]) @ basis
    lifted = np.column_stack([y, coeffs])
    scale = 1 + np.max(np.abs(lifted))

    candidates = []
    if np.linalg.matrix_rank(lifted[1:] - lifted[0], tol=eps * scale) == dim + 1:
        try:
            hull = ConvexHull(lifted)
        except QhullError:
            hull = ConvexHull(lifted, qhull_options='QJ')  # Nearly degenerate: joggle the input
        for equation in hull.equations:  # equation = (outer normal, offset)
            if equation[dim] > 1e-9:  # upper facet
                candidates.append(equation[:dim] / equation[dim])
    if not candidates:
        # The lifted points lie in a hyperplane, i.e. there is only one cell of minimal dimension,
        # the set where all monomials are maximal: Solve c_a + <y_a, q> = t for all a
        M = np.column_stack([y, -np.ones(len(y))])
        candidates = [np.linalg.lstsq(M, -coeffs, rcond=None)[0][:dim]]

    # Back to R^3, check that the cell really has minimal dimension, delete points that occur twice
    points = []
    for q in candidates:
        p = basis @ q
        vals = coeffs + exps @ p
        monoms = np.flatnonzero(vals >= vals.max() - 1e-6 * (1 + np.max(np.abs(vals))))
        if affineDim(exps[monoms]) < dim:
            continue
        if not any(np.max(np.abs(p - p2)) <= 1e-6 * (1 + np.max(np.abs(p))) for p2 in points):
            points.append(p)
    return np.array(points), dim


def onSegment(point, p, q, tol):
    '''
    True if point has distance at most tol to the segment [p, q]
    '''
    d = q - p
    length_sq = np.dot(d, d)
    t = 0.0 if length_sq == 0 else min(max(np.dot(point - p, d) / length_sq, 0.0), 1.0)
    return np.linalg.norm(p + t * d - point) <= tol


''' Code - do stuff '''

''' Define the polynomial '''
A = parse_polynomial_trop(the_polynomial, use_min)

if (len(A) < 2):
    raise ValueError("Please enter a polynomial with at least two monomials.")

exps = np.array(list(A.keys()), dtype=float)      # Row l: exponent of the l-th monomial
coeffs = np.array(list(A.values()), dtype=float)  # Entry l: coefficient of the l-th monomial
num_monoms = len(A)

''' Choose the plotting range automatically '''
if fit_range:
    vertices, newton_dim = tropicalVertices(exps, coeffs, eps)
    low = vertices.min(axis=0)
    high = vertices.max(axis=0)
    extent = np.max(high - low)
    half = extent / 2 + max(1.0, 0.25 * extent)  # A cube, with some space around the vertices
    mid = (low + high) / 2
    step = 10.0 ** (np.floor(np.log10(2 * half)) - 1)  # Round the boundaries outwards to two significant digits
    a_x, a_y, a_z = (round(float(np.floor(v / step + 1e-9) * step), 10) + 0.0 for v in mid - half)
    b_x, b_y, b_z = (round(float(np.ceil(v / step - 1e-9) * step), 10) + 0.0 for v in mid + half)
    if newton_dim == 3:
        print(f"{len(vertices)} vertices.", end=" ")
    else:
        print(f"The Newton polytope is {newton_dim}-dimensional, so there are no vertices; {len(vertices)} cells of minimal dimension.", end=" ")
    print(f"Plotting range: x in [{a_x:.3g}, {b_x:.3g}], y in [{a_y:.3g}, {b_y:.3g}], z in [{a_z:.3g}, {b_z:.3g}]")

if (a_x >= b_x) or (a_y >= b_y) or (a_z >= b_z):
    raise ValueError("Your plotting range has too few dimensions or the boundarys dont make sense (e.g. intervall [2,1]). If you want to plot something 2-dim there is an extra script for that.")

lower = np.array([a_x, a_y, a_z], dtype=float)
upper = np.array([b_x, b_y, b_z], dtype=float)

# Tolerances: for values of monomials and for lengths
degF = int(np.max(exps.sum(axis=1)))
tol_val = eps * (1 + np.max(np.abs(coeffs)) + degF * np.max(np.abs(np.concatenate([lower, upper]))))
tol_len = eps * (1 + np.linalg.norm(upper - lower) / 2 + np.linalg.norm((lower + upper) / 2))

''' Compute the 2-cells '''
# For each pair of monomials (m1, m2): The set {p in plotting range : m1(p) = m2(p) >= m(p) for all m} is convex.
# If it is 2-dimensional, it is (the part in the plotting range of) the 2-cell dual to the edge [m1,m2]
# of the subdivision of the Newton polytope.
# If it is a segment or a point, it is either part of a lower dimensional cell (if [m1,m2] is not an edge of the
# subdivision) or a 2-cell that only touches the boundary of the plotting range.

cells = []       # In here go dicts with the vertices (in order) of a polygon, the maximal monomials on it and the normal vector
seen_cells = set()
small_pieces = [] # In here go the segments and points, together with the 2-cell they belong to (None if [m1,m2] is no edge of the subdivision)

for l1, l2 in itertools.combinations(range(num_monoms), 2):
    dim, points = cellPiece(l1, l2, exps, coeffs, lower, upper, tol_val, tol_len)
    if dim == -1:
        continue
    normal = exps[l1] - exps[l2]

    if dim == 2:
        # All monomials that are maximal on the polygon. If m1, m2, m3 lie on a line and are maximal on the same polygon,
        # the pairs (m1,m2), (m1,m3), (m2,m3) give the same polygon. Draw it only once.
        monoms = maxMonomsAt(np.mean(points, axis=0), exps, coeffs, tol_val)
        key = tuple(monoms)
        if key in seen_cells:
            continue
        seen_cells.add(key)
        cells.append({'vertices': points, 'monoms': monoms, 'normal': normal})
    else:
        # Is [m1,m2] an edge of the subdivision, i.e. is the 2-cell 2-dimensional in a slightly bigger box?
        touching_cell = None
        big_dim, big_points = cellPiece(l1, l2, exps, coeffs, lower - 1, upper + 1, tol_val, tol_len)
        if big_dim == 2:
            monoms = maxMonomsAt(np.mean(big_points, axis=0), exps, coeffs, tol_val)
            touching_cell = {'monoms': monoms, 'normal': normal}
        small_pieces.append((dim, points, touching_cell))

''' Compute the edges '''
# A side of a polygon is an edge of the surface iff at least three affinely independent monomials are maximal on it.
# The other sides lie on the boundary of the plotting range.

edges = [] # In here go triples (point, point, colour)
seen_edges = set()
polygon_sides = [] # All sides of all polygons
for cell in cells:
    verts = cell['vertices']
    for i in range(len(verts)):
        p, q = verts[i], verts[(i + 1) % len(verts)]
        polygon_sides.append((p, q))
        monoms = maxMonomsAt((p + q) / 2, exps, coeffs, tol_val)
        if affineDim(exps[monoms]) >= 2:
            key = frozenset([tuple(np.round(p, 7)), tuple(np.round(q, 7))])
            if key not in seen_edges:
                seen_edges.add(key)
                edges.append((p, q, 'black'))

''' The parts that only touch the boundary of the plotting range '''
# A segment or point that is not contained in a side of a polygon is drawn extra:
# black if it lies on an edge or is a vertex of the surface, blue if it lies on a 2-cell.
lone_points = [] # In here go pairs (point, colour)
extra_segments = []
touching_cells = [] # 2-cells that only touch the boundary of the plotting range in such a segment or point. Only needed for the dual.
tol_cover = 10 * tol_len
for dim, points, touching_cell in sorted(small_pieces, key=lambda piece: -piece[0]):  # First the segments, then the points
    if dim == 1:
        p, q = points
        test_points = [p, (p + q) / 2, q]
    else:
        test_points = points
    # Already visible as part of a polygon
    if all(any(onSegment(x, s[0], s[1], tol_cover) for s in polygon_sides) for x in test_points):
        continue
    if touching_cell is not None:
        touching_cells.append(touching_cell)
    # Already drawn
    if all(any(onSegment(x, s[0], s[1], tol_cover) for s in extra_segments) for x in test_points):
        continue
    if dim == 0 and any(np.linalg.norm(points[0] - pt) <= tol_cover for pt, _ in lone_points):
        continue
    monoms = maxMonomsAt(np.mean(test_points, axis=0), exps, coeffs, tol_val)
    colour = 'black' if affineDim(exps[monoms]) >= 2 else 'b'
    if dim == 1:
        edges.append((p, q, colour))
        extra_segments.append((p, q))
    else:
        lone_points.append((points[0], colour))

''' The dual '''
# Each 2-cell is dual to an edge of the subdivision of the Newton polytope.
# Its maximal monomials lie on this edge; the end points are the first and the last of them in direction of the normal.
x_vals_edge_dual = []
y_vals_edge_dual = []
z_vals_edge_dual = []
seen_dual = set()
for cell in cells + touching_cells:
    key = tuple(cell['monoms'])
    if key in seen_dual:
        continue
    seen_dual.add(key)
    monom_exps = exps[cell['monoms']]
    proj = monom_exps @ cell['normal']
    first = monom_exps[np.argmin(proj)]
    last = monom_exps[np.argmax(proj)]
    x_vals_edge_dual.append([first[0], last[0]])
    y_vals_edge_dual.append([first[1], last[1]])
    z_vals_edge_dual.append([first[2], last[2]])

if not cells and not edges and not lone_points:
    print("The surface does not intersect your plotting range")
else:
    print(f"{len(cells)} 2-cells, {len(edges)} edges/segments and {len(lone_points)} isolated points in the plotting range")

''' Plotting'''
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12,6))

if cells:
    poly = Poly3DCollection([cell['vertices'] for cell in cells], facecolors='b', linewidths=0, alpha=0.5)
    ax1.add_collection3d(poly)

for p, q, colour in edges:
    ax1.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], color=colour)

for p, colour in lone_points:
    ax1.scatter([p[0]], [p[1]], [p[2]], color=colour)

ax1.set_xlim([a_x, b_x])
ax1.set_ylim([a_y, b_y])
ax1.set_zlim([a_z, b_z])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')

''' Draw the dual '''
# Create the coordinates for \Delta_deg
x_coords_dual = []
y_coords_dual = []
z_coords_dual = []

# Draw the points
for k in range(degF+1):
    for i in range(degF+1-k):
        for j in range(degF+1-k-i):
            x_coords_dual.append(i)
            y_coords_dual.append(j)
            z_coords_dual.append(k)

ax2.scatter(x_coords_dual, y_coords_dual, z_coords_dual, marker='o', color='b', alpha=0.5)

for i in range(len(x_vals_edge_dual)):
    ax2.plot(x_vals_edge_dual[i], y_vals_edge_dual[i], z_vals_edge_dual[i], marker='o', color='b')

plt.show()
