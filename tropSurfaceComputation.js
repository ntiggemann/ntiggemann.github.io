/*
Computation of tropical surfaces in R^3, given by a tropical polynomial in x, y, z.
This is a port of 3dTropPolyPlotLGSv1.4.py. It does not touch the website, so it can also be used with Node.js
(for tests). Internally everything is computed in the max-convention; the min-convention is handled by sign changes.

Idea of the algorithm:
For each pair of monomials m1, m2 the set
    {p in plotting range : m1(p) = m2(p) >= m(p) for all monomials m}
is convex. It is computed by cutting the plane m1 = m2 with half planes (Sutherland-Hodgman).
If it is 2-dimensional, it is the part in the plotting range of the 2-cell dual to the edge [m1,m2] of the
subdivision of the Newton polytope. If it is a segment or a point, it either lies on a lower dimensional cell or it is
a 2-cell that only touches the boundary of the plotting range.
*/

/*
Small vector helpers (vectors in R^3 are arrays of length 3)
*/

function dot3(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function add3(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

function sub3(a, b) {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function scale3(a, s) {
    return [s * a[0], s * a[1], s * a[2]];
}

function cross3(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

function norm3(a) {
    return Math.sqrt(dot3(a, a));
}

function maxAbs3(a) {
    return Math.max(Math.abs(a[0]), Math.abs(a[1]), Math.abs(a[2]));
}

function mean3(points) {
    const s = [0, 0, 0];
    for (const p of points) {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    return scale3(s, 1 / points.length);
}

function gcd(a, b) {
    a = Math.abs(a);
    b = Math.abs(b);

    while (b !== 0) {
        const r = a % b;
        a = b;
        b = r;
    }

    return a;
}

/*
Parsing
*/

function parsePolynomialTrop(polyStr, useMin) {
    /**
     * Parses a tropical polynomial in the variables x, y and z.
     *
     * Terms are separated by '+'. A term is an optional coefficient followed by variables with optional exponents.
     * A missing coefficient means 0 (the tropical 1). A '-' in front of a coefficient makes it negative,
     * so "x+-1" and "x-1" mean the same. Repeated variables are added up: xyx = x2y.
     *
     * Example input:
     *   "3.2x2y + yx + -1z3 + 9"
     *
     * @param {string} polyStr - The polynomial.
     * @param {boolean} useMin - If a monomial occurs twice, keep the min of the coefficients (else the max).
     * @returns {Array<[number, number, number, number]>} Array of [xExponent, yExponent, zExponent, coefficient].
     * @throws {Error} If the string cannot be read.
     */
    let s = polyStr.replace(/\s+/g, "");
    for (const char of ["*", "^", "(", ")"]) {
        s = s.replaceAll(char, "");
    }
    s = s.replaceAll("+-", "-");

    if (!s) {
        throw new Error("The polynomial is empty.");
    }

    // sign, coefficient, variables with exponents. Sticky: matches exactly at lastIndex.
    const termPattern = /([+-]?)(\d+\.?\d*|\.\d+)?((?:[xyz]\d*)*)/y;

    // "i,j,k" -> coefficient
    const coefficients = new Map();
    let pos = 0;

    while (pos < s.length) {
        termPattern.lastIndex = pos;
        const [, sign, number, variables] = termPattern.exec(s);

        if ((pos > 0 && !sign) || (!number && !variables)) {
            throw new Error(`Cannot read the polynomial at position ${pos}: '${s.slice(pos)}'`);
        }
        if (sign === "-" && !number) {
            throw new Error(`A '-' has to be followed by a number (the coefficient), at position ${pos}: '${s.slice(pos)}'`);
        }

        const coefficient = number ? parseFloat(sign + number) : 0;

        const exponents = { x: 0, y: 0, z: 0 };
        for (const [, variable, exponent] of variables.matchAll(/([xyz])(\d*)/g)) {
            exponents[variable] += exponent ? parseInt(exponent, 10) : 1;
        }
        const key = `${exponents.x},${exponents.y},${exponents.z}`;

        // Tropical sum if the same monomial occurs more than once
        if (coefficients.has(key)) {
            const old = coefficients.get(key);
            coefficients.set(key, useMin ? Math.min(coefficient, old) : Math.max(coefficient, old));
        } else {
            coefficients.set(key, coefficient);
        }

        pos = termPattern.lastIndex;
    }

    const monomials = [];
    for (const [key, coefficient] of coefficients) {
        const [i, j, k] = key.split(",").map(Number);
        monomials.push([i, j, k, coefficient]);
    }

    // Sort descending by x exponent, then y exponent, then z exponent
    monomials.sort((a, b) => (b[0] - a[0]) || (b[1] - a[1]) || (b[2] - a[2]));

    return monomials;
}

/*
Tropical helpers
*/

function monomialValues(p, exps, coeffs) {
    // The values of all monomials at the point p
    return coeffs.map((c, l) => c + dot3(exps[l], p));
}

function maxMonomsAt(p, exps, coeffs, tol) {
    /**
     * @returns {number[]} The indices of the monomials attaining the maximum at p (up to tol).
     */
    const values = monomialValues(p, exps, coeffs);
    const max = Math.max(...values);
    const result = [];
    for (let l = 0; l < values.length; l++) {
        if (values[l] >= max - tol) {
            result.push(l);
        }
    }
    return result;
}

function integerRank(vectors) {
    /**
     * The rank of a set of integer vectors in Z^3. Exact, since all computations are with integers.
     */
    const nonzero = vectors.filter(v => v[0] !== 0 || v[1] !== 0 || v[2] !== 0);
    if (nonzero.length === 0) {
        return 0;
    }
    const v1 = nonzero[0];
    let normal = null;
    for (const v of nonzero) {
        const c = cross3(v1, v);
        if (c[0] !== 0 || c[1] !== 0 || c[2] !== 0) {
            normal = c;
            break;
        }
    }
    if (normal === null) {
        return 1;
    }
    for (const v of nonzero) {
        if (dot3(normal, v) !== 0) {
            return 3;
        }
    }
    return 2;
}

function affineDim(points) {
    // The dimension of the affine span of integer points (e.g. exponents of monomials)
    if (points.length < 2) {
        return 0;
    }
    return integerRank(points.map(p => sub3(p, points[0])));
}

/*
Polygons in a plane
*/

function planeBasis(normal) {
    // An orthonormal basis u, v of the plane orthogonal to normal
    const n = scale3(normal, 1 / norm3(normal));
    const absN = n.map(Math.abs);
    const e = [0, 0, 0];
    e[absN.indexOf(Math.min(...absN))] = 1;
    let u = sub3(e, scale3(n, dot3(e, n)));
    u = scale3(u, 1 / norm3(u));
    const v = cross3(n, u);
    return [u, v];
}

function clipPolygon(poly, a, b, tol) {
    /**
     * One step of Sutherland-Hodgman: intersects a convex polygon with the half plane {y : a.y <= b}.
     * Also works for degenerate polygons (a segment given by two points, a single point).
     *
     * @param {number[][]} poly - Vertices of a convex polygon in the plane, in order.
     * @param {number[]} a - Normal vector of the half plane, of length 1.
     * @param {number} b - Offset of the half plane.
     * @param {number} tol - Points with a.y <= b + tol count as inside.
     * @returns {number[][]} Vertices of the intersection, in order.
     */
    const result = [];
    const d = poly.map(p => a[0] * p[0] + a[1] * p[1] - b);
    for (let i = 0; i < poly.length; i++) {
        const j = (i + 1) % poly.length;
        const insideI = d[i] <= tol;
        const insideJ = d[j] <= tol;
        if (insideI) {
            result.push(poly[i]);
        }
        if (insideI !== insideJ) { // The side crosses the line a.y = b
            const t = Math.min(Math.max(d[i] / (d[i] - d[j]), 0), 1);
            result.push([poly[i][0] + t * (poly[j][0] - poly[i][0]), poly[i][1] + t * (poly[j][1] - poly[i][1])]);
        }
    }
    return result;
}

function cleanPolygon(poly, tol) {
    // Removes vertices that occur twice and vertices that lie on a side of the polygon (up to tol)
    const close = (p, q) => Math.max(Math.abs(p[0] - q[0]), Math.abs(p[1] - q[1])) <= tol;
    const pts = [];
    for (const p of poly) {
        if (pts.length === 0 || !close(p, pts[pts.length - 1])) {
            pts.push(p);
        }
    }
    while (pts.length > 1 && close(pts[0], pts[pts.length - 1])) {
        pts.pop();
    }
    let changed = true;
    while (changed && pts.length > 2) {
        changed = false;
        for (let i = 0; i < pts.length; i++) {
            const p = pts[(i - 1 + pts.length) % pts.length];
            const q = pts[i];
            const r = pts[(i + 1) % pts.length];
            const pr = [r[0] - p[0], r[1] - p[1]];
            const length = Math.hypot(pr[0], pr[1]);
            // distance of q to the line through p and r
            if (length <= tol || Math.abs(pr[0] * (q[1] - p[1]) - pr[1] * (q[0] - p[0])) / length <= tol) {
                pts.splice(i, 1);
                changed = true;
                break;
            }
        }
    }
    return pts;
}

function polygonArea(poly) {
    // Shoelace formula
    let area = 0;
    for (let i = 0; i < poly.length; i++) {
        const p = poly[i];
        const q = poly[(i + 1) % poly.length];
        area += p[0] * q[1] - p[1] * q[0];
    }
    return Math.abs(area) / 2;
}

function polygonPerimeter(poly) {
    let perimeter = 0;
    for (let i = 0; i < poly.length; i++) {
        const p = poly[i];
        const q = poly[(i + 1) % poly.length];
        perimeter += Math.hypot(q[0] - p[0], q[1] - p[1]);
    }
    return perimeter;
}

function cellPiece(l1, l2, exps, coeffs, lower, upper, tolVal, tolLen) {
    /**
     * Computes the set {p in the box [lower, upper] : m1(p) = m2(p) >= m(p) for all monomials m}, where m1, m2 are
     * the monomials with indices l1, l2. This set is convex: a polygon, a segment, a point or empty.
     *
     * @returns {{dim: number, points: number[][]}} dim 2: vertices of the polygon in order; dim 1: the two end points
     *          of the segment; dim 0: the point; dim -1: empty.
     */
    const empty = { dim: -1, points: [] };
    const normal = sub3(exps[l1], exps[l2]);
    const offset = coeffs[l2] - coeffs[l1]; // The plane m1 = m2 is {p : normal.p = offset}
    const normalSq = dot3(normal, normal);
    const center = scale3(add3(lower, upper), 0.5);
    const radius = norm3(sub3(upper, lower)) / 2; // Every point of the box has at most this distance to center

    // The plane does not even meet the ball around the box
    if (Math.abs(offset - dot3(normal, center)) / Math.sqrt(normalSq) > radius + tolLen) {
        return empty;
    }

    // Coordinates in the plane: p = p0 + s*u + t*v
    const p0 = add3(center, scale3(normal, (offset - dot3(normal, center)) / normalSq));
    const [u, v] = planeBasis(normal);

    // Half spaces g.p <= h: the box and m1 >= m for all other monomials m
    const halfSpaces = [
        [[1, 0, 0], upper[0]], [[-1, 0, 0], -lower[0]],
        [[0, 1, 0], upper[1]], [[0, -1, 0], -lower[1]],
        [[0, 0, 1], upper[2]], [[0, 0, -1], -lower[2]]
    ];
    for (let l = 0; l < coeffs.length; l++) {
        if (l !== l1 && l !== l2) {
            halfSpaces.push([sub3(exps[l], exps[l1]), coeffs[l1] - coeffs[l]]);
        }
    }

    // Start with a square containing the intersection of the plane with the box
    const R = radius + 1;
    let poly = [[-R, -R], [R, -R], [R, R], [-R, R]];

    for (const [g, h] of halfSpaces) {
        const a = [dot3(g, u), dot3(g, v)];
        const b = h - dot3(g, p0);
        const normA = Math.hypot(a[0], a[1]);
        if (normA < 1e-9) {
            // The half space contains the whole plane or nothing of it
            if (b < -tolVal) {
                return empty;
            }
            continue;
        }
        poly = clipPolygon(poly, [a[0] / normA, a[1] / normA], b / normA, tolLen);
        if (poly.length === 0) {
            return empty;
        }
    }

    poly = cleanPolygon(poly, tolLen);
    const to3d = q => add3(p0, add3(scale3(u, q[0]), scale3(v, q[1])));

    // 2-dimensional, if it is not just a thin sliver (width about 2*area/perimeter)
    if (poly.length >= 3 && 2 * polygonArea(poly) / polygonPerimeter(poly) > 10 * tolLen) {
        return { dim: 2, points: poly.map(to3d) };
    }

    // A segment or a point: take the two points with the largest distance
    let best = [0, poly[0], poly[0]];
    for (let i = 0; i < poly.length; i++) {
        for (let j = i + 1; j < poly.length; j++) {
            const dist = Math.hypot(poly[i][0] - poly[j][0], poly[i][1] - poly[j][1]);
            if (dist > best[0]) {
                best = [dist, poly[i], poly[j]];
            }
        }
    }
    if (best[0] > 10 * tolLen) {
        return { dim: 1, points: [to3d(best[1]), to3d(best[2])] };
    }
    return { dim: 0, points: [to3d([(best[1][0] + best[2][0]) / 2, (best[1][1] + best[2][1]) / 2])] };
}

function onSegment(point, p, q, tol) {
    // True if point has distance at most tol to the segment [p, q]
    const d = sub3(q, p);
    const lengthSq = dot3(d, d);
    const t = lengthSq === 0 ? 0 : Math.min(Math.max(dot3(sub3(point, p), d) / lengthSq, 0), 1);
    return norm3(sub3(add3(p, scale3(d, t)), point)) <= tol;
}

function tolerances(exps, coeffs, lower, upper, eps) {
    // Tolerances for values of monomials and for lengths
    const degree = Math.max(...exps.map(e => e[0] + e[1] + e[2]));
    const maxCoeff = Math.max(...coeffs.map(Math.abs));
    const maxBound = Math.max(...lower.map(Math.abs), ...upper.map(Math.abs));
    return {
        tolVal: eps * (1 + maxCoeff + degree * maxBound),
        tolLen: eps * (1 + norm3(sub3(upper, lower)) / 2 + norm3(scale3(add3(lower, upper), 0.5)))
    };
}

/*
The surface in a box
*/

function computeSurfaceInBox(exps, coeffs, lower, upper, eps) {
    /**
     * Computes the part of the tropical surface inside the box [lower, upper] (max-convention).
     *
     * @returns {{cells: Object[], edges: Object[], lonePoints: Object[], dualEdges: Object[]}}
     *   cells: {vertices (in order), monoms (indices of the maximal monomials), normal, weight}
     *   edges: {p, q, onCell}, onCell = true if the segment lies on a 2-cell (and not on an edge of the surface)
     *   lonePoints: {p, onCell}
     *   dualEdges: {from, to, weight}, exponents of the end points of the edges of the subdivision dual to the 2-cells
     */
    const { tolVal, tolLen } = tolerances(exps, coeffs, lower, upper, eps);
    const n = coeffs.length;

    // The 2-cells
    const cells = [];
    const seenCells = new Set();
    const smallPieces = []; // segments and points, with the 2-cell they belong to (null if [m1,m2] is no edge of the subdivision)
    const lowerBig = lower.map(v => v - 1);
    const upperBig = upper.map(v => v + 1);

    for (let l1 = 0; l1 < n; l1++) {
        for (let l2 = l1 + 1; l2 < n; l2++) {
            const { dim, points } = cellPiece(l1, l2, exps, coeffs, lower, upper, tolVal, tolLen);
            if (dim === -1) {
                continue;
            }
            const normal = sub3(exps[l1], exps[l2]);

            if (dim === 2) {
                // If m1, m2, m3 lie on a line and are maximal on the same polygon, the pairs (m1,m2), (m1,m3), (m2,m3)
                // give the same polygon. Keep it only once.
                const monoms = maxMonomsAt(mean3(points), exps, coeffs, tolVal);
                const key = monoms.join(",");
                if (seenCells.has(key)) {
                    continue;
                }
                seenCells.add(key);
                cells.push({ vertices: points, monoms, normal });
            } else {
                // Is [m1,m2] an edge of the subdivision, i.e. is the 2-cell 2-dimensional in a slightly bigger box?
                let touchingCell = null;
                const big = cellPiece(l1, l2, exps, coeffs, lowerBig, upperBig, tolVal, tolLen);
                if (big.dim === 2) {
                    touchingCell = { monoms: maxMonomsAt(mean3(big.points), exps, coeffs, tolVal), normal };
                }
                smallPieces.push({ dim, points, touchingCell });
            }
        }
    }

    // The edges: sides of polygons on which at least three affinely independent monomials are maximal
    const edges = [];
    const seenEdges = new Set();
    const polygonSides = [];
    const roundedKey = p => p.map(t => (Math.round(t * 1e7) / 1e7 + 0).toString()).join(",");

    for (const cell of cells) {
        const verts = cell.vertices;
        for (let i = 0; i < verts.length; i++) {
            const p = verts[i];
            const q = verts[(i + 1) % verts.length];
            polygonSides.push([p, q]);
            const monoms = maxMonomsAt(scale3(add3(p, q), 0.5), exps, coeffs, tolVal);
            if (affineDim(monoms.map(l => exps[l])) >= 2) {
                const keys = [roundedKey(p), roundedKey(q)].sort();
                const key = keys.join("|");
                if (!seenEdges.has(key)) {
                    seenEdges.add(key);
                    edges.push({ p, q, onCell: false });
                }
            }
        }
    }

    // Segments and points that only touch the boundary of the box
    const lonePoints = [];
    const extraSegments = [];
    const touchingCells = [];
    const tolCover = 10 * tolLen;
    smallPieces.sort((a, b) => b.dim - a.dim); // First the segments, then the points

    for (const { dim, points, touchingCell } of smallPieces) {
        const testPoints = dim === 1 ? [points[0], scale3(add3(points[0], points[1]), 0.5), points[1]] : points;
        const coveredBy = segments => testPoints.every(x => segments.some(([p, q]) => onSegment(x, p, q, tolCover)));

        // Already visible as part of a polygon
        if (coveredBy(polygonSides)) {
            continue;
        }
        if (touchingCell !== null) {
            touchingCells.push(touchingCell);
        }
        // Already drawn
        if (coveredBy(extraSegments)) {
            continue;
        }
        if (dim === 0 && lonePoints.some(pt => norm3(sub3(pt.p, points[0])) <= tolCover)) {
            continue;
        }
        const monoms = maxMonomsAt(mean3(testPoints), exps, coeffs, tolVal);
        const onCell = affineDim(monoms.map(l => exps[l])) < 2;
        if (dim === 1) {
            edges.push({ p: points[0], q: points[1], onCell });
            extraSegments.push([points[0], points[1]]);
        } else {
            lonePoints.push({ p: points[0], onCell });
        }
    }

    // The dual: Each 2-cell is dual to an edge of the subdivision. Its maximal monomials lie on this edge;
    // the end points are the first and the last of them in direction of the normal.
    const dualEdges = [];
    const seenDual = new Set();
    for (const cell of cells.concat(touchingCells)) {
        const key = cell.monoms.join(",");
        const monomExps = cell.monoms.map(l => exps[l]);
        const proj = monomExps.map(e => dot3(e, cell.normal));
        const from = monomExps[proj.indexOf(Math.min(...proj))];
        const to = monomExps[proj.indexOf(Math.max(...proj))];
        const diff = sub3(to, from);
        const weight = gcd(gcd(diff[0], diff[1]), diff[2]);
        cell.weight = weight;
        if (!seenDual.has(key)) {
            seenDual.add(key);
            dualEdges.push({ from, to, weight });
        }
    }

    return { cells, edges, lonePoints, dualEdges };
}

/*
Vertices and the plotting range
*/

function minimalCellPoints(exps, coeffs, eps) {
    /**
     * Computes the vertices of the tropical surface. If the Newton polytope is not 3-dimensional, the surface has no
     * vertices: the cells of minimal dimension are lines or planes. Then for each of them the point closest to the
     * origin is returned.
     *
     * Method: The surface is computed in a box that surely contains all these points (bound by Cramer's rule and
     * Hadamard's inequality). The vertices are then the vertices of the polygons in the interior of the box.
     *
     * @returns {{points: number[][], newtonDim: number}}
     */
    const newtonDim = affineDim(exps);
    const deltaC = Math.max(...coeffs) - Math.min(...coeffs);
    let deg = 0;
    for (let m = 0; m < 3; m++) {
        const column = exps.map(e => e[m]);
        deg = Math.max(deg, Math.max(...column) - Math.min(...column));
    }
    const R = 3 * Math.sqrt(3) * deltaC * deg * deg + 1;
    const lower = [-R, -R, -R];
    const upper = [R, R, R];
    const { tolVal, tolLen } = tolerances(exps, coeffs, lower, upper, eps);

    const points = [];
    const addPoint = p => {
        const tol = Math.max(1e-6 * (1 + maxAbs3(p)), 100 * tolLen);
        if (!points.some(q => maxAbs3(sub3(p, q)) <= tol)) {
            points.push(p);
        }
    };

    const n = coeffs.length;
    for (let l1 = 0; l1 < n; l1++) {
        for (let l2 = l1 + 1; l2 < n; l2++) {
            const { dim, points: verts } = cellPiece(l1, l2, exps, coeffs, lower, upper, tolVal, tolLen);
            if (dim !== 2) {
                continue;
            }
            if (newtonDim === 3) {
                // Vertices of the polygon inside the box are vertices of the surface
                for (const p of verts) {
                    if (maxAbs3(p) < R - 0.5) {
                        addPoint(p);
                    }
                }
            } else if (newtonDim === 2) {
                // The cells of minimal dimension are lines: sides of polygons on which three monomials are maximal
                for (let i = 0; i < verts.length; i++) {
                    const p = verts[i];
                    const q = verts[(i + 1) % verts.length];
                    const monoms = maxMonomsAt(scale3(add3(p, q), 0.5), exps, coeffs, tolVal);
                    if (affineDim(monoms.map(l => exps[l])) >= 2) {
                        const d = sub3(q, p);
                        addPoint(sub3(p, scale3(d, dot3(p, d) / dot3(d, d))));
                    }
                }
            } else {
                // The cells of minimal dimension are the planes m1 = m2
                const normal = sub3(exps[l1], exps[l2]);
                addPoint(scale3(normal, (coeffs[l2] - coeffs[l1]) / dot3(normal, normal)));
            }
        }
    }
    return { points, newtonDim };
}

function chooseRange(points, dist, cube) {
    /**
     * A box containing all points with distance at least dist to its boundary, rounded outwards.
     *
     * @param {number[][]} points - The points (at least one).
     * @param {number} dist - Minimal distance of the points to the boundary.
     * @param {boolean} cube - If true, the box is a cube.
     * @returns {{lower: number[], upper: number[]}}
     */
    let lower = [0, 1, 2].map(m => Math.min(...points.map(p => p[m])));
    let upper = [0, 1, 2].map(m => Math.max(...points.map(p => p[m])));

    if (cube) {
        const size = Math.max(...[0, 1, 2].map(m => upper[m] - lower[m]));
        const mid = [0, 1, 2].map(m => (lower[m] + upper[m]) / 2);
        lower = mid.map(c => c - size / 2 - dist);
        upper = mid.map(c => c + size / 2 + dist);
    } else {
        lower = lower.map(c => c - dist);
        upper = upper.map(c => c + dist);
    }

    // Avoid empty intervals (all points in a plane and dist = 0)
    for (let m = 0; m < 3; m++) {
        if (upper[m] - lower[m] < 1e-9) {
            lower[m] -= 1;
            upper[m] += 1;
        }
    }

    // Round outwards to two significant digits
    const step = Math.pow(10, Math.floor(Math.log10(Math.max(...[0, 1, 2].map(m => upper[m] - lower[m])))) - 1);
    const clean = v => Number(v.toFixed(10)) + 0;
    lower = lower.map(v => clean(Math.floor(v / step + 1e-9) * step));
    upper = upper.map(v => clean(Math.ceil(v / step - 1e-9) * step));

    return { lower, upper };
}

/*
Everything together
*/

function computeTropicalSurface(monomials, options) {
    /**
     * @param {Array<[number, number, number, number]>} monomials - Output of parsePolynomialTrop.
     * @param {Object} options - {useMin, autoAdjust, autoCube, autoBoundaryDist, lower, upper, eps}
     *        lower, upper: the manual plotting range (ignored if autoAdjust).
     * @returns {Object} The output of computeSurfaceInBox, in the convention of the user, together with
     *          lower, upper (the plotting range), vertices (see minimalCellPoints), newtonDim,
     *          verticesOutside (number of vertices outside the plotting range) and the monomials.
     */
    const eps = options.eps ?? 1e-9;
    const sign = options.useMin ? -1 : 1; // min(c + <a,p>) = -max(-c + <a,-p>)
    const exps = monomials.map(m => [m[0], m[1], m[2]]);
    const coeffs = monomials.map(m => sign * m[3]);

    const { points, newtonDim } = minimalCellPoints(exps, coeffs, eps);

    let lower, upper; // in the max-convention
    if (options.autoAdjust) {
        ({ lower, upper } = chooseRange(points.length > 0 ? points : [[0, 0, 0]], options.autoBoundaryDist, options.autoCube));
    } else if (options.useMin) {
        lower = options.upper.map(v => -v);
        upper = options.lower.map(v => -v);
    } else {
        lower = options.lower.slice();
        upper = options.upper.slice();
    }

    const surface = computeSurfaceInBox(exps, coeffs, lower, upper, eps);

    const verticesOutside = points.filter(p => [0, 1, 2].some(m => p[m] < lower[m] || p[m] > upper[m])).length;

    // Back to the convention of the user
    const flip = p => p.map(t => sign * t + 0);
    if (options.useMin) {
        [lower, upper] = [upper.map(v => -v + 0), lower.map(v => -v + 0)];
        for (const cell of surface.cells) {
            cell.vertices = cell.vertices.map(flip);
        }
        for (const edge of surface.edges) {
            edge.p = flip(edge.p);
            edge.q = flip(edge.q);
        }
        for (const point of surface.lonePoints) {
            point.p = flip(point.p);
        }
    }

    return {
        ...surface,
        lower,
        upper,
        vertices: points.map(flip),
        newtonDim,
        verticesOutside,
        monomials
    };
}

// For tests with Node.js
if (typeof module !== "undefined" && module.exports) {
    module.exports = {
        parsePolynomialTrop,
        computeTropicalSurface,
        computeSurfaceInBox,
        minimalCellPoints,
        chooseRange,
        affineDim
    };
}
