/*
Plotting the intersection of two tropical curves V(F) (curve 1, blue) and V(G) (curve 2, red).

Idea: The tropical curve of the product F*G is the union of V(F) and V(G), and its dual
subdivision is the mixed subdivision of the Newton polytopes of F and G. At every point
the monomials of F*G attaining the maximum are exactly the sums of a maximal monomial
of F and a maximal monomial of G. Hence:
- The vertices and edges of the union are the vertices and edges of V(F*G).
- Whether an edge belongs to V(F) or V(G) is decided by evaluating F and G on it.
- At a transversal intersection point the maximal exponents of F lie on an edge of the
  dual subdivision of F, those of G on an edge of the one of G, and these edges are not
  parallel. The dual cell is the parallelogram spanned by them, its area is the
  multiplicity of the intersection point.
*/

/*
Eventlisteners
*/

const plotForm = document.getElementById("plotForm");

plotForm.addEventListener("submit", function(event) {
    event.preventDefault();
    plotte();
});

const copyLatexButton = document.getElementById("copyLatexButton");
copyLatexButton.addEventListener("click", copyLatex);

const copyLatexButtonDS = document.getElementById("copyLatexButtonDS");
copyLatexButtonDS.addEventListener("click", copyLatexDS);

const latexButton = document.getElementById("latexButton");
latexButton.addEventListener("click", showLatex);

const latexButtonDS = document.getElementById("latexButtonDS");
latexButtonDS.addEventListener("click", showLatexDS);

const closeLatexButton = document.getElementById("closeLatexButton");
closeLatexButton.addEventListener("click", closeLatex);

const closeLatexButtonDS = document.getElementById("closeLatexButtonDS");
closeLatexButtonDS.addEventListener("click", closeLatexDS);

addShiftListeners("polyF", "shiftXF", "shiftYF", "shiftButtonF");
addShiftListeners("polyG", "shiftXG", "shiftYG", "shiftButtonG");

// Changing one of the checkbox options plots, if both polynomials are entered
for (const checkbox of plotForm.querySelectorAll('input[type="checkbox"]')) {
    checkbox.addEventListener("change", plotIfComplete);
}

/*
The initial LaTeX code: Color declarations for easier changes
*/
let latexCode = "";
let latexCodeDS = "";

/*
eps - the tolerance how close to numbers must be to be considered equal.
Larger than in the plotter for single curves, since the product of the polynomials
and solving linear systems accumulate rounding errors.
*/
const eps = 1e-9;

/*
Colors and TikZ color macros for the two curves, and for edges where the curves overlap
(only possible if the intersection is not transversal).
The fill colors are blue!25 and red!25.
*/
const curveStyles = {
    one: { color: "blue", fill: "#bfbfff", tikzCurve: "\\colorCurveOne", tikzDS: "\\colorDSOne" },
    two: { color: "red", fill: "#ffbfbf", tikzCurve: "\\colorCurveTwo", tikzDS: "\\colorDSTwo" },
    both: { color: "#800080", fill: null, tikzCurve: "\\colorOverlap", tikzDS: "\\colorDSOverlap" }
};

const curveWidth = 0.03; // stroke width of the curves, in SVG units (the plot is 20 units wide)
const labelSize = 0.5; // font size of weights and multiplicities, in SVG units
const multLabelDist = 0.35; // distance of the multiplicities to their intersection point, in SVG units
const circleRadius = 0.3; // radius of the circles around non-transversal points, in SVG units
const circleRadiusTikz = "4pt"; // the same for the TikZ code

/*
Computations & helpers for that (taken from the plotter for single curves)
*/

function parseTermsTrop(polyStr) {
    /*
     * Splits a tropical polynomial in the variables x and y into its terms,
     * in the order of the input and without combining terms with equal exponents.
     *
     * Example input:
     *   "3.2x2y + yx + 9 + 2xy"
     *
     * Returns:
     *   Array of triples [xExponent, yExponent, coefficient]
     *
     * Example output:
     * [
     *   [2, 1, 3.2],
     *   [1, 1, 0],
     *   [0, 0, 9],
     *   [1, 1, 2]
     * ]
     */

    // Remove whitespace and symbols that should be ignored.
    polyStr = polyStr
        .replaceAll(" ", "")
        .replaceAll("*", "")
        .replaceAll("^", "")
        .replaceAll("(", "")
        .replaceAll(")", "")
        .replaceAll("+-", "-");

    // Split the polynomial into individual terms.
    const terms = polyStr.match(/[+-]?\d*\.?\d*[xy]?\d*[xy]?\d*/g) || [];

    const result = [];

    for (const term of terms) {
        if (!term) {
            continue;
        }

        const match = term.match(
            /^([+-]?\d*\.?\d*)([xy]\d*)?([xy]\d*)?$/
        );

        if (!match) {
            throw new Error(`Invalid term: ${term}`);
        }

        const [, coeffStr, var1, var2] = match;

        let coefficient;

        if (coeffStr === "" || coeffStr === "+" || coeffStr === "-") {
            coefficient = Number(coeffStr + "0");
        } else {
            coefficient = parseFloat(coeffStr);
        }

        let xExp = 0;
        let yExp = 0;

        for (const variable of [var1, var2]) {
            if (!variable) {
                continue;
            }

            if (variable[0] === "x") {
                xExp = variable.length > 1
                    ? parseInt(variable.slice(1), 10)
                    : 1;
            } else if (variable[0] === "y") {
                yExp = variable.length > 1
                    ? parseInt(variable.slice(1), 10)
                    : 1;
            }
        }

        result.push([xExp, yExp, coefficient]);
    }

    return result;
}

function parsePolynomialTrop(polyStr, useMin) {
    /*
     * Parses a tropical polynomial in the variables x and y.
     * Terms with equal exponents are combined (taking the maximum or minimum).
     *
     * Example input:
     *   "3.2x2y + yx + 9"
     *
     * Returns:
     *   Array of triples [xExponent, yExponent, coefficient]
     *
     * Example output:
     * [
     *   [2, 1, 3.2],
     *   [1, 1, 0],
     *   [0, 0, 9]
     * ]
     */

    // Temporary map:
    // "xExp,yExp" -> coefficient
    const coefficients = new Map();

    for (const [xExp, yExp, coefficient] of parseTermsTrop(polyStr)) {
        const key = `${xExp},${yExp}`;

        if (coefficients.has(key)) {
            const oldCoefficient = coefficients.get(key);

            coefficients.set(
                key,
                useMin
                    ? Math.min(coefficient, oldCoefficient)
                    : Math.max(coefficient, oldCoefficient)
            );
        } else {
            coefficients.set(key, coefficient);
        }
    }

    // Convert to list of triples.
    const monomials = [];

    for (const [key, coefficient] of coefficients) {
        const [xExp, yExp] = key.split(",").map(Number);
        monomials.push([xExp, yExp, coefficient]);
    }

    // Sort descending by x exponent, then by y exponent.
    monomials.sort((a, b) => {
        if (a[0] !== b[0]) {
            return b[0] - a[0];
        }
        return b[1] - a[1];
    });

    return monomials;
}

function formatPolynomialTrop(terms) {
    /*
     * Writes a tropical polynomial as a string that parseTermsTrop reads back.
     *
     * Example input:
     *   [[2, 0, -2], [1, 1, 0], [0, 0, 1.5]]
     *
     * Example output:
     *   "-2x2+xy+1.5"
     *
     * A coefficient 0 is omitted in front of non-constant monomials, as the parser reads it as 0.
     */
    let result = "";

    for (const [xExp, yExp, coefficient] of terms) {
        const variables =
            (xExp > 0 ? "x" + (xExp > 1 ? xExp : "") : "") +
            (yExp > 0 ? "y" + (yExp > 1 ? yExp : "") : "");

        // Round away floating point noise like 0.30000000000000004, without exponential notation
        let coeffStr = coefficient.toFixed(10).replace(/\.?0+$/, "");

        if (coeffStr === "-0") {
            coeffStr = "0";
        }
        if (coeffStr === "0" && variables !== "") {
            coeffStr = "";
        }

        // Terms are separated by their signs
        const sign = (result !== "" && !coeffStr.startsWith("-")) ? "+" : "";
        result = result + sign + coeffStr + variables;
    }

    return result;
}

function changeSignsOfPolynomial(monomials) {
    /*
     * Input:
     *   Array of triples [xExponent, yExponent, coefficient]
     *
     * Output:
     *   A new array with all coefficients multiplied by -1.
     *
     * The input array is not modified.
     */

    const result = [];

    for (const [xExp, yExp, coeff] of monomials) {
        result.push([xExp, yExp, -coeff]);
    }

    return result;
}

function evaluateTropMonomials(x, y, A) {
    /**
     * Evaluates all tropical monomials of a tropical polynomial at the point (x, y).
     *
     * The polynomial is stored as an array of triples:
     *   [xExponent, yExponent, coefficient]
     *
     * @param {number} x - x-coordinate.
     * @param {number} y - y-coordinate.
     * @param {Array<[number, number, number]>} A - List of monomials.
     * @returns {number[]} An array containing the value of each tropical monomial.
     */

    const values = [];

    for (const [xExp, yExp, coeff] of A) {
        values.push(xExp * x + yExp * y + coeff);
    }

    return values;
}

function evaluateTropPolynomial(x,y,A){
    /**
     * Evaluates a tropical polynomial at the point (x, y).
     *
     * The polynomial is stored as an array of triples:
     *   [xExponent, yExponent, coefficient]
     *
     * @param {number} x - x-coordinate.
     * @param {number} y - y-coordinate.
     * @param {Array<[number, number, number]>} A - List of monomials.
     * @returns {number} the value of the tropical polynomial at (x,y).
     */
    return Math.max(...evaluateTropMonomials(x,y,A));
}

function tupleIsInList(x, L, eps) {
    /**
     * Checks whether a given 2-dimensional point occurs in a list of points,
     * up to a given tolerance in the infinity norm.
     *
     * @param {number[]} x - A point represented as [x, y].
     * @param {number[][]} L - An array of points, each represented as [x, y].
     * @param {number} eps - Numerical tolerance.
     * @returns {number} The index of the matching point, or -1 if no match is found.
     */

    // Iterate over all points in the list.
    for (let i = 0; i < L.length; i++) {

        // Check whether both coordinates agree up to the given tolerance.
        if (
            Math.abs(L[i][0] - x[0]) < eps &&
            Math.abs(L[i][1] - x[1]) < eps
        ) {
            return i;
        }
    }

    // No matching point was found.
    return -1;
}

function solve2x2(A, b) {
    /**
     * Solves A*x = b for an invertible 2x2 matrix A.
     *
     * @param {number[][]} A - [[a11, a12], [a21, a22]]
     * @param {number[]} b - [b1, b2]
     * @returns {number[] | null} [x, y], or null if the system has no unique solution.
     */
    const [[a11, a12], [a21, a22]] = A;
    const [b1, b2] = b;

    const det = a11 * a22 - a12 * a21;

    if (det === 0) {
        return null;
    }

    return [
        (b1 * a22 - a12 * b2) / det,
        (a11 * b2 - b1 * a21) / det
    ];
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

function showError(errorText){
    const error = document.getElementById("error");
    error.textContent = errorText;
    return;
}

function xSVGCoords(x,aX,bX){
    return (-10+(x-aX)/(bX-aX) * 20);
}

function ySVGCoords(y,aY,bY){
    return (10 - (y - aY)/(bY-aY) * 20);
}

function plotFrame(aX, bX, aY, bY) {
    /**
     * Places the plotting range [aX, bX] x [aY, bY] in the SVG (20 x 20 units) without distortion:
     * Both axes get the same scale, the longer side of the range fills the 20 units,
     * and the range is centered.
     *
     * @param {number} aX, bX, aY, bY - The plotting range.
     * @returns {Object} {left, right, top, bottom, toSVG}: the bounds of the range in SVG
     *          coordinates and the map from points [x, y] to SVG coordinates.
     */
    const scale = 20 / Math.max(bX - aX, bY - aY);
    const halfWidth = scale * (bX - aX) / 2;
    const halfHeight = scale * (bY - aY) / 2;

    return {
        left: -halfWidth,
        right: halfWidth,
        top: -halfHeight,
        bottom: halfHeight,
        toSVG: ([x, y]) => [-halfWidth + (x - aX) * scale, halfHeight - (y - aY) * scale]
    };
}

/*
Computations for the intersection
*/

function multiplyTropPolynomials(F, G) {
    /**
     * Computes the tropical product of two tropical polynomials: Exponents and
     * coefficients of all pairs of monomials are added. Of the monomials with
     * equal exponents only the one with maximal coefficient is kept, since the
     * others are never maximal.
     *
     * @param {Array<[number, number, number]>} F - List of monomials.
     * @param {Array<[number, number, number]>} G - List of monomials.
     * @returns {Array<[number, number, number]>} The monomials of F*G.
     */
    const coefficients = new Map();

    for (const [i1, j1, c1] of F) {
        for (const [i2, j2, c2] of G) {
            const key = `${i1 + i2},${j1 + j2}`;

            if (!coefficients.has(key) || coefficients.get(key) < c1 + c2) {
                coefficients.set(key, c1 + c2);
            }
        }
    }

    const monomials = [];

    for (const [key, coefficient] of coefficients) {
        const [xExp, yExp] = key.split(",").map(Number);
        monomials.push([xExp, yExp, coefficient]);
    }

    return monomials;
}

function maxIndices(x, y, A) {
    /**
     * Returns the indices of the monomials of A attaining the maximum at (x, y),
     * up to the tolerance eps.
     *
     * @param {number} x - x-coordinate.
     * @param {number} y - y-coordinate.
     * @param {Array<[number, number, number]>} A - List of monomials.
     * @returns {number[]} The indices of the maximal monomials.
     */
    const values = evaluateTropMonomials(x, y, A);
    const max = Math.max(...values);
    const indices = [];

    for (let i = 0; i < values.length; i++) {
        if (max - values[i] < eps) {
            indices.push(i);
        }
    }

    return indices;
}

function computeVertices(P) {
    /**
     * Computes the vertices of the tropical curve of P, i.e. the points where
     * three monomials with affinely independent exponents attain the maximum.
     *
     * @param {Array<[number, number, number]>} P - List of monomials.
     * @returns {number[][]} The vertices as [x, y], each listed once.
     */
    const vertices = [];

    // Walk through the monomials of P, check all triples
    for (let a = 0; a < P.length; a++) {
        for (let b = a + 1; b < P.length; b++) {
            for (let c = b + 1; c < P.length; c++) {
                const [i1, j1, c1] = P[a];
                const [i2, j2, c2] = P[b];
                const [i3, j3, c3] = P[c];

                // The three monomials agree iff (i1-i2)x + (j1-j2)y = c2-c1 and (i1-i3)x + (j1-j3)y = c3-c1.
                // This has no unique solution iff the exponents are affinely dependent.
                const pt = solve2x2([[i1 - i2, j1 - j2], [i1 - i3, j1 - j3]], [c2 - c1, c3 - c1]);

                if (pt === null) {
                    continue;
                }

                // Check if it is actually a maximum
                if (evaluateTropPolynomial(pt[0], pt[1], P) - (c1 + i1 * pt[0] + j1 * pt[1]) < eps &&
                    tupleIsInList(pt, vertices, eps) === -1) {
                    vertices.push(pt);
                }
            }
        }
    }

    return vertices;
}

function computeLinePoints(P) {
    /**
     * For a tropical curve without vertices (which consists of parallel lines):
     * Computes points on the lines, used to determine the plotting range.
     *
     * @param {Array<[number, number, number]>} P - List of monomials.
     * @returns {number[][]} Points [x, y] on the curve.
     */
    const points = [];

    for (let a = 0; a < P.length; a++) {
        for (let b = a + 1; b < P.length; b++) {
            const [i1, j1, c1] = P[a];
            const [i2, j2, c2] = P[b];

            // The monomials agree on the line (i1-i2)x + (j1-j2)y = c2-c1.
            // Take its point with y = 0, or with x = 0 if the line is horizontal.
            const pt = (i1 !== i2)
                ? [(c2 - c1) / (i1 - i2), 0]
                : [0, (c2 - c1) / (j1 - j2)];

            // If it is a point of the curve
            if (evaluateTropPolynomial(pt[0], pt[1], P) - (c1 + i1 * pt[0] + j1 * pt[1]) < eps) {
                points.push(pt);
            }
        }
    }

    return points;
}

function computeBoundaryPoints(P, aX, bX, aY, bY) {
    /**
     * Computes the points where the tropical curve of P meets the boundary of
     * the plotting range [aX, bX] x [aY, bY].
     *
     * @param {Array<[number, number, number]>} P - List of monomials.
     * @param {number} aX, bX, aY, bY - The plotting range.
     * @returns {number[][]} The points as [x, y], each listed once.
     */
    const points = [];

    // Check all pairs of monomials
    for (let a = 0; a < P.length; a++) {
        for (let b = a + 1; b < P.length; b++) {
            const [i1, j1, c1] = P[a];
            const [i2, j2, c2] = P[b];

            // The monomials agree on the line (i1-i2)x + (j1-j2)y = c2-c1. Intersect it with the boundary.
            const candidates = [];

            if (j1 !== j2) {
                for (const x of [aX, bX]) {
                    const y = (c2 - c1 - (i1 - i2) * x) / (j1 - j2);
                    if (y >= aY - eps && y <= bY + eps) {
                        candidates.push([x, y]);
                    }
                }
            }

            if (i1 !== i2) {
                for (const y of [aY, bY]) {
                    const x = (c2 - c1 - (j1 - j2) * y) / (i1 - i2);
                    if (x >= aX - eps && x <= bX + eps) {
                        candidates.push([x, y]);
                    }
                }
            }

            // Keep the points where the two monomials are maximal
            for (const pt of candidates) {
                if (evaluateTropPolynomial(pt[0], pt[1], P) - (c1 + i1 * pt[0] + j1 * pt[1]) < eps &&
                    tupleIsInList(pt, points, eps) === -1) {
                    points.push(pt);
                }
            }
        }
    }

    return points;
}

function computeEdges(points, maxSets, P, F, G) {
    /**
     * Computes the edges of the union of the curves: Two points are connected
     * iff they share at least two maximal monomials of F*G.
     *
     * @param {number[][]} points - Vertices and points on the boundary, as [x, y].
     * @param {number[][]} maxSets - For each point the indices of the maximal monomials of P.
     * @param {Array<[number, number, number]>} P - The monomials of F*G.
     * @param {Array<[number, number, number]>} F, G - The monomials of F and G.
     * @returns {Object[]} Edges {start, end, curve, weight, dual}: start and end are indices
     *          in points, curve is "one", "two" or "both" (if the curves overlap there),
     *          dual are the end points of the dual edge in the dual subdivision.
     */
    const edges = [];

    for (let s = 0; s < points.length; s++) {
        for (let e = s + 1; e < points.length; e++) {
            const common = maxSets[s].filter(k => maxSets[e].includes(k));

            if (common.length < 2) {
                continue;
            }

            // The curve(s) the edge belongs to: Evaluate F and G in its middle
            const xMid = (points[s][0] + points[e][0]) / 2;
            const yMid = (points[s][1] + points[e][1]) / 2;
            const inF = maxIndices(xMid, yMid, F).length > 1;
            const inG = maxIndices(xMid, yMid, G).length > 1;

            // The dual edge connects the common exponents furthest apart, its lattice length is the weight
            let weight = 0;
            let dual = null;

            for (let a = 0; a < common.length; a++) {
                for (let b = a + 1; b < common.length; b++) {
                    const [i1, j1] = P[common[a]];
                    const [i2, j2] = P[common[b]];
                    const w = gcd(i1 - i2, j1 - j2);

                    if (w > weight) {
                        weight = w;
                        dual = [[i1, j1], [i2, j2]];
                    }
                }
            }

            edges.push({
                start: s,
                end: e,
                curve: (inF && inG) ? "both" : (inF ? "one" : "two"),
                weight: weight,
                dual: dual
            });
        }
    }

    return edges;
}

function collinearExtent(points) {
    /**
     * Checks whether integer points lie on a line.
     *
     * @param {number[][]} points - At least two points [i, j] with integer coordinates.
     * @returns {number[] | null} The vector between the two outermost points if all
     *          points lie on a line, otherwise null.
     */
    let p = points[0];
    let q = points[1];
    let maxDist = -1;

    for (let a = 0; a < points.length; a++) {
        for (let b = a + 1; b < points.length; b++) {
            const dist = (points[a][0] - points[b][0]) ** 2 + (points[a][1] - points[b][1]) ** 2;

            if (dist > maxDist) {
                maxDist = dist;
                p = points[a];
                q = points[b];
            }
        }
    }

    const dx = q[0] - p[0];
    const dy = q[1] - p[1];

    for (const r of points) {
        if (dx * (r[1] - p[1]) - dy * (r[0] - p[0]) !== 0) {
            return null;
        }
    }

    return [dx, dy];
}

function convexHull(points) {
    /**
     * Computes the convex hull of integer points (monotone chain algorithm).
     *
     * @param {number[][]} points - Points [i, j] with integer coordinates.
     * @returns {number[][]} The vertices of the convex hull in counterclockwise order.
     */
    const sorted = [...points].sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]));
    const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);

    const lower = [];
    for (const pt of sorted) {
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], pt) <= 0) {
            lower.pop();
        }
        lower.push(pt);
    }

    const upper = [];
    for (const pt of [...sorted].reverse()) {
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], pt) <= 0) {
            upper.pop();
        }
        upper.push(pt);
    }

    // The last point of each chain is the first point of the other one
    return lower.slice(0, -1).concat(upper.slice(0, -1));
}

function doubleArea(points) {
    // Twice the area of the convex hull of integer points (an integer)
    const hull = convexHull(points);
    let sum = 0;

    for (let k = 0; k < hull.length; k++) {
        const [i1, j1] = hull[k];
        const [i2, j2] = hull[(k + 1) % hull.length];
        sum = sum + i1 * j2 - i2 * j1;
    }

    return Math.abs(sum);
}

function mixedArea(A, B) {
    /**
     * Computes the mixed area MV(conv A, conv B) = Area(conv A + conv B) - Area(conv A) - Area(conv B).
     * For the dual cells of two curves at a common point, this is the stable intersection
     * multiplicity of the point (and the usual multiplicity if the intersection is transversal).
     *
     * @param {number[][]} A, B - Points [i, j] with integer coordinates.
     * @returns {number} The mixed area.
     */
    const sums = [];

    for (const [i1, j1] of A) {
        for (const [i2, j2] of B) {
            sums.push([i1 + i2, j1 + j2]);
        }
    }

    return (doubleArea(sums) - doubleArea(A) - doubleArea(B)) / 2;
}

function classifyVertices(vertices, maxSets, P, F, G) {
    /**
     * Classifies the vertices of the union of the curves by evaluating F and G there.
     *
     * @param {number[][]} vertices - The vertices as [x, y].
     * @param {number[][]} maxSets - For each vertex the indices of the maximal monomials of P.
     * @param {Array<[number, number, number]>} P - The monomials of F*G.
     * @param {Array<[number, number, number]>} F, G - The monomials of F and G.
     * @returns {Object[]} For each vertex {type, multiplicity, cell}: type is
     *          "one" / "two" (vertex of curve 1 / curve 2 only),
     *          "intersection" (transversal intersection point, with its multiplicity) or
     *          "nonTransversal" (vertex of one curve lying on the other curve, with its
     *          stable intersection multiplicity).
     *          cell are the vertices of the dual cell in the dual subdivision.
     */
    return vertices.map(([x, y], index) => {
        const A = maxIndices(x, y, F).map(k => [F[k][0], F[k][1]]);
        const B = maxIndices(x, y, G).map(k => [G[k][0], G[k][1]]);
        const cell = convexHull(maxSets[index].map(k => [P[k][0], P[k][1]]));

        if (A.length > 1 && B.length > 1) {
            // The point lies on both curves. The intersection is transversal iff the maximal exponents
            // of F and of G each lie on a line (i.e. the point lies on an edge of each curve)
            // and these lines are not parallel.
            const edgeA = collinearExtent(A);
            const edgeB = collinearExtent(B);

            if (edgeA !== null && edgeB !== null) {
                const det = edgeA[0] * edgeB[1] - edgeA[1] * edgeB[0];

                if (det !== 0) {
                    // The dual cell is the parallelogram spanned by the two dual edges, its area is the multiplicity
                    return { type: "intersection", multiplicity: Math.abs(det), cell: cell };
                }
            }

            // The stable intersection multiplicity is the mixed area of the dual cells of the two curves
            return { type: "nonTransversal", multiplicity: mixedArea(A, B), cell: cell };
        }

        return { type: (A.length > 1) ? "one" : "two", multiplicity: null, cell: cell };
    });
}

function reportNonTransversality(points, edges, vertexInfo) {
    /**
     * Shows an error message if the intersection is not transversal.
     *
     * @param {number[][]} points - Vertices (first) and points on the boundary, as [x, y].
     * @param {Object[]} edges - The edges, see computeEdges.
     * @param {Object[]} vertexInfo - The classification of the vertices, see classifyVertices.
     */
    const badPoints = points.filter((pt, index) => index < vertexInfo.length && vertexInfo[index].type === "nonTransversal");
    const overlap = edges.some(edge => edge.curve === "both");

    if (badPoints.length === 0 && !overlap) {
        return;
    }

    const problems = [];

    if (badPoints.length > 0) {
        const list = badPoints.map(([x, y]) => `(${Number(x.toFixed(4))}, ${Number(y.toFixed(4))})`).join(", ");
        problems.push(badPoints.length === 1
            ? `a vertex of one curve lies on the other curve at ${list} (circled in the plot)`
            : `vertices of one curve lie on the other curve at ${list} (circled in the plot)`);
    }

    if (overlap) {
        problems.push("the curves overlap along the purple edges");
    }

    const note = (badPoints.length > 0)
        ? " The numbers at the circled points are their stable intersection multiplicities."
        : "";

    showError(`The intersection is not transversal: ${problems.join(", and ")}.${note}`);
}

function freeDirection(angles) {
    /**
     * Returns the direction bisecting the largest angle between the given directions.
     *
     * @param {number[]} angles - Directions as angles in radians.
     * @returns {number} The bisecting direction in radians.
     */
    if (angles.length === 0) {
        return Math.PI / 4;
    }

    const sorted = [...angles].sort((a, b) => a - b);
    let bestGap = -1;
    let bestDirection = 0;

    for (let k = 0; k < sorted.length; k++) {
        const next = (k + 1 < sorted.length) ? sorted[k + 1] : sorted[0] + 2 * Math.PI;
        const gap = next - sorted[k];

        if (gap > bestGap) {
            bestGap = gap;
            bestDirection = sorted[k] + gap / 2;
        }
    }

    return bestDirection;
}

function cleanNumber(v) {
    // Rounds away floating point noise like 4.000000000000001 for the TikZ code
    return Number(v.toFixed(10));
}

/*
Coordinate lines in the plot of the curves (taken from the plotter for single curves)
*/

// The number of intervals between coordinate lines the plotting range is divided into (approximately)
const gridTargetIntervals = 8;

function niceGridStep(a, b) {
    /**
     * Computes a "nice" distance between coordinate lines for the interval [a, b]:
     * 1, 2 or 5 times a power of ten, such that [a, b] is divided into roughly
     * gridTargetIntervals pieces.
     *
     * @param {number} a - Lower bound of the interval.
     * @param {number} b - Upper bound of the interval.
     * @returns {number} The distance, or NaN if the interval is degenerate.
     */
    const range = b - a;

    if (!(range > 0) || !Number.isFinite(range)) {
        return NaN;
    }

    const power = Math.pow(10, Math.floor(Math.log10(range / gridTargetIntervals)));

    let bestStep = NaN;
    let bestRatio = Infinity;

    for (const factor of [1, 2, 5, 10]) {
        const step = factor * power;
        const intervals = range / step;
        // How far the number of intervals is off from the target, as a factor >= 1
        const ratio = Math.max(intervals / gridTargetIntervals, gridTargetIntervals / intervals);

        if (ratio < bestRatio) {
            bestRatio = ratio;
            bestStep = step;
        }
    }

    return bestStep;
}

function gridValues(a, b, step) {
    /**
     * Returns all integer multiples of step in the interval [a, b].
     * Since 0 is a multiple of every step, the origin always lies on
     * coordinate lines if it is in the plotting range.
     *
     * @param {number} a - Lower bound of the interval.
     * @param {number} b - Upper bound of the interval.
     * @param {number} step - Distance between coordinate lines.
     * @returns {number[]} The values, rounded to remove floating point noise.
     */
    const tolerance = 1e-9;
    // Number of decimal places needed to represent multiples of step exactly
    const decimals = Math.max(0, -Math.floor(Math.log10(step) + tolerance));
    const values = [];

    for (let k = Math.ceil(a / step - tolerance); k <= Math.floor(b / step + tolerance); k++) {
        values.push(Number((k * step).toFixed(decimals)));
    }

    return values;
}

function createSVGElement(tag, attributes) {
    // Creates an SVG element with the given attributes
    const element = document.createElementNS("http://www.w3.org/2000/svg", tag);

    for (const [name, value] of Object.entries(attributes)) {
        element.setAttribute(name, value);
    }

    return element;
}

function drawCoordinateLines(svg, aX, bX, aY, bY, showLabels) {
    /**
     * Draws coordinate lines (and optionally their labels) into the plot
     * of the curves and appends the corresponding TikZ code to latexCode.
     * Since the plot is not distorted, both directions get the same distance,
     * chosen according to the longer side of the plotting range.
     *
     * @param {SVGElement} svg - The plot of the curves.
     * @param {number} aX, bX, aY, bY - The plotting range.
     * @param {boolean} showLabels - Whether to label the lines with their values.
     */
    const frame = plotFrame(aX, bX, aY, bY);
    const step = niceGridStep(0, Math.max(bX - aX, bY - aY));

    if (Number.isNaN(step)) {
        return;
    }

    const xValues = gridValues(aX, bX, step);
    const yValues = gridValues(aY, bY, step);

    // Rounded plotting range for the TikZ code
    const [aXTikz, bXTikz, aYTikz, bYTikz] = [aX, bX, aY, bY].map(v => Number(v.toFixed(10)));

    const gridColor = "#d0d0d0";
    const labelColor = "#808080";
    const labelSize = 0.4; // in SVG units, the plot is 20 units wide
    const labelGap = 0.1; // distance of labels to their line and to the boundary

    latexCode = latexCode + "\n\\def\\colorGrid{black!20}";
    if (showLabels) {
        latexCode = latexCode + "\n\\def\\colorGridLabels{black!60}";
    }

    // Vertical lines, labels at the bottom (right of the line, left of it close to the right boundary)
    for (const x of xValues) {
        const xSVG = frame.toSVG([x, aY])[0];

        svg.appendChild(createSVGElement("line", {
            "x1": xSVG, "y1": frame.top, "x2": xSVG, "y2": frame.bottom,
            "stroke": gridColor, "stroke-width": 0.02
        }));

        latexCode = latexCode + `\n\\draw[\\colorGrid, very thin] (${x},${aYTikz}) -- (${x},${bYTikz});`;

        if (showLabels) {
            const nearRightBoundary = xSVG > frame.right - 3 * labelSize;
            const label = createSVGElement("text", {
                "x": nearRightBoundary ? xSVG - labelGap : xSVG + labelGap,
                "y": frame.bottom - 1.5 * labelGap,
                "text-anchor": nearRightBoundary ? "end" : "start",
                "font-size": labelSize,
                "fill": labelColor
            });
            label.textContent = `${x}`;
            svg.appendChild(label);

            latexCode = latexCode + `\n\\node[\\colorGridLabels, below, font=\\scriptsize] at (${x},${aYTikz}) {$${x}$};`;
        }
    }

    // Horizontal lines, labels at the left (above the line, below it close to the upper boundary)
    for (const y of yValues) {
        const ySVG = frame.toSVG([aX, y])[1];

        svg.appendChild(createSVGElement("line", {
            "x1": frame.left, "y1": ySVG, "x2": frame.right, "y2": ySVG,
            "stroke": gridColor, "stroke-width": 0.02
        }));

        latexCode = latexCode + `\n\\draw[\\colorGrid, very thin] (${aXTikz},${y}) -- (${bXTikz},${y});`;

        if (showLabels) {
            // In the SVG, skip labels that would collide with the labels of the vertical lines at the bottom
            if (ySVG < frame.bottom - 2 * labelSize) {
                const nearUpperBoundary = ySVG - labelGap - labelSize < frame.top;
                const label = createSVGElement("text", {
                    "x": frame.left + 1.5 * labelGap,
                    "y": nearUpperBoundary ? ySVG + labelGap + 0.75 * labelSize : ySVG - labelGap,
                    "text-anchor": "start",
                    "font-size": labelSize,
                    "fill": labelColor
                });
                label.textContent = `${y}`;
                svg.appendChild(label);
            }

            // In TikZ, the labels are outside the plotting range and do not collide
            latexCode = latexCode + `\n\\node[\\colorGridLabels, left, font=\\scriptsize] at (${aXTikz},${y}) {$${y}$};`;
        }
    }
}

/*
Drawing
*/

function drawCurves(svg, points, edges, vertexInfo, aX, bX, aY, bY) {
    /**
     * Draws the two curves with their weights, the multiplicities of the intersection
     * points and circles around points where the intersection is not transversal
     * (with their stable intersection multiplicities).
     * Appends the corresponding TikZ code to latexCode.
     *
     * @param {SVGElement} svg - The plot of the curves.
     * @param {number[][]} points - Vertices (first) and points on the boundary, as [x, y].
     * @param {Object[]} edges - The edges, see computeEdges.
     * @param {Object[]} vertexInfo - The classification of the vertices, see classifyVertices.
     * @param {number} aX, bX, aY, bY - The plotting range.
     */
    const toSVG = plotFrame(aX, bX, aY, bY).toSVG;

    if (edges.some(edge => edge.curve === "both")) {
        latexCode = latexCode + "\n\\def\\colorOverlap{violet}";
    }
    if (vertexInfo.some(info => info.type === "nonTransversal")) {
        latexCode = latexCode + "\n\\def\\colorNonTransversal{black}";
    }

    // Edges and their weights
    for (const edge of edges) {
        const style = curveStyles[edge.curve];
        const [x1, y1] = points[edge.start];
        const [x2, y2] = points[edge.end];
        const [x1SVG, y1SVG] = toSVG(points[edge.start]);
        const [x2SVG, y2SVG] = toSVG(points[edge.end]);

        svg.appendChild(createSVGElement("line", {
            "x1": x1SVG, "y1": y1SVG, "x2": x2SVG, "y2": y2SVG,
            "stroke": style.color, "stroke-width": curveWidth
        }));

        latexCode = latexCode + `\n\\draw[${style.tikzCurve}] (${cleanNumber(x1)},${cleanNumber(y1)}) -- (${cleanNumber(x2)},${cleanNumber(y2)});`;

        if (edge.weight > 1) {
            const xMid = (x1 + x2) / 2;
            const yMid = (y1 + y2) / 2;
            const [xMidSVG, yMidSVG] = toSVG([xMid, yMid]);

            const weight = createSVGElement("text", {
                "x": xMidSVG, "y": yMidSVG,
                "font-size": labelSize, "fill": style.color
            });
            weight.textContent = `${edge.weight}`;
            svg.appendChild(weight);

            latexCode = latexCode + `\n\\node[${style.tikzCurve}] at (${cleanNumber(xMid)},${cleanNumber(yMid)}) {${edge.weight}};`;
        }
    }

    // Multiplicities of the intersection points, circles around non-transversal points
    vertexInfo.forEach((info, index) => {
        if (info.type !== "intersection" && info.type !== "nonTransversal") {
            return;
        }

        const [x, y] = points[index];
        const [xSVG, ySVG] = toSVG(points[index]);
        const circled = (info.type === "nonTransversal");

        if (circled) {
            svg.appendChild(createSVGElement("circle", {
                "cx": xSVG, "cy": ySVG, "r": circleRadius,
                "fill": "none", "stroke": "black", "stroke-width": curveWidth
            }));

            latexCode = latexCode + `\n\\draw[\\colorNonTransversal] (${cleanNumber(x)},${cleanNumber(y)}) circle (${circleRadiusTikz});`;
        }

        // Put the label into the largest angle between the edges at the point.
        // The SVG is not distorted, so the angles agree with TikZ (up to the y-axis pointing down).
        const angles = [];

        for (const edge of edges) {
            if (edge.start !== index && edge.end !== index) {
                continue;
            }

            const other = points[(edge.start === index) ? edge.end : edge.start];
            angles.push(Math.atan2(other[1] - y, other[0] - x));
        }

        const angle = freeDirection(angles);
        // Labels of circled points are placed outside of the circle
        const dist = circled ? circleRadius + multLabelDist : multLabelDist;

        const label = createSVGElement("text", {
            "x": xSVG + dist * Math.cos(angle),
            "y": ySVG - dist * Math.sin(angle),
            "text-anchor": "middle", "dominant-baseline": "central",
            "font-size": labelSize, "fill": "black"
        });
        label.textContent = `${info.multiplicity}`;
        svg.appendChild(label);

        // The node is placed in direction angle by anchoring it at the opposite side
        const anchor = Math.round(angle * 180 / Math.PI + 180) % 360;
        const outerSep = circled ? `, outer sep=${circleRadiusTikz}` : "";
        latexCode = latexCode + `\n\\node[\\colorMult, anchor=${anchor}${outerSep}] at (${cleanNumber(x)},${cleanNumber(y)}) {${info.multiplicity}};`;
    });
}

function drawDualSubdivision(svgDS, P, edges, vertexInfo) {
    /**
     * Draws the dual subdivision of the union of the curves: Edges are colored
     * like their dual edges, cells dual to vertices of only one curve are filled
     * in its color, the parallelograms dual to intersection points stay white.
     * Appends the corresponding TikZ code to latexCodeDS.
     *
     * @param {SVGElement} svgDS - The plot of the dual subdivision.
     * @param {Array<[number, number, number]>} P - The monomials of F*G.
     * @param {Object[]} edges - The edges, see computeEdges.
     * @param {Object[]} vertexInfo - The classification of the vertices, see classifyVertices.
     */
    let maxX = 0;
    let maxY = 0;
    let maxDeg = 0;

    for (const [i, j] of P) {
        maxX = Math.max(maxX, i);
        maxY = Math.max(maxY, j);
        maxDeg = Math.max(maxDeg, i + j);
    }

    const sizeDS = Math.max(maxX, maxY);
    // If the Newton polytope is wider than high, move it to the top of the plot
    const shiftY = sizeDS - maxY;
    const toSVG = ([i, j]) => [xSVGCoords(i, -0.15, sizeDS + 0.15), ySVGCoords(j + shiftY, -0.15, sizeDS + 0.15)];

    if (edges.some(edge => edge.curve === "both")) {
        latexCodeDS = latexCodeDS + "\n\\def\\colorDSOverlap{violet}";
    }

    // Fill the cells dual to vertices of only one curve (first, such that edges and points lie on top)
    for (const info of vertexInfo) {
        if (info.type !== "one" && info.type !== "two") {
            continue;
        }

        const style = curveStyles[info.type];

        svgDS.appendChild(createSVGElement("polygon", {
            "points": info.cell.map(pt => toSVG(pt).join(",")).join(" "),
            "fill": style.fill
        }));

        latexCodeDS = latexCodeDS + `\n\\fill[fill=${style.tikzDS}!25] ${info.cell.map(([i, j]) => `(${i},${j})`).join(" -- ")} -- cycle;`;
    }

    // Dual edges
    for (const edge of edges) {
        const style = curveStyles[edge.curve];
        const [[i1, j1], [i2, j2]] = edge.dual;
        const [x1SVG, y1SVG] = toSVG([i1, j1]);
        const [x2SVG, y2SVG] = toSVG([i2, j2]);

        svgDS.appendChild(createSVGElement("line", {
            "x1": x1SVG, "y1": y1SVG, "x2": x2SVG, "y2": y2SVG,
            "stroke": style.color, "stroke-width": 0.02
        }));

        latexCodeDS = latexCodeDS + `\n\\draw[${style.tikzDS}] (${i1},${j1}) -- (${i2},${j2});`;
    }

    // Lattice points: As in the TikZ figures of the lecture notes, all (i,j) with i <= maxX, j <= maxY and i + j <= maxDeg
    for (let i = 0; i <= maxX; i++) {
        for (let j = 0; j <= maxY && i + j <= maxDeg; j++) {
            const [xSVG, ySVG] = toSVG([i, j]);

            svgDS.appendChild(createSVGElement("circle", {
                "cx": xSVG, "cy": ySVG, "r": 0.15, "fill": "black"
            }));

            latexCodeDS = latexCodeDS + `\n\\filldraw[color = \\colorDSDots, fill=\\colorDSDots] (${i},${j}) circle (1.5pt);`;
        }
    }
}

/*
Main function
*/

function plotTropIntersection(polyStrF, polyStrG, aX, bX, aY, bY, useMin, autoAdjust, autoBoundaryDist, autoSquare, showGrid, showGridLabels) {
    const svg = document.getElementById("plot");
    const svgDS = document.getElementById("plotDS");

    // Delete error messages and old plots
    showError("");
    svg.replaceChildren();
    svgDS.replaceChildren();

    // Get polynomials
    let F;
    let G;

    try {
        F = parsePolynomialTrop(polyStrF, useMin);
        G = parsePolynomialTrop(polyStrG, useMin);
    } catch (error) {
        showError(error.message);
        return;
    }

    // Catch case of only one monomial (then the curve is empty)
    if (F.length < 2 || G.length < 2) {
        showError("Please enter at least two different monomials for each curve.");
        return;
    }

    // With automatic range, vertices must not lie on the boundary of the plot
    if (autoAdjust && !(autoBoundaryDist > 0)) {
        showError("The minimal distance of vertices to the boundary must be positive.");
        return;
    }

    // Compute everything in the max-convention: For min, change all signs and mirror the
    // plotting range. Before plotting, all points are mirrored back.
    const sign = useMin ? -1 : 1;

    if (useMin) {
        F = changeSignsOfPolynomial(F);
        G = changeSignsOfPolynomial(G);
        [aX, bX, aY, bY] = [-bX, -aX, -bY, -aY];
    }

    const P = multiplyTropPolynomials(F, G);

    //// Vertices of the union of the curves, including the intersection points
    const vertices = computeVertices(P);

    //// Plotting range
    if (vertices.length > 0) {
        if (autoAdjust) {
            aX = Math.min(...vertices.map(pt => pt[0]));
            bX = Math.max(...vertices.map(pt => pt[0]));
            aY = Math.min(...vertices.map(pt => pt[1]));
            bY = Math.max(...vertices.map(pt => pt[1]));
        } else {
            if (!(aX < bX) || !(aY < bY)) {
                showError("Invalid boundaries.");
                return;
            }

            // All vertices have to be in the plotting range
            if (vertices.some(([x, y]) => x <= aX || x >= bX || y <= aY || y >= bY)) {
                // Suggest a range, in the coordinates of the chosen convention
                const xs = vertices.map(pt => sign * pt[0]);
                const ys = vertices.map(pt => sign * pt[1]);

                showError(`All vertices and intersection points must be in the plotting range. Suggested range: [${Math.floor(Math.min(...xs)) - 1},${Math.ceil(Math.max(...xs)) + 1}] x [${Math.floor(Math.min(...ys)) - 1},${Math.ceil(Math.max(...ys)) + 1}]`);
                return;
            }
        }
    } else if (autoAdjust) {
        // We only have parallel lines: Use points on them
        const linePoints = computeLinePoints(P);

        aX = Math.min(...linePoints.map(pt => pt[0]));
        bX = Math.max(...linePoints.map(pt => pt[0]));
        aY = Math.min(...linePoints.map(pt => pt[1]));
        bY = Math.max(...linePoints.map(pt => pt[1]));
    }

    if (autoAdjust && autoSquare) {
        const size = Math.max(bX - aX, bY - aY);
        const centerX = (aX + bX) / 2;
        const centerY = (aY + bY) / 2;

        aX = centerX - size / 2 - autoBoundaryDist;
        bX = centerX + size / 2 + autoBoundaryDist;
        aY = centerY - size / 2 - autoBoundaryDist;
        bY = centerY + size / 2 + autoBoundaryDist;
    } else if (autoAdjust) {
        aX = aX - autoBoundaryDist;
        bX = bX + autoBoundaryDist;
        aY = aY - autoBoundaryDist;
        bY = bY + autoBoundaryDist;
    }

    if (!(aX < bX) || !(aY < bY)) {
        showError("Invalid boundaries.");
        return;
    }

    //// All points at which edges start and end: The vertices (first) and the points on the boundary of the plot
    const points = [...vertices];

    for (const pt of computeBoundaryPoints(P, aX, bX, aY, bY)) {
        if (tupleIsInList(pt, points, eps) === -1) {
            points.push(pt);
        }
    }

    // The maximal monomials of F*G at these points
    const maxSets = points.map(([x, y]) => maxIndices(x, y, P));

    const edges = computeEdges(points, maxSets, P, F, G);
    const vertexInfo = classifyVertices(vertices, maxSets, P, F, G);

    //// Back to the coordinates of the chosen convention
    const plotPoints = points.map(([x, y]) => [sign * x, sign * y]);

    if (useMin) {
        [aX, bX, aY, bY] = [-bX, -aX, -bY, -aY];
    }

    reportNonTransversality(plotPoints, edges, vertexInfo);

    //// Plotting
    // Coordinate lines first, such that the curves lie on top of them
    if (showGrid) {
        drawCoordinateLines(svg, aX, bX, aY, bY, showGridLabels);
    }

    drawCurves(svg, plotPoints, edges, vertexInfo, aX, bX, aY, bY);
    drawDualSubdivision(svgDS, P, edges, vertexInfo);
}

/*
Plot the computed curves
*/

function plotte(){
    const polyStrF = document.getElementById("polyF").value;
    const polyStrG = document.getElementById("polyG").value;
    const aX = document.getElementById("aX").valueAsNumber;
    const bX = document.getElementById("bX").valueAsNumber;
    const aY = document.getElementById("aY").valueAsNumber;
    const bY = document.getElementById("bY").valueAsNumber;
    const useMin = !document.getElementById("useMax").checked;
    const autoAdjust = document.getElementById("autoAdjust").checked;
    const autoBoundaryDist = document.getElementById("autoBoundaryDist").valueAsNumber;
    const autoSquare = document.getElementById("autoSquare").checked;
    const showGrid = document.getElementById("showGrid").checked;
    const showGridLabels = document.getElementById("showGridLabels").checked;

    // Reset the initial LaTeX code: Color declarations for easier changes
    latexCode = "\\def\\colorCurveOne{blue}\n\\def\\colorCurveTwo{red}\n\\def\\colorMult{black}";
    latexCodeDS = "\\def\\colorDSOne{blue}\n\\def\\colorDSTwo{red}\n\\def\\colorDSDots{black}";

    plotTropIntersection(polyStrF, polyStrG, aX, bX, aY, bY, useMin, autoAdjust, autoBoundaryDist, autoSquare, showGrid, showGridLabels);
}

/*
Functions that handle the websites behaviour
*/

function shiftCurve(polyId, shiftXId, shiftYId) {
    /*
     * Shifts a curve by the vector (s, t) from the shift fields by changing its polynomial
     * in the input field: F(x, y) is replaced by F(x - s, y - t), i.e. the coefficient of
     * x^a y^b changes by -a*s - b*t. This works for the max- and the min-convention.
     * The order of the terms in the input is kept.
     */
    const polyInput = document.getElementById(polyId);
    // Empty shift fields count as 0
    const s = document.getElementById(shiftXId).valueAsNumber || 0;
    const t = document.getElementById(shiftYId).valueAsNumber || 0;

    let terms;

    try {
        terms = parseTermsTrop(polyInput.value);
    } catch (error) {
        showError(error.message);
        return;
    }

    if (terms.length === 0) {
        return;
    }

    polyInput.value = formatPolynomialTrop(terms.map(([a, b, c]) => [a, b, c - a * s - b * t]));

    // Reset the shift fields
    document.getElementById(shiftXId).value = 0;
    document.getElementById(shiftYId).value = 0;

    // Show the shifted curve
    plotIfComplete();
}

function plotIfComplete() {
    // Plots, if both polynomials are entered
    if (document.getElementById("polyF").value.trim() !== "" && document.getElementById("polyG").value.trim() !== "") {
        plotte();
    }
}

function addShiftListeners(polyId, shiftXId, shiftYId, buttonId) {
    // The shift button and Enter in the shift fields shift the curve (instead of submitting the form)
    const shift = () => shiftCurve(polyId, shiftXId, shiftYId);

    document.getElementById(buttonId).addEventListener("click", shift);

    for (const id of [shiftXId, shiftYId]) {
        document.getElementById(id).addEventListener("keydown", function(event) {
            if (event.key === "Enter") {
                event.preventDefault();
                shift();
            }
        });
    }
}

function showLatex() {
    document.getElementById("latexOutput").value = latexCode;
    document.getElementById("latexPopup").style.display = "block";
}

function closeLatex() {
    document.getElementById("latexPopup").style.display = "none";
}

function copyLatex() {
    const latex = document.getElementById("latexOutput").value;
    navigator.clipboard.writeText(latex);
}

function showLatexDS() {
    document.getElementById("latexOutputDS").value = latexCodeDS;
    document.getElementById("latexPopupDS").style.display = "block";
}

function closeLatexDS() {
    document.getElementById("latexPopupDS").style.display = "none";
}

function copyLatexDS() {
    const latex = document.getElementById("latexOutputDS").value;
    navigator.clipboard.writeText(latex);
}
