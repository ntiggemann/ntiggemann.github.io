/*
Plotting of tropical surfaces: website behaviour, 3D view, TikZ code.
The computations are in tropSurfaceComputation.js.
*/

/*
Eventlisteners
*/

const plotForm = document.getElementById("plotForm");

plotForm.addEventListener("submit", function(event) {
    event.preventDefault();
    plotte();
});

document.getElementById("copyLatexButton").addEventListener("click", copyLatex);
document.getElementById("copyLatexButtonDS").addEventListener("click", copyLatexDS);
document.getElementById("latexButton").addEventListener("click", showLatex);
document.getElementById("latexButtonDS").addEventListener("click", showLatexDS);
document.getElementById("closeLatexButton").addEventListener("click", closeLatex);
document.getElementById("closeLatexButtonDS").addEventListener("click", closeLatexDS);

// Re-plot whenever one of the checkbox options is toggled (only if a polynomial has been entered)
plotForm.addEventListener("change", function(event) {
    if (event.target.type === "checkbox" && document.getElementById("poly").value.trim() !== "") {
        plotte();
    }
});
for (const id of ["viewTheta", "viewPhi"]) {
    document.getElementById(id).addEventListener("input", function() {
        readView();
        render();
    });
}
document.getElementById("resetViewButton").addEventListener("click", function() {
    setView(defaultView.theta, defaultView.phi);
    render();
});

/*
State
*/

const eps = 1e-9; // Relative tolerance: how close two numbers must be to be considered equal

const defaultView = { theta: 60, phi: 30 }; // About the default view of matplotlib
const view = { theta: defaultView.theta, phi: defaultView.phi };

let result = null; // The last computed surface

/*
Projection
*/

function viewVectors(theta, phi) {
    /**
     * The same projection as \tdplotsetmaincoords{theta}{phi} in tikz-3dplot.
     *
     * @returns {{ex: number[], ey: number[], ez: number[]}} A point p is drawn at (ex.p, ey.p) (y pointing up);
     *          ez points towards the viewer.
     */
    const t = theta * Math.PI / 180;
    const f = phi * Math.PI / 180;
    return {
        ex: [Math.cos(f), Math.sin(f), 0],
        ey: [-Math.cos(t) * Math.sin(f), Math.cos(t) * Math.cos(f), Math.sin(t)],
        ez: [Math.sin(t) * Math.sin(f), -Math.sin(t) * Math.cos(f), Math.cos(t)]
    };
}

function makeProjection(lower, upper, size) {
    /**
     * Returns a function mapping a point in R^3 to SVG coordinates {x, y, depth}. The box [lower, upper] fits into a
     * disc of radius size around the center of the SVG, for every view. Larger depth means closer to the viewer.
     */
    const { ex, ey, ez } = viewVectors(view.theta, view.phi);
    const center = lower.map((a, m) => (a + upper[m]) / 2);
    const radius = Math.max(Math.hypot(...upper.map((b, m) => b - lower[m])) / 2, 1e-12);
    const s = size / radius;
    return function(p) {
        const q = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
        return {
            x: s * (ex[0] * q[0] + ex[1] * q[1] + ex[2] * q[2]),
            y: -s * (ey[0] * q[0] + ey[1] * q[1] + ey[2] * q[2]),
            depth: ez[0] * q[0] + ez[1] * q[1] + ez[2] * q[2]
        };
    };
}

/*
Helpers for drawing
*/

function createSVGElement(tag, attributes) {
    // Creates an SVG element with the given attributes
    const element = document.createElementNS("http://www.w3.org/2000/svg", tag);

    for (const [name, value] of Object.entries(attributes)) {
        element.setAttribute(name, value);
    }

    return element;
}

function svgLine(svg, a, b, color, width, extra = {}) {
    svg.appendChild(createSVGElement("line", {
        "x1": a.x, "y1": a.y, "x2": b.x, "y2": b.y,
        "stroke": color, "stroke-width": width, "stroke-linecap": "round", ...extra
    }));
}

function svgText(svg, a, text, size, color) {
    const label = createSVGElement("text", {
        "x": a.x, "y": a.y,
        "text-anchor": "middle", "dominant-baseline": "central",
        "font-size": size, "fill": color
    });
    label.textContent = text;
    svg.appendChild(label);
}

function fmt(v) {
    // Short number format for TikZ code
    return `${Number(v.toFixed(4)) + 0}`;
}

function tikzPoint(p) {
    return `(${fmt(p[0])},${fmt(p[1])},${fmt(p[2])})`;
}

function tikzAnchor(dx, dy) {
    // The TikZ anchor of a label that is shifted in direction (dx, dy) (y pointing up) from its point
    const angle = Math.atan2(dy, dx) * 180 / Math.PI;
    const anchors = ["west", "south west", "south", "south east", "east", "north east", "north", "north west"];
    return anchors[((Math.round(angle / 45) % 8) + 8) % 8];
}

/*
Coordinate lines (as in the curve plotter)
*/

// The number of intervals between coordinate lines the plotting range is divided into (approximately)
const gridTargetIntervals = 5;

function niceGridStep(a, b) {
    /**
     * Computes a "nice" distance between coordinate lines for the interval [a, b]:
     * 1, 2 or 5 times a power of ten, such that [a, b] is divided into roughly gridTargetIntervals pieces.
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
    // All integer multiples of step in the interval [a, b], rounded to remove floating point noise
    const tolerance = 1e-9;
    const decimals = Math.max(0, -Math.floor(Math.log10(step) + tolerance));
    const values = [];

    for (let k = Math.ceil(a / step - tolerance); k <= Math.floor(b / step + tolerance); k++) {
        values.push(Number((k * step).toFixed(decimals)) + 0);
    }

    return values;
}

/*
The plot of the surface
*/

function boxDecorations(lower, upper, project, showGrid, showGridLabels) {
    /**
     * The back sides of the box, coordinate lines, labels and names of the axes (for the current view).
     *
     * @returns {{lines: Array, labels: Array}} lines: {a, b, kind} with kind "box" or "grid" (points in R^3);
     *          labels: {p, text, kind, dx, dy} with kind "grid" or "axis", (dx, dy) the shift in SVG coordinates.
     */
    const { ez } = viewVectors(view.theta, view.phi);
    // The coordinate of the back side in direction m (the side facing away from the viewer) and of the front side
    const back = [0, 1, 2].map(m => ez[m] > 0 ? lower[m] : upper[m]);
    const front = [0, 1, 2].map(m => ez[m] > 0 ? upper[m] : lower[m]);
    const lines = [];
    const labels = [];

    // Edges of the box that lie on a back side
    for (let m = 0; m < 3; m++) {
        const [r, s] = [0, 1, 2].filter(k => k !== m);
        for (const vr of [lower[r], upper[r]]) {
            for (const vs of [lower[s], upper[s]]) {
                if (vr === back[r] || vs === back[s]) {
                    const a = [0, 0, 0];
                    a[r] = vr;
                    a[s] = vs;
                    const b = a.slice();
                    a[m] = lower[m];
                    b[m] = upper[m];
                    lines.push({ a, b, kind: "box" });
                }
            }
        }
    }

    // Coordinate lines on the back sides
    const steps = [0, 1, 2].map(m => niceGridStep(lower[m], upper[m]));
    if (showGrid && steps.every(Number.isFinite)) {
        for (let m = 0; m < 3; m++) { // the back side x_m = back[m]
            for (const r of [0, 1, 2].filter(k => k !== m)) { // lines x_r = const
                const s = 3 - m - r;
                for (const value of gridValues(lower[r], upper[r], steps[r])) {
                    const a = [0, 0, 0];
                    a[m] = back[m];
                    a[r] = value;
                    const b = a.slice();
                    a[s] = lower[s];
                    b[s] = upper[s];
                    lines.push({ a, b, kind: "grid" });
                }
            }
        }
    }

    // Where the axes are labelled: x and y along the bottom (back side in z) at the front, z at the rightmost vertical edge
    const labelEdges = [];
    labelEdges.push({ m: 0, base: [0, front[1], back[2]] });
    labelEdges.push({ m: 1, base: [front[0], 0, back[2]] });
    let best = null;
    for (const vx of [lower[0], upper[0]]) {
        for (const vy of [lower[1], upper[1]]) {
            const x = project([vx, vy, (lower[2] + upper[2]) / 2]).x;
            if (best === null || x > best.x) {
                best = { x, base: [vx, vy, 0] };
            }
        }
    }
    labelEdges.push({ m: 2, base: best.base });

    const center = project(lower.map((a, m) => (a + upper[m]) / 2));
    const names = ["x", "y", "z"];

    for (const { m, base } of labelEdges) {
        // Shift the labels in SVG coordinates orthogonally to the projected edge, away from the center of the box
        const start = base.slice();
        const end = base.slice();
        start[m] = lower[m];
        end[m] = upper[m];
        const a = project(start);
        const b = project(end);
        let normal = [a.y - b.y, b.x - a.x];
        const length = Math.hypot(normal[0], normal[1]);
        normal = length > 1e-9 ? [normal[0] / length, normal[1] / length] : [0, 1];
        if (normal[0] * ((a.x + b.x) / 2 - center.x) + normal[1] * ((a.y + b.y) / 2 - center.y) < 0) {
            normal = [-normal[0], -normal[1]];
        }

        const withGridLabels = showGrid && showGridLabels && Number.isFinite(steps[m]);
        if (withGridLabels) {
            for (const value of gridValues(lower[m], upper[m], steps[m])) {
                const p = base.slice();
                p[m] = value;
                labels.push({ p, text: `${value}`, kind: "grid", dx: 0.5 * normal[0], dy: 0.5 * normal[1] });
            }
        }
        const p = base.slice();
        p[m] = (lower[m] + upper[m]) / 2;
        const shift = withGridLabels ? 1.4 : 0.7;
        labels.push({ p, text: names[m], kind: "axis", dx: shift * normal[0], dy: shift * normal[1], behindGridLabels: withGridLabels });
    }

    return { lines, labels };
}

function renderSurface() {
    const svg = document.getElementById("plot");
    svg.replaceChildren(); // Clears the plot
    if (result === null) {
        return;
    }

    const { lower, upper } = result;
    const project = makeProjection(lower, upper, 8.3);
    const showBox = document.getElementById("showBox").checked;
    const showGrid = document.getElementById("showGrid").checked;
    const showGridLabels = document.getElementById("showGridLabels").checked;

    // Box and coordinate lines first, such that the surface lies on top of them
    if (showBox || showGrid) {
        const { lines, labels } = boxDecorations(lower, upper, project, showGrid, showGridLabels);
        for (const { a, b, kind } of lines) {
            if (kind === "box" && !showBox) {
                continue;
            }
            svgLine(svg, project(a), project(b), kind === "box" ? "#808080" : "#d0d0d0", kind === "box" ? 0.03 : 0.02);
        }
        for (const { p, text, kind, dx, dy } of labels) {
            if (kind === "axis" && !showBox) {
                continue;
            }
            const q = project(p);
            svgText(svg, { x: q.x + dx, y: q.y + dy }, text, kind === "axis" ? 0.6 : 0.4, kind === "axis" ? "black" : "#808080");
        }
    }

    // The polygons, from back to front
    const polygons = result.cells.map(cell => {
        const pts = cell.vertices.map(project);
        return { pts, depth: pts.reduce((sum, q) => sum + q.depth, 0) / pts.length };
    });
    polygons.sort((a, b) => a.depth - b.depth);
    for (const { pts } of polygons) {
        svg.appendChild(createSVGElement("polygon", {
            "points": pts.map(q => `${q.x},${q.y}`).join(" "),
            "fill": "blue", "fill-opacity": 0.5, "stroke": "none"
        }));
    }

    // Edges and isolated points
    for (const { p, q, onCell } of result.edges) {
        svgLine(svg, project(p), project(q), onCell ? "blue" : "black", 0.04);
    }
    for (const { p, onCell } of result.lonePoints) {
        const q = project(p);
        svg.appendChild(createSVGElement("circle", { "cx": q.x, "cy": q.y, "r": 0.12, "fill": onCell ? "blue" : "black" }));
    }

    // Weights of the 2-cells
    for (const cell of result.cells) {
        if (cell.weight > 1) {
            const q = project(cell.vertices.reduce((s, p) => p.map((t, m) => s[m] + t / cell.vertices.length), [0, 0, 0]));
            svgText(svg, q, `${cell.weight}`, 0.5, "black");
        }
    }
}

/*
The plot of the dual subdivision
*/

function dualData() {
    // Lattice points, monomials of the polynomial and the box [0, max exponents]
    const maxExp = [0, 1, 2].map(m => Math.max(...result.monomials.map(mon => mon[m])));
    const support = new Set(result.monomials.map(mon => `${mon[0]},${mon[1]},${mon[2]}`));
    const lattice = [];
    for (let i = 0; i <= maxExp[0]; i++) {
        for (let j = 0; j <= maxExp[1]; j++) {
            for (let k = 0; k <= maxExp[2]; k++) {
                lattice.push({ p: [i, j, k], inSupport: support.has(`${i},${j},${k}`) });
            }
        }
    }
    // The box that has to fit into the plot: lattice points and the axes with their names
    return { maxExp, lattice, lower: [-0.3, -0.3, -0.3], upper: maxExp.map(v => v + 1) };
}

function renderDual() {
    const svg = document.getElementById("plotDS");
    svg.replaceChildren(); // Clears the plot
    if (result === null) {
        return;
    }

    const { maxExp, lattice, lower, upper } = dualData();
    const project = makeProjection(lower, upper, 8.3);
    const names = ["x", "y", "z"];

    // Axes through the origin
    const origin = project([0, 0, 0]);
    for (let m = 0; m < 3; m++) {
        const end = [0, 0, 0];
        end[m] = maxExp[m] + 0.6;
        svgLine(svg, origin, project(end), "#b0b0b0", 0.03);
        end[m] = maxExp[m] + 0.9;
        svgText(svg, project(end), names[m], 0.6, "black");
    }

    // Lattice points not in the support, then the edges, then the monomials of the polynomial
    for (const { p, inSupport } of lattice) {
        if (!inSupport) {
            const q = project(p);
            svg.appendChild(createSVGElement("circle", { "cx": q.x, "cy": q.y, "r": 0.08, "fill": "#c0c0c0" }));
        }
    }
    for (const { from, to } of result.dualEdges) {
        svgLine(svg, project(from), project(to), "black", 0.04);
    }
    for (const { p, inSupport } of lattice) {
        if (inSupport) {
            const q = project(p);
            svg.appendChild(createSVGElement("circle", { "cx": q.x, "cy": q.y, "r": 0.15, "fill": "black" }));
        }
    }
}

function render() {
    renderSurface();
    renderDual();
}

/*
Rotating the plots with the mouse (or a finger)
*/

let drag = null;
let renderPending = false;

function setView(theta, phi) {
    view.theta = Math.min(Math.max(theta, 0), 180);
    view.phi = ((phi % 360) + 360) % 360;
    document.getElementById("viewTheta").value = Math.round(view.theta);
    document.getElementById("viewPhi").value = Math.round(view.phi);
}

function readView() {
    const theta = document.getElementById("viewTheta").valueAsNumber;
    const phi = document.getElementById("viewPhi").valueAsNumber;
    if (Number.isFinite(theta)) {
        view.theta = Math.min(Math.max(theta, 0), 180);
    }
    if (Number.isFinite(phi)) {
        view.phi = phi;
    }
}

for (const id of ["plot", "plotDS"]) {
    const svg = document.getElementById(id);
    svg.addEventListener("pointerdown", function(event) {
        drag = { x: event.clientX, y: event.clientY };
        svg.setPointerCapture(event.pointerId);
    });
    svg.addEventListener("pointermove", function(event) {
        if (drag === null) {
            return;
        }
        // Degrees per pixel, such that dragging over the whole plot turns by about 180 degrees
        const degreesPerPixel = 180 / Math.max(svg.getBoundingClientRect().width, 100);
        setView(view.theta - (event.clientY - drag.y) * degreesPerPixel, view.phi - (event.clientX - drag.x) * degreesPerPixel);
        drag = { x: event.clientX, y: event.clientY };
        if (!renderPending) {
            renderPending = true;
            requestAnimationFrame(function() {
                renderPending = false;
                render();
            });
        }
    });
    for (const type of ["pointerup", "pointercancel"]) {
        svg.addEventListener(type, function() {
            drag = null;
        });
    }
}

/*
TikZ code
*/

function tikzHeader(lower, upper) {
    const scale = Number((6 / Math.max(...upper.map((b, m) => b - lower[m]))).toPrecision(2));
    return [
        "% Needs \\usepackage{tikz-3dplot}. Put this into a tikzpicture. The view is the one in the browser.",
        `\\tdplotsetmaincoords{${fmt(view.theta)}}{${fmt(view.phi)}}`,
        `\\begin{scope}[tdplot_main_coords, scale=${scale}]`
    ];
}

function latexSurface() {
    if (result === null) {
        return "";
    }
    const { lower, upper } = result;
    const project = makeProjection(lower, upper, 8.3);
    const showBox = document.getElementById("showBox").checked;
    const showGrid = document.getElementById("showGrid").checked;
    const showGridLabels = document.getElementById("showGridLabels").checked;

    const code = [
        "\\def\\colorSurface{blue}",
        "\\def\\colorCurve{black}",
        "\\def\\colorWeights{black}"
    ];
    if (showBox) {
        code.push("\\def\\colorBox{black!50}");
    }
    if (showGrid) {
        code.push("\\def\\colorGrid{black!20}");
        if (showGridLabels) {
            code.push("\\def\\colorGridLabels{black!60}");
        }
    }
    code.push(...tikzHeader(lower, upper));

    if (showBox || showGrid) {
        const { lines, labels } = boxDecorations(lower, upper, project, showGrid, showGridLabels);
        for (const { a, b, kind } of lines) {
            if (kind === "box" && showBox) {
                code.push(`\\draw[\\colorBox, thin] ${tikzPoint(a)} -- ${tikzPoint(b)};`);
            } else if (kind === "grid") {
                code.push(`\\draw[\\colorGrid, very thin] ${tikzPoint(a)} -- ${tikzPoint(b)};`);
            }
        }
        for (const { p, text, kind, dx, dy, behindGridLabels } of labels) {
            if (kind === "grid") {
                code.push(`\\node[\\colorGridLabels, anchor=${tikzAnchor(dx, -dy)}, font=\\scriptsize] at ${tikzPoint(p)} {$${text}$};`);
            } else if (showBox) {
                // Names of the axes: further outside, if there are labels of coordinate lines
                let shift = "";
                if (behindGridLabels) {
                    const d = Math.hypot(dx, dy);
                    shift = `, xshift=${fmt(1.5 * dx / d)}em, yshift=${fmt(-1.5 * dy / d)}em`;
                }
                code.push(`\\node[anchor=${tikzAnchor(dx, -dy)}${shift}] at ${tikzPoint(p)} {$${text}$};`);
            }
        }
    }

    // Polygons from back to front, as in the browser
    const polygons = result.cells.map(cell => ({
        cell,
        depth: cell.vertices.reduce((sum, p) => sum + project(p).depth, 0) / cell.vertices.length
    }));
    polygons.sort((a, b) => a.depth - b.depth);
    for (const { cell } of polygons) {
        code.push(`\\fill[\\colorSurface, fill opacity=0.5] ${cell.vertices.map(tikzPoint).join(" -- ")} -- cycle;`);
    }
    for (const { p, q, onCell } of result.edges) {
        code.push(`\\draw[${onCell ? "\\colorSurface" : "\\colorCurve"}] ${tikzPoint(p)} -- ${tikzPoint(q)};`);
    }
    for (const { p, onCell } of result.lonePoints) {
        code.push(`\\fill[${onCell ? "\\colorSurface" : "\\colorCurve"}] ${tikzPoint(p)} circle (1.5pt);`);
    }
    for (const cell of result.cells) {
        if (cell.weight > 1) {
            const c = cell.vertices.reduce((s, p) => p.map((t, m) => s[m] + t / cell.vertices.length), [0, 0, 0]);
            code.push(`\\node[\\colorWeights] at ${tikzPoint(c)} {${cell.weight}};`);
        }
    }
    code.push("\\end{scope}");
    return code.join("\n");
}

function latexDual() {
    if (result === null) {
        return "";
    }
    const { maxExp, lattice, lower, upper } = dualData();
    const code = [
        "\\def\\colorDSLines{black}",
        "\\def\\colorDSDots{black}",
        "\\def\\colorDSLattice{black!25}",
        "\\def\\colorDSAxes{black!30}"
    ];
    code.push(...tikzHeader(lower, upper));
    const names = ["x", "y", "z"];
    for (let m = 0; m < 3; m++) {
        const end = [0, 0, 0];
        end[m] = maxExp[m] + 0.6;
        code.push(`\\draw[\\colorDSAxes] (0,0,0) -- ${tikzPoint(end)} node[black, pos=1.25] {$${names[m]}$};`);
    }
    for (const { p, inSupport } of lattice) {
        if (!inSupport) {
            code.push(`\\fill[\\colorDSLattice] ${tikzPoint(p)} circle (1pt);`);
        }
    }
    for (const { from, to } of result.dualEdges) {
        code.push(`\\draw[\\colorDSLines] ${tikzPoint(from)} -- ${tikzPoint(to)};`);
    }
    for (const { p, inSupport } of lattice) {
        if (inSupport) {
            code.push(`\\filldraw[color=\\colorDSDots, fill=\\colorDSDots] ${tikzPoint(p)} circle (1.5pt);`);
        }
    }
    code.push("\\end{scope}");
    return code.join("\n");
}

/*
Plot the computed surface
*/

function showError(errorText) {
    document.getElementById("error").textContent = errorText;
}

function showInfo(infoText) {
    document.getElementById("info").textContent = infoText;
}

function readNumber(id) {
    return document.getElementById(id).valueAsNumber;
}

function plotte() {
    const polyStr = document.getElementById("poly").value;
    const useMin = !document.getElementById("useMax").checked;
    const autoAdjust = document.getElementById("autoAdjust").checked;
    const autoCube = document.getElementById("autoCube").checked;
    const autoBoundaryDist = readNumber("autoBoundaryDist");
    const lower = [readNumber("aX"), readNumber("aY"), readNumber("aZ")];
    const upper = [readNumber("bX"), readNumber("bY"), readNumber("bZ")];

    // Delete messages
    showError("");
    showInfo("");

    let monomials;
    try {
        monomials = parsePolynomialTrop(polyStr, useMin);
    } catch (error) {
        showError(error.message);
        return;
    }

    // Catch case of only one monomial
    if (monomials.length < 2) {
        showError("Please enter at least two different monomials.");
        return;
    }

    if (autoAdjust && !(autoBoundaryDist >= 0)) {
        showError("The minimal distance of vertices to the boundary must be a number >= 0.");
        return;
    }
    if (!autoAdjust && ![0, 1, 2].every(m => Number.isFinite(lower[m]) && Number.isFinite(upper[m]) && lower[m] < upper[m])) {
        showError("Invalid boundaries.");
        return;
    }

    readView();
    result = computeTropicalSurface(monomials, { useMin, autoAdjust, autoCube, autoBoundaryDist, lower, upper, eps });

    // Information about the result
    const { cells, edges, lonePoints, vertices, newtonDim, verticesOutside } = result;
    const range = [0, 1, 2].map(m => `[${result.lower[m]}, ${result.upper[m]}]`).join(" x ");
    let info = `Plotting range ${range}. `;
    if (newtonDim === 3) {
        if (vertices.length === 1) {
            info += verticesOutside > 0 ? "The surface has 1 vertex, outside the plotting range." : "The surface has 1 vertex.";
        } else {
            info += `The surface has ${vertices.length} vertices`;
            info += verticesOutside > 0 ? `, ${verticesOutside} of them outside the plotting range.` : ".";
        }
    } else {
        info += `The Newton polytope is ${newtonDim}-dimensional, so the surface has no vertices.`;
    }
    if (cells.length === 0 && edges.length === 0 && lonePoints.length === 0) {
        showError("The surface does not intersect the plotting range.");
    }
    showInfo(info);

    render();
}

/*
Functions that handle the websites behaviour
*/

function showLatex() {
    document.getElementById("latexOutput").value = latexSurface();
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
    document.getElementById("latexOutputDS").value = latexDual();
    document.getElementById("latexPopupDS").style.display = "block";
}

function closeLatexDS() {
    document.getElementById("latexPopupDS").style.display = "none";
}

function copyLatexDS() {
    const latex = document.getElementById("latexOutputDS").value;
    navigator.clipboard.writeText(latex);
}
