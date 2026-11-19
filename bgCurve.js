// bgCurve.js
// Draws a random smooth plane tropical curve (max-convention) into the page background.
// Loaded with "defer" on the content pages, not on the plotters.
(function () {
    "use strict";

    const EPS = 1e-9;
    const MAX_TRIES = 300;
    const STROKE = "#c4c4c4";
    const STROKE_WIDTH = 2;

    // Random tropical polynomial of degree d: list of monomials [i, j, c] with i + j <= d.
    // A strictly concave lifting plus noise gives a regular, mostly unimodular subdivision.
    function randomPolynomial(d) {
        const A = [];
        for (let i = 0; i <= d; i++) {
            for (let j = 0; i + j <= d; j++) {
                const c = -(i * i + i * j + j * j) + 1.6 * (Math.random() - 0.5);
                A.push([i, j, c]);
            }
        }
        return A;
    }

    function evaluate(A, x, y) {
        let m = -Infinity;
        for (const [i, j, c] of A) m = Math.max(m, c + i * x + j * y);
        return m;
    }

    // Vertices of the curve: points where exactly three monomials are maximal.
    // Returns null unless the dual subdivision is a unimodular triangulation (smooth curve).
    function computeCurve(A, d) {
        const vertices = [];
        for (let a = 0; a < A.length; a++) {
            for (let b = a + 1; b < A.length; b++) {
                for (let c = b + 1; c < A.length; c++) {
                    const [i1, j1, c1] = A[a], [i2, j2, c2] = A[b], [i3, j3, c3] = A[c];
                    const det = (i1 - i2) * (j1 - j3) - (i1 - i3) * (j1 - j2);
                    if (Math.abs(det) < EPS) continue;
                    // Solve (i1-i2)x + (j1-j2)y = c2-c1 and (i1-i3)x + (j1-j3)y = c3-c1
                    const r1 = c2 - c1, r2 = c3 - c1;
                    const x = (r1 * (j1 - j3) - r2 * (j1 - j2)) / det;
                    const y = ((i1 - i2) * r2 - (i1 - i3) * r1) / det;
                    const value = c1 + i1 * x + j1 * y;
                    if (Math.abs(evaluate(A, x, y) - value) > EPS) continue;
                    // A fourth maximal monomial means a non-smooth vertex
                    let maxCount = 0;
                    for (const [i, j, cc] of A) if (Math.abs(cc + i * x + j * y - value) < EPS) maxCount++;
                    if (maxCount !== 3 || Math.abs(det) !== 1) return null;
                    vertices.push({ x, y, mons: [a, b, c] });
                }
            }
        }
        if (vertices.length !== d * d) return null;

        // Each pair of maximal monomials of a vertex gives an edge:
        // shared by two vertices -> bounded edge, otherwise -> unbounded ray.
        const edgeMap = new Map();
        vertices.forEach((v, k) => {
            const [a, b, c] = v.mons;
            for (const [p, q, r] of [[a, b, c], [a, c, b], [b, c, a]]) {
                const key = p + "," + q;
                if (!edgeMap.has(key)) edgeMap.set(key, []);
                edgeMap.get(key).push({ k, p, q, r });
            }
        });
        const segments = [];
        const rays = [];
        for (const list of edgeMap.values()) {
            if (list.length === 2) {
                segments.push([vertices[list[0].k], vertices[list[1].k]]);
            } else {
                const { k, p, q, r } = list[0];
                const [ip, jp] = A[p], [iq, jq] = A[q], [ir, jr] = A[r];
                // Direction along which p and q stay equal; orient it so that p beats r.
                let dx = -(jp - jq), dy = ip - iq;
                if ((ip - ir) * dx + (jp - jr) * dy < 0) { dx = -dx; dy = -dy; }
                const len = Math.hypot(dx, dy);
                rays.push({ from: vertices[k], dx: dx / len, dy: dy / len });
            }
        }
        return { vertices, segments, rays };
    }

    // "Not too boring": smooth, and no bounded edge tiny compared to the whole curve.
    function randomCurve() {
        let fallback = null;
        for (let t = 0; t < MAX_TRIES; t++) {
            const d = Math.random() < 0.5 ? 3 : 4;
            const curve = computeCurve(randomPolynomial(d), d);
            if (!curve) continue;
            const xs = curve.vertices.map(v => v.x), ys = curve.vertices.map(v => v.y);
            const diag = Math.hypot(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys));
            const minEdge = Math.min(...curve.segments.map(([u, v]) => Math.hypot(u.x - v.x, u.y - v.y)));
            if (minEdge > 0.06 * diag) return curve;
            fallback = fallback || curve;
        }
        return fallback;
    }

    function draw(svg, curve) {
        const W = window.innerWidth, H = window.innerHeight;
        svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
        svg.replaceChildren();

        const xs = curve.vertices.map(v => v.x), ys = curve.vertices.map(v => v.y);
        const minX = Math.min(...xs), maxX = Math.max(...xs);
        const minY = Math.min(...ys), maxY = Math.max(...ys);
        const scale = 0.8 * Math.min(W / (maxX - minX || 1), H / (maxY - minY || 1));
        const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
        const toX = x => W / 2 + (x - cx) * scale;
        const toY = y => H / 2 - (y - cy) * scale; // SVG y-axis points down

        let path = "";
        for (const [u, v] of curve.segments) {
            path += `M${toX(u.x)},${toY(u.y)}L${toX(v.x)},${toY(v.y)}`;
        }
        const far = 2 * (W + H);
        for (const { from, dx, dy } of curve.rays) {
            const x0 = toX(from.x), y0 = toY(from.y);
            path += `M${x0},${y0}L${x0 + far * dx},${y0 - far * dy}`;
        }
        const p = document.createElementNS("http://www.w3.org/2000/svg", "path");
        p.setAttribute("d", path);
        p.setAttribute("fill", "none");
        p.setAttribute("stroke", STROKE);
        p.setAttribute("stroke-width", STROKE_WIDTH);
        p.setAttribute("stroke-linejoin", "round");
        svg.appendChild(p);
    }

    const curve = randomCurve();
    if (!curve) return;
    const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    svg.setAttribute("class", "bg-curve");
    svg.setAttribute("aria-hidden", "true");
    document.body.prepend(svg);
    draw(svg, curve);

    let pending = false;
    window.addEventListener("resize", () => {
        if (pending) return;
        pending = true;
        requestAnimationFrame(() => { pending = false; draw(svg, curve); });
    });
})();
