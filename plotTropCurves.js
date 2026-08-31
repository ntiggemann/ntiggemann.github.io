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

/*
The initial LaTeX code: Color declarations for easier changes
*/
let latexCode = "";
let latexCodeDS = "";

/*
eps - the tolerance how close to numbers must be to be considered equal
*/
const eps = 0.00000000000001;

/*
Compuations & helpers for that
*/

function parsePolynomialTrop(polyStr, useMin) {
    /*
     * Parses a tropical polynomial in the variables x and y.
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

    // Temporary map:
    // "xExp,yExp" -> coefficient
    const coefficients = new Map();

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

function getValuesIndices(a, b, value, eps) {
    /**
     * Returns the indices of all entries in two arrays whose value is equal
     * to a given value up to a specified tolerance.
     *
     * @param {number[]} a - First array of values.
     * @param {number[]} b - Second array of values.
     * @param {number} value - Target value.
     * @param {number} eps - Numerical tolerance.
     * @returns {{aValueIndices: number[], bValueIndices: number[]}}
     *          Object containing the matching indices for both arrays.
     */
    const aValueIndices = [];
    const bValueIndices = [];

    // Find all indices in the first array whose value matches the target.
    for (let i = 0; i < a.length; i++) {
        if (Math.abs(a[i] - value) < eps) {
            aValueIndices.push(i);
        }
    }

    // Find all indices in the second array whose value matches the target.
    for (let i = 0; i < b.length; i++) {
        if (Math.abs(b[i] - value) < eps) {
            bValueIndices.push(i);
        }
    }

    return {
        aValueIndices,
        bValueIndices
    };
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

function commonSublists(list1, list2) {
    /**
     * Finds all 2-dimensional integer arrays that are present in both input arrays.
     *
     * @param {number[][]} list1 - First array of 2-element arrays.
     * @param {number[][]} list2 - Second array of 2-element arrays.
     * @returns {number[][]} An array containing all common 2-element arrays.
     */

    // Convert the first list into a set of strings of the form "a,b".
    const set1 = new Set(list1.map(sublist => `${sublist[0]},${sublist[1]}`));

    // Convert the second list into a set of strings.
    const set2 = new Set(list2.map(sublist => `${sublist[0]},${sublist[1]}`));

    const common = [];

    // Iterate over the elements of the first set and keep those
    // that also occur in the second set.
    for (const key of set1) {
        if (set2.has(key)) {
            // Convert the string back into a pair of numbers.
            common.push(key.split(",").map(Number));
        }
    }

    return common;
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

function containsArray(arrayOfArrays, target) {
/*     Checks if target is an element of the Array of arrays */
    return arrayOfArrays.some(
        arr =>
            arr.length === target.length &&
            arr.every((value, i) => value === target[i])
    );
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

function plotTropPolynomial(polyStr,aX,bX,aY,bY,useMin,autoAdjust,autoBoundaryDist,autoSquare){

    // Get polynomial
    let A = parsePolynomialTrop(polyStr,useMin);

    // Delete error messages
    showError("");

    // Catch case of only one monomial
    if (A.length < 2){
        showError("Please enter at least two different monomials.");
        return;
    }
    
    if (useMin) {
        A = changeSignsOfPolynomial(A);
        let temp = aX
        aX = (-1)*bX;
        bX = (-1)*temp;
        temp = aY;
        aY = (-1)*bY;
        bY = (-1)*temp;
    }
    // Compute the vertices of the curve, aka triple max points
    let xValsTriple = []; // List of the x values for triple-max points
    let yValsTriple = []; // List of the y values for triple-max points
    let isTriple = []; // List of length 3 arrays, the x-exponents of the maximal monomials
    let jsTriple = []; // List of length 3 arrays, the y-exponents of the maximal monomials

    // Walk through the coefficients A, check all triples of monomials
    for (let a = 0; a < A.length; a++) {
        for (let b = a + 1; b < A.length; b++) {
            for (let c = b + 1; c < A.length; c++) {

                const [i1, j1, c1] = A[a];
                const [i2, j2, c2] = A[b];
                const [i3, j3, c3] = A[c];

                //// Solve for points where a value is attained thrice
                
                // The coeffs of the lineqs we have to solve
                const B = [[i1-i2, j1-j2],
                            [i1-i3, j1-j3]];
                // proceed only if the system of eqs is uniquely solvable
                if (!(Math.abs((i1-i2)*(j1-j3) - (i1-i3)*(j1-j2)) < eps)) {
                    const b = [c2-c1,c3-c1];
                    const x = solve2x2(B,b);
                    // Check if it is actually a maximum
                    if (Math.abs(evaluateTropPolynomial(x[0],x[1],A) - c1 - i1*x[0] - j1*x[1])<eps){
                        xValsTriple.push(x[0]);
                        yValsTriple.push(x[1]);
                        isTriple.push([i1,i2,i3]);
                        jsTriple.push([j1,j2,j3]);
                    }
                }
            }
        }
    }
    // There exist vertices:
    if (xValsTriple.length > 0) {
        if (autoAdjust){
            aX = Math.min(...xValsTriple);
            bX = Math.max(...xValsTriple);
            aY = Math.min(...yValsTriple);
            bY = Math.max(...yValsTriple);
        } else {
            if ((aX >= bX) || (aY >= bY)){

console.log(aX,bX,aY,bY);
                showError("Invalid boundaries.");
                return;
            }

            let missingPoints = []; // Store tripple max points not contained in plotting range
            for (let i = 0; i < xValsTriple.length; i++){
                let pt = [xValsTriple[i],yValsTriple[i]];
                // If vertex not in plot
                if ((pt[0] <= aX) || (pt[0] >= bX) || (pt[1] <= aY) || (pt[1] >= bY)){
                    missingPoints.push(pt);
                }
            }

            if (missingPoints.length > 0){
                aX = Math.min(...xValsTriple) - 1;
                bX = Math.max(...xValsTriple) + 1;
                aY = Math.min(...yValsTriple) - 1;
                bY = Math.max(...yValsTriple) + 1;

                showError(`All vertices must be in the plotting range. Suggested range: [${aX},${bX}] x [${aY},${bY}]`);
                return;
            }
        }
    } else { // We only have parallel lines
        if (autoAdjust){
            let pseudoVertices = [];
            for (let a = 0; a < A.length; a++) {
                for (let b = a + 1; b < A.length; b++) {
                    // The equation where they agree is given by (i1-i2)x+(j1-j2)y+c1-c2 = 0
                    const [i1, j1, c1] = A[a];
                    const [i2, j2, c2] = A[b];
                    // If i1-i2 not= 0: one point is: (...,0)
                    if (i1-i2 !== 0){
                        let xValue = (c2-c1)/(i1-i2);
                        let yValue = 0;
                        // If its a point of the "curve"
                        if ((Math.abs(evaluateTropPolynomial(xValue,yValue,A) - c1 - i1*xValue) < eps)){
                            pseudoVertices.push([xValue,yValue]);
                        }
                    } else {
                        // Else: (0,...)
                        let xValue = 0;
                        let yValue = (c2-c1)/(j1-j2);
                        // If its a point of the "curve"
                        if ((Math.abs(evaluateTropPolynomial(xValue,yValue,A) - c1 - j1*yValue)) < eps){
                            pseudoVertices.push([xValue,yValue]);
                        }
                    }
                }
            }
            // Min/Max for x values
            aX = Math.min(...pseudoVertices.map(x => x[0]));
            bX = Math.max(...pseudoVertices.map(x => x[0]));
            // Min/Max for y values
            aY = Math.min(...pseudoVertices.map(x => x[1]));
            bY = Math.max(...pseudoVertices.map(x => x[1]));
        }
    }

    if (autoAdjust && autoSquare){
        const size = Math.max(bX - aX, bY - aY);
        const aXNew = ((aX+bX)/2) - (size/2) - autoBoundaryDist;
        const bXNew = ((aX+bX)/2) + (size/2) + autoBoundaryDist;
        const aYNew = ((aY+bY)/2) - (size/2) - autoBoundaryDist;
        const bYNew = ((aY+bY)/2) + (size/2) + autoBoundaryDist;

        aX = aXNew;
        bX = bXNew;
        aY = aYNew;
        bY = bYNew;
    } else if (autoAdjust) {
        aX = aX - autoBoundaryDist;
        bX = bX + autoBoundaryDist;
        aY = aY - autoBoundaryDist;
        bY = bY + autoBoundaryDist;
    }
    // Compute the intersection points of the curve with the boundary of the plot
    // Line 383

    let xValsBdry = [];// List of the x values for points in the trop variety on the boundary
    let yValsBdry = [];// List of the y values for points in the trop variety on the boundary
    let isBdry = [];// List of length 2 arrays, the x-exponents of the maximal monomials
    let jsBdry = [];// List of length 2 arrays, the y-exponents of the maximal monomials

    // Check all pairs of monomials
    for (let a = 0; a < A.length; a++) {
        for (let b = a + 1; b < A.length; b++) {
            const [i1, j1, c1] = A[a];
            const [i2, j2, c2] = A[b];

            // Solve for points on boundary where a value is attained twice
            const rhsAX = c2 + (i2 - i1)* aX - c1;
            const rhsBX = c2 + (i2 - i1)* bX - c1;
            const rhsAY = c2 + (j2 - j1)* aY - c1;
            const rhsBY = c2 + (j2 - j1)* bY - c1;

            // If equation solvable
            if (Math.abs(j1-j2)>eps) {
                // If max actually attained
                let valBdry = (1/(j1-j2))*rhsAX;
                if (Math.abs(evaluateTropPolynomial(aX,valBdry,A) - c1 - i1*aX - j1*valBdry)<eps){
                    // if point in plotting range
                    if ((valBdry >= aY) && (valBdry <= bY)){
                        xValsBdry.push(aX);
                        yValsBdry.push(valBdry);
                        isBdry.push([i1,i2]);
                        jsBdry.push([j1,j2]);
                    }
                }
                valBdry = (1/(j1-j2))*rhsBX;
                if (Math.abs(evaluateTropPolynomial(bX,valBdry,A) - c1 - i1*bX - j1*valBdry)<eps){
                    // if point in plotting range
                    if ((valBdry >= aY) && (valBdry <= bY)){
                        xValsBdry.push(bX);
                        yValsBdry.push(valBdry);
                        isBdry.push([i1,i2]);
                        jsBdry.push([j1,j2]);
                    }
                }
            }

            // If equation solvable
            if (Math.abs(i1-i2)>eps) {
                // If max actually attained
                let valBdry = (1/(i1-i2))*rhsAY;
                if (Math.abs(evaluateTropPolynomial(valBdry,aY,A) - c1 - i1*valBdry - j1*aY)<eps){
                    // if point in plotting range
                    if ((valBdry >= aX) && (valBdry <= bX)){
                        xValsBdry.push(valBdry);
                        yValsBdry.push(aY);
                        isBdry.push([i1,i2]);
                        jsBdry.push([j1,j2]);
                    }
                }
                valBdry = (1/(i1-i2))*rhsBY;
                if (Math.abs(evaluateTropPolynomial(valBdry,bY,A) - c1 - i1*valBdry - j1*bY)<eps){
                    // if point in plotting range
                    if ((valBdry >= aX) && (valBdry <= bX)){
                        xValsBdry.push(valBdry);
                        yValsBdry.push(bY);
                        isBdry.push([i1,i2]);
                        jsBdry.push([j1,j2]);
                    }
                }
            }
        }
    }

    //// Now we have all points at which edges start and end
    // 437ff
    // Now we check which points we should connect with a line

    allVertices = [];// List of all vertices: Tuples [x,y] of coordinates of vertices of the curve (including boundarys of the plot)
    allVerticesExp = [];// And a list of lists of their exponents [i,j], all listed only one time

    for (let i = 0; i < xValsTriple.length; i++) {
        const indexOfTriplePoint = tupleIsInList([xValsTriple[i], yValsTriple[i]],allVertices,eps);
        if (indexOfTriplePoint === -1) {
            allVertices.push([xValsTriple[i], yValsTriple[i]]);
            allVerticesExp.push([[isTriple[i][0],jsTriple[i][0]],[isTriple[i][1],jsTriple[i][1]],[isTriple[i][2],jsTriple[i][2]]]);
        } else {
            for (let d = 0; d < 3; d++) {
                if (!containsArray(allVerticesExp[indexOfTriplePoint],[isTriple[i][d],jsTriple[i][d]])) {
                    allVerticesExp[indexOfTriplePoint].push([isTriple[i][d],jsTriple[i][d]]);
                }
            } 
        }
    }

    for (let i = 0; i < xValsBdry.length; i++) {
        const indexOfBdryPoint = tupleIsInList([xValsBdry[i], yValsBdry[i]],allVertices,eps);
        if (indexOfBdryPoint === -1) {
            allVertices.push([xValsBdry[i], yValsBdry[i]]);
            allVerticesExp.push([[isBdry[i][0],jsBdry[i][0]],[isBdry[i][1],jsBdry[i][1]]]);
        } else {
            for (let d = 0; d < 2; d++) {
                if (!containsArray(allVerticesExp[indexOfBdryPoint],[isBdry[i][d],jsBdry[i][d]])) {
                    allVerticesExp[indexOfBdryPoint].push([isBdry[i][d],jsBdry[i][d]]);
                }
            } 
        }
    }

    //// Check which vertices to connect 
    // 469ff
    const numberOfVertices = allVertices.length; // Number of vertices
    let xEdgeVals = []; // List of the x values, as length 2 arrays, of start and endpoint of edges of the trop variety
    let yEdgeVals = []; // List of the y values, as length 2 arrays, of start and endpoint of edges of the trop variety 
    let edgeWeights = []; // The weights of the edges

    let xValsDual = []; // List of the x values of start and endpoint of edges in the dual
    let yValsDual = []; // List of the y values of start and endpoint of edges in the dual

    // Go through pairs of vertices
    let edgeCount = 0;
    for (let i = 0; i < numberOfVertices; i++) {
        for (let j = i + 1; j < numberOfVertices; j++) {
            // If the points should get connected
            const commonExps = commonSublists(allVerticesExp[i],allVerticesExp[j]);
            if (commonExps.length > 1){
                // Add the start and end points of the edge
                xEdgeVals.push([allVertices[i][0],allVertices[j][0]]);
                yEdgeVals.push([allVertices[i][1],allVertices[j][1]]);

                // And its dual + weight
                // The exponents which give the weight are the ones that define the dual edge. So compute all possible weights.
                edgeWeights.push(0);
                xValsDual.push([]);
                yValsDual.push([]);
                for (let a = 0; a < commonExps.length;a++){
                    for (let b = a+1; b < commonExps.length;b++){
                        const [i1, j1] = commonExps[a];
                        const [i2, j2] = commonExps[b];
                        const w = gcd(i1-i2,j1-j2);
                        if (edgeWeights[edgeCount] < w){
                            edgeWeights[edgeCount] = w;
                            xValsDual[edgeCount] = [i1,i2];
                            yValsDual[edgeCount] = [j1,j2];
                        }
                    }
                }
                edgeCount++;
            }
        }
    }
    //// Computations for the dual
    let maxX = 0;
    let maxY = 0;
    for (let a = 0; a < A.length; a++){
        const [i1,j1,c1] = A[a];
        if (i1>maxX){
            maxX = i1;
        }
        if (j1>maxY){
            maxY = j1;
        }
    }
    // Then we can create the coordinates for the points where we draw the DS
    let xCoordsDual = [];
    let yCoordsDual = [];

    for (let i = 0; i < maxX + 1; i++){
        for (let j = 0; j < maxY + 1; j++){
            xCoordsDual.push(i);
            yCoordsDual.push(j);
        }
    }
    //// Plotting

    //// Plotting the curve
    const svg = document.getElementById("plot");

    svg.replaceChildren();// Clears the plot

    if (useMin) {
        A = changeSignsOfPolynomial(A);
        let temp = aX
        aX = (-1)*bX;
        bX = (-1)*temp;
        temp = aY;
        aY = (-1)*bY;
        bY = (-1)*temp;

        for (let i = 0; i < xEdgeVals.length; i++){
            xEdgeVals[i][0] = (-1)*xEdgeVals[i][0];
            xEdgeVals[i][1] = (-1)*xEdgeVals[i][1];
            yEdgeVals[i][0] = (-1)*yEdgeVals[i][0];
            yEdgeVals[i][1] = (-1)*yEdgeVals[i][1];
        }
    }

    // 530ff: Draw Curve
    for (let i = 0; i<xEdgeVals.length; i++){
        let line = document.createElementNS(
        "http://www.w3.org/2000/svg",
        "line"
        );

        line.setAttribute("x1",xSVGCoords(xEdgeVals[i][0],aX,bX));
        line.setAttribute("y1",ySVGCoords(yEdgeVals[i][0],aY,bY));
        line.setAttribute("x2",xSVGCoords(xEdgeVals[i][1],aX,bX));
        line.setAttribute("y2",ySVGCoords(yEdgeVals[i][1],aY,bY));
        line.setAttribute("stroke", "black");
        line.setAttribute("stroke-width", 0.02);

        svg.appendChild(line);

        latexCode = latexCode + `\n\\draw[\\colorCurve] (${xEdgeVals[i][0]},${yEdgeVals[i][0]}) -- (${xEdgeVals[i][1]},${yEdgeVals[i][1]});`
        
        

        if (edgeWeights[i] > 1){            
            let weight = document.createElementNS(
                "http://www.w3.org/2000/svg",
                "text"
            );
            let xCoordWeight = (xEdgeVals[i][0] + xEdgeVals[i][1]) / 2;
            let yCoordWeight = (yEdgeVals[i][0] + yEdgeVals[i][1]) / 2;

            weight.setAttribute("x", xSVGCoords(xCoordWeight,aX,bX));
            weight.setAttribute("y", ySVGCoords(yCoordWeight,aY,bY));
            
            weight.textContent = `${edgeWeights[i]}`;
            weight.setAttribute("font-size", "0.5");

            svg.appendChild(weight);

            latexCode = latexCode + `\n\\node[\\colorWeights] at (${xCoordWeight},${yCoordWeight}) {${edgeWeights[i]}};`
        }
    }

    
    //// Plotting the dual
    const svgDS = document.getElementById("plotDS");

    svgDS.replaceChildren();// Clears the plot
    sizeDS = Math.max(maxX,maxY);

    if (sizeDS > maxY){
        for (let i = 0; i < yCoordsDual.length;i++){
            yCoordsDual[i] = yCoordsDual[i] - maxY + maxX;
        }
        for (let i = 0; i < yValsDual.length;i++){
            yValsDual[i][0] = yValsDual[i][0] - maxY + maxX;
            yValsDual[i][1] = yValsDual[i][1] - maxY + maxX;
        }
    }
    // Draw lines dual
    for (let i=0; i < xValsDual.length; i++) {
        let line = document.createElementNS(
        "http://www.w3.org/2000/svg",
        "line"
        );

        line.setAttribute("x1",xSVGCoords(xValsDual[i][0],-0.15,sizeDS+0.15));
        line.setAttribute("y1",ySVGCoords(yValsDual[i][0],-0.15,sizeDS+0.15));
        line.setAttribute("x2",xSVGCoords(xValsDual[i][1],-0.15,sizeDS+0.15));
        line.setAttribute("y2",ySVGCoords(yValsDual[i][1],-0.15,sizeDS+0.15));
        line.setAttribute("stroke", "black");
        line.setAttribute("stroke-width", 0.02);

        svgDS.appendChild(line);
        
        latexCodeDS = latexCodeDS + `\n\\draw[\\colorDSLines] (${xValsDual[i][0]},${yValsDual[i][0]}) -- (${xValsDual[i][1]},${yValsDual[i][1]});`
        
    }
    // Draw all points of \Delta_deg
    for (let i = 0; i < xCoordsDual.length; i++) {
        let circle = document.createElementNS(
            "http://www.w3.org/2000/svg",
            "circle"
        );

        circle.setAttribute("cx", xSVGCoords(xCoordsDual[i],-0.15,sizeDS+0.15));
        circle.setAttribute("cy", ySVGCoords(yCoordsDual[i],-0.15,sizeDS+0.15));
        circle.setAttribute("r", 0.15);
        circle.setAttribute("fill", "black");

        svgDS.appendChild(circle);
        
        latexCodeDS = latexCodeDS + `\n\\filldraw[color = \\colorDSDots, fill=\\colorDSDots] (${xCoordsDual[i]},${yCoordsDual[i]}) circle (1.5pt);`
    }
}

/*
Plot the computed curve
*/

function plotte(){
    const polyStr = document.getElementById("poly").value;
    const aX = document.getElementById("aX").valueAsNumber;
    const bX = document.getElementById("bX").valueAsNumber;
    const aY = document.getElementById("aY").valueAsNumber;
    const bY = document.getElementById("bY").valueAsNumber;
    const useMin = !document.getElementById("useMax").checked;
    const autoAdjust = document.getElementById("autoAdjust").checked;
    const autoBoundaryDist = document.getElementById("autoBoundaryDist").valueAsNumber;
    const autoSquare = document.getElementById("autoSquare").checked;
    
    // Reset the initial LaTeX code: Color declarations for easier changes
    latexCode = "\\def\\colorCurve{black}\n\\def\\colorWeights{black}";
    latexCodeDS = "\\def\\colorDSLines{black}\n\\def\\colorDSDots{black}";

    plotTropPolynomial(polyStr,aX,bX,aY,bY,useMin,autoAdjust,autoBoundaryDist,autoSquare);
}

/*
Functions that handle the websites behaviour
*/

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

