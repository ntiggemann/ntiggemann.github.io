///
This script checks the polytopes in listSmooth3D and listSmooth2D for combinatorial orientability. Further, it does the computations for the example on the blow up of P^3
///

///
Import packages
///

needsPackage "LatticePolytopes";

///
Define some functions
///
-- Function to compute lattice points on boundary of a polyhedron
-- Input: Polyhedron P
-- Output: The lattice points on the boundary as a list of 3x1 matrices
boundaryLatticePoints = P -> (
    inner := interiorLatticePoints(P);
    latPts := latticePoints(P);
    return select(latPts, a -> not member(a,inner));
);

-- Function to compute matrix mod 2
-- Input: Matrix M with integer entries
-- Output: Matrix with entries defined by the entries of M modulo 2
-- Note: Output still matrix over same ring as input
mod2Matrix = M -> matrix apply(entries M, row -> apply(row, x -> x % 2));

-- Function to compute matrix-mod-2 image of list of matrices as a set
-- Input: List of matrices
-- Output: Set of images of the matrices under taking mod 2
-- Note: Output still matrices over same ring as input
mod2Set = L -> (set apply(L, mod2Matrix));

-- Checks a necessary condition for comb orient - Only on partial P.
-- Input: Polyhedron P
-- Output: Boolean: True iff the polyhedron satisfies: the image of the coordinates of the lattice points on the boundary mod 2 is NOT whole (Z/2Z)^n
couldBeCombOrient = P -> (
	(#(mod2Set(boundaryLatticePoints(P))) != 2^(ambDim P))
);

-- Function to get the vertices of a polyhedron as a list of nx1 Matrices
-- I.p. the output is like at the other get-points-commands; this saves some headaches
-- Input: Polyhedron P
-- Output: List of nx1 matrices that are the vertices of P
verticesMat = P -> (
    apply(entries transpose substitute(vertices(P),ZZ), c -> transpose matrix{c})
);

-- In our context: Given vertex i, returns the edges adjacent to vertex i
-- Input: Matrix M and integer i
-- Output: List of indices for which M_(i,j) = 1 and j>0.
oneColumns = (M,i) -> (
    select(toList(1 .. numColumns(M)-1), j -> M_(i,j) == 1)
);

-- In our context: Searches the other (i.e. != i) vertex adjacent to edge j
-- Input: Matrix M, indices i,j. 
-- Output: The index k other than i in column j with entry 1 and k > 0
otherOneRow = (M,i,j) -> (
    for k from 1 to (numRows(M)-1) do (
        if i!=k and M_(k,j) == 1 then (
            return k;
        )
    );
    return null;
);

-- Input: Integer vector as nx1 matrix. 
-- Output: This vector divided by the gcd of its entries, i.p. a primitive vector
primitiveVector = v -> (
    g := gcd flatten entries v;
    v // g
);

-- Function: Checks if for a smooth polytope one of its NP_sigma has an odd boundary point
-- Note: This function does NOT check for smoothness; to check if P is smooth, use: isSmooth(normalFan(P)) 
-- Input: A smooth polytope 
-- Output: A boolean, true iff all NP_sigma have NO odd boundary point
hasNoOddBdrySigmas = testPolyh -> (
    -- Make the faster first check for the necessary condition
    if not couldBeCombOrient(testPolyh) then (
        return false;
    );
    -- Get the vertices
    testPolyhVerts := verticesMat(testPolyh);
    -- Get the data which vertices share an edge
    testPolyhEdgeMatrix := vertexEdgeMatrix(testPolyh);
    -- Define the all-one-vector
    ones := matrix apply(ambDim(testPolyh), i -> {1});
    -- Iterate over the vertices
    for i from 0 to #testPolyhVerts-1 do (
        -- This safes some headaches, as vertices are counted starting with 1 in vertexEdgeMatrix but indices generally start at 0
        vertexInd := i+1;
        -- Get the edges containing vertex #vertexInd
        edgeInds := oneColumns(testPolyhEdgeMatrix,vertexInd);
        -- Empty list to store basis vectors in
        basisVectors := {};
        -- Compute basis vectors
        for edgeInd in edgeInds do (
                --Get the other vertex
                otherVertex := testPolyhVerts#(otherOneRow(testPolyhEdgeMatrix,vertexInd,edgeInd)-1);
                -- Compute difference and normalize, i.e. compute the basis vector in this direction
                basisVectors = append(basisVectors,primitiveVector(otherVertex-testPolyhVerts#i));
        );
        -- Compute the basechange matrix: Write the basisvectors as columns and invert the matrix
        testPolyhBaseChange := inverse(matrix{basisVectors});
        -- Compute the shift of the boundary, shift the vertex to 0
        shiftedPolyhBdry := apply(boundaryLatticePoints(testPolyh), A -> (A-testPolyhVerts#i));
        -- Transform the boundary by applying basechange matrix
        bdryNPSigma := apply(shiftedPolyhBdry, A -> (testPolyhBaseChange * A));
        -- Compute image modulo 2 as a set
        bdryNPSigmaMod2 := mod2Set(bdryNPSigma);
        -- Check if the image contains (1,...,1), i.e. the polytope relative to the corner has an odd boundary point
        if (member(ones,bdryNPSigmaMod2)) then (
                return false;
        );
    );
    return true;
);

-- Function: Computes all points of (Z/2Z)^n
-- Input: A natural number
-- Output: A list containing all points of (Z/2Z)^n as n x 1 matrices over ZZ/2
allPointsF2n = n -> (
    apply(subsets(toList(0..n-1)), S ->
        substitute(matrix apply(n, i -> if member(i,S) then {1} else {0}),ZZ/2)
    )
);

-- Function: Checks if a polytope is comb orient
-- Input: A polytope 
-- Output: A boolean, true iff the polytope is combinatorially orientable
isCombOrient = P -> (
    n := ambDim(P);
    -- Two dimensional case is easy, in dim 2 checking boundary points suffices
    if n == 2 then (
        return couldBeCombOrient(P);
    );
    if couldBeCombOrient(P) then (
        -- Get the supporting hyperspaces
        -- Gets sequence of length two: 0th entry is coeff matrix of equalities, 1st entry is vector with the affine part 
        -- Then take these entries mod 2
        (eqMod2,affMod2) := apply(ambientHalfspaces(P), A -> (substitute(A,ZZ/2)));
        -- Get number of equations
        m := numRows(eqMod2);
        -- Get all points of (Z/2Z)^n
        pnts := allPointsF2n(n);
        -- Iterate over the points
        for pt in pnts do (
            B := eqMod2*pt;
            -- If a point does not solve any of the equations, we are done
            if not isMember(true,apply(m, i -> (B_(i,0) == affMod2_(i,0)))) then (
                return true;
            );
        );
    );
    -- All points solved at least one equation or already a problem for boundary points
    return false;
);

///
Computations
///

///
For small two dimensional smooth lattice polygons
///
-- Get list of smooth 2D polytopes
smoothies = listSmooth2D();
-- Initiate counters
combOrientCount = 0;
notCombOrientCount = 0;
for testPolyh in smoothies do (
    if couldBeCombOrient(testPolyh) then (
        combOrientCount = combOrientCount+1;
    ) else  (
        notCombOrientCount = notCombOrientCount+1;
    );
);
print("Of the checked 41 smooth two dimensional polytopes");
print(combOrientCount);
print("were combinatorially oriented and");
print(notCombOrientCount);
print("were not.");



///
For small three dimensional smooth lattice polytopes
///

-- Get list of smooth 3D polytopes
smoothies = listSmooth3D();
-- Initiate counters
combOrientCount = 0;
notCombOrientCount = 0;
noOddBdryCount = 0;
couldBeCombOrientCount = 0;
for testPolyh in smoothies do (
    if isCombOrient(testPolyh) then (
        combOrientCount = combOrientCount+1;
    ) else  (
        notCombOrientCount = notCombOrientCount+1;
    );
    if hasNoOddBdrySigmas(testPolyh) then  (
        noOddBdryCount = noOddBdryCount+1;
    );
    if couldBeCombOrient(testPolyh) then (
        couldBeCombOrientCount = couldBeCombOrientCount+1;
    );
);

print("Of the checked 103 smooth three dimensional polytopes");
print(combOrientCount);
print("were combinatorially oriented and");
print(notCombOrientCount);
print("were not,");
print(noOddBdryCount);
print("have no odd point on the boundary for every vertex and for");
print(couldBeCombOrientCount);
print("the boundary does not have all parities.");

///
For the example
///
print("For the polytope of the divisor on the blowup of P^3:");
-- The polytope
P = convexHull matrix {
    {1,0,0,4,0,0},
    {0,1,0,0,4,0},
    {0,0,1,0,0,4}
};

-- Check that smooth
if not isSmooth(normalFan(P)) then (
    error "No smoothy :(";
);
-- Get the vertices
PVerts := verticesMat(P);
-- Get the data which vertices share an edge
PEdgeMatrix := vertexEdgeMatrix(P);
-- Define the all-one-vector
ones := matrix apply(ambDim(P), i -> {1});
-- Iterate over the vertices
for i from 0 to #PVerts-1 do (
    -- This safes some headaches, as vertices are counted starting with 1 in vertexEdgeMatrix but indices generally start at 0
    vertexInd := i+1;
    -- Get the edges containing vertex #vertexInd
    edgeInds := oneColumns(PEdgeMatrix,vertexInd);
    -- Empty list to store basis vectors in
    basisVectors := {};
    -- Compute basis vectors
    for edgeInd in edgeInds do (
            --Get the other vertex
            otherVertex := PVerts#(otherOneRow(PEdgeMatrix,vertexInd,edgeInd)-1);
            -- Compute difference and normalize, i.e. compute the basis vector in this direction
            basisVectors = append(basisVectors,primitiveVector(otherVertex-PVerts#i));
    );
    -- Compute the basechange matrix: Write the basisvectors as columns and invert the matrix
    PBaseChange := inverse(matrix{basisVectors});
    -- Compute the shift of the boundary, shift the vertex to 0
    shiftedPolyhBdry := apply(boundaryLatticePoints(P), A -> (A-PVerts#i));
    -- Transform the boundary by applying basechange matrix
    bdryNPSigma := apply(shiftedPolyhBdry, A -> (PBaseChange * A));
    -- Compute image modulo 2 as a set
    bdryNPSigmaMod2 := mod2Set(bdryNPSigma);
    -- Check if the image contains (1,...,1), i.e. the polytope relative to the corner has an odd boundary point
    hasOddBdry = member(ones,bdryNPSigmaMod2);

    -- The Output:
    print("For vertex");
    print(PVerts#i);
    print("there exists an odd boundary point:");
    print(hasOddBdry);
);
print("This polytope does not have all parities in its boundary:");
print(couldBeCombOrient(P));
print("This polytope is combinatorially orientable:");
print(isCombOrient(P));