###############################
# This script plots the tropical curve defined by a trivariate tropical polynomial
# The tropical numbers are defined with max as addition
# v1.1
###############################

''' Preamble '''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from matplotlib.ticker import MaxNLocator

infty = -9999999999999999 # Auxiliar-minus-infty - just some very negative number

''' INPUT '''
# This is the only stuff you should change

''' Define the polynomial '''
# Highest occuring exponent of a variable (e.g. for x^3y^2 + y^2 set d = 3)
# Setting this to sth like 10 and just forgetting about it doesn't really harm the computation time.
d = 10

# A[i,j,k] is the coefficient that belongs to x^iy^jz^k
A = infty * np.ones((d+1,d+1,d+1))# Initialising A - this sould not be deleted

try:    
    A[0,0,0] = -2
    A[0,1,0] = 0
    A[1,0,0] = 0
    A[2,0,0] = -1
    A[0,2,0] = -1
    A[1,1,0] = 1
    A[2,1,0] = -1
    A[1,2,0] = -1
    A[3,0,0] = -4
    A[0,3,0] = -4
    A[0,0,1] = 0
    
    '''
    A[0,0,0] = 1
    A[0,0,1] = -1
    A[0,0,2] = 0
    A[0,0,3] = 0
    '''
    '''
    A[0,0,0] = 0
    A[1,0,0] = 0
    '''
    '''
    A[2,0,0] = 0
    A[0,2,0] = 0
    A[0,0,2] = 0
    A[0,0,0] = 0
    A[1,0,0] = 0
    '''
    '''
    A[0,0,0] = 0
    A[1,0,0] = 0
    A[0,1,0] = 0
    A[0,0,1] = 0
    '''
    '''
    A[0,0,0] = 0
    A[1,0,1] = 0
    A[1,1,0] = 0
    A[2,1,1] = 0 # I = K = 0; J = [1,1,0], L = [1,0,1]
    '''
except:
    raise ValueError("Your highest occuring degree of a variable, d, is too small. (Or you have a syntax error)")

''' Plotting range '''
# As long as you don't do stuff near infty you should be fine

a_x = -6  # Boundary x-axis left
b_x = 6   # Boundary x-axis right
a_y = -6  # Boundary y-axis left
b_y = 6   # Boundary y-axis right
a_z = -6
b_z = 6

''' If something weird happens, read this:'''
## Non-symbolic computing problems: One does not simply check floats for equality
# If you don't do sth weird, this part should not concern you. Some words:
# This eps can usually be chosen very small (~10^{-14}). If the polynomial is of very high degree you might need a bigger eps
# So if the sanity check of a plot fails, this is the first thing to change and check again.

eps = 0.00000000000001  # Epsilon - because testing for equality of floats is bad

''' Code - functions '''
# Don't change stuff after here - unless you really know what you are doing

def linToExp(i):
    '''
        Args:
            Integer between 0 and (d+1)^3
        Returns:
            The exponent of the i-th monomial
    '''
    return (i // (d+1)**2), ((i // (d+1))%(d+1)), (i % (d+1))

def tropPolySet(x,y,z,A,inds):
    '''
    Args:
        x (int of float): x-coordinate
        y (int of float): y-coordinate
        z (int of float): z-coordinate
        A (3-dim Array of int or float): The coefficientmatrix of a tropical polynomial
        inds (list of ont): The entries of A not infty
    Returns:
        array of floats: The values in the max of the tropical polynomial, has length len(inds)
    '''
    n = len(inds)
    v = np.zeros(n)
    for i in range(n):
        i1,j1,k1 = linToExp(inds[i])
        v[i] = i1*x + j1*y + k1*z + A[i1,j1,k1]
    return v

def threeEqual(a,b,val,eps):
    '''
    Args:
        a(list or array of int or float): First tuple
        b(list or array of int or float): Second tuple
        val (int or float): A value
        eps (int or float): error tolerance
    Returns:
        True if at least three entries of both a and b are closer than eps to val
    '''
    count = 0

    for i in range(0,len(a)):
        if (abs(val-a[i]) < eps):
            count = count + 1
    if count >= 3:
        count = 0
        for i in range(0,len(b)):
            if (abs(val-b[i]) < eps):
                count = count + 1
    else:
        return False
    if count >= 3:
        return True
    else:
        return False

def isSubSet3(a1,a2,a3,b1,b2,b3):
    '''
    Args:
        a1,a2,a3,b1,b2,a3 (list of int or float): The coordinates of two sets of points
    Returns:
        True if {(a1,a2,a3)} is a subset of {(b1,b2,b3)}, else False
    '''
    a = []
    b = []
    for i in range(0,len(a1)):
        a.append(np.array([a1[i],a2[i],a3[i]]))
    for i in range(0,len(b1)):
        b.append(np.array([b1[i],b2[i],b3[i]]))
    for i in range(0,len(a)):
        if not any(np.array_equal(a[i], bi) for bi in b):
            return False
    return True

def bothAreIn(i1,j1,k1,i2,j2,k2,list_is,list_js,list_ks):
    '''
    Args:
        i1,j1,k1,i2,j2,k2 (int): Two exponents of monomials
    Returns:
        true if (i1,j1,k1), (i2,j2,k2) both are in {(list_is[l],list_js[l],list_ks[l])|l}
    '''
    exp1 = [i1,j1,k1]
    exp2 = [i2,j2,k2]
    exps = []
    for i in range(0,len(list_is)):
        exps.append([list_is[i],list_js[i],list_ks[i]])
    if (exp1 in exps) and (exp2 in exps):
        return True
    else:
        return False

def getTheEquals(exp1,exp2):
    '''
        Args:
            Two lists of length 3 containing lists of int of length 3
        Returns:
            The Entries that occur twice
    '''
    sol = []
    for i in range(0,len(exp1)):
        if (exp1[i] in exp2):
            sol.append(exp1[i])
    return sol

def tropMonomVal(x,y,z,exp,A):
    '''
        Args:
            Point (x,y,z), exponent of monomial [i,j,k], coeff matrix A
        Returns:
            Value of the tropical monomial at (x,y,z)
    '''
    return (A[exp[0],exp[1],exp[2]] + exp[0]*x + exp[1]*y + exp[2]*z)

def noDouble(arr, eps):
    '''
    Entfernt doppelte Elemente aus der Eingabeliste basierend auf einer Max-Norm-Toleranz.

    Args:
        arr: List of numpy arrays (length 3)
        eps: float - error tolerance in max norm

    Returns:
        Eine Liste mit einzigartigen Elementen aus `arr`, wobei zwei Elemente als gleich betrachtet 
        werden, wenn ihre Max-Norm-Distanz kleiner als `eps` ist.
    '''
    unique = []
    for vector in arr:
        if not isInSet(vector, unique, eps):
            unique.append(vector)
    return unique


def isInSet(vector, unique_list, eps):
    """
    Überprüft, ob ein Vektor bereits in einer Liste enthalten ist, basierend auf der Max-Norm.

    Args:
        vector: numpy array der Länge 3
        unique_list: Liste von numpy arrays (Länge 3)
        eps: float - Fehler-Toleranz für die Max-Norm

    Returns:
        True, wenn `vector` in `unique_list` enthalten ist (bis auf eine Toleranz von `eps`), 
        ansonsten False.
    """
    for existing_vector in unique_list:
        if np.all(np.abs(vector - existing_vector) < eps):
            return True
    return False

def noDoublePolygon(vert_poly, monoms_poly):
    '''
        Args:
            List of lists of three-dim np arrays, List of Lists of ints (exponents of two of the maximal monoms at each polygon)
        Returns:
            Maximal list of lists of three-dim np arrays, such that no set of three dim arrays occurs twice. For each polygon the exponents of all maximal monoms
            Gets rid of polygons with only two vertices
    '''
    temp_vert_poly = [] # Start with no polygon
    temp_monoms = []  # And no max monoms
    for i in range(len(vert_poly)): # For each polygon: If not already in list, add it to the list
        if (len(vert_poly[i]) > 2): # Actually a polygon, i.e. more than two vertices
            index_eq = notAnyEqual(vert_poly[i],temp_vert_poly)
            if index_eq == -1:
                temp_vert_poly.append(vert_poly[i])
                temp_monoms.append(monoms_poly[i])
            else: # If already in list, add the monomials exponents
                temp_monoms[index_eq] = concatListsNoDouble(temp_monoms[index_eq],monoms_poly[i])
    return temp_vert_poly, temp_monoms

def notAnyEqual(vert_poly_i,temp_vert_poly):
    '''
        Args:
            vert_poly_i (list of length 3 arrays of floats): A list of vertics of a polygon
            temp_vert_poly (list of lists of length 3 arrays of floats): A list of the above
        Returns:
            int: -1 iff vert_poly_i not in temp_vert_poly; i iff vert_poly_i in temp_vert_poly, i is the index at which the equality occured
    '''
    for i in range(len(temp_vert_poly)):
        if equalTwoListsOfArr(vert_poly_i, temp_vert_poly[i]):
            return i
    return -1

def concatListsNoDouble(list1, list2):
    '''
        Args:
            Two lists with entries lists of ints of length 3
        Returns:
            List containing all entries of the two lists, but each entry exactly once
    '''
    temp = list2
    for el1 in list1:
        if not(el1 in temp):
            temp.append(el1)
    return temp

def equalTwoListsOfArr(larr1,larr2):
    '''
        Args:
            Two lists of np arrays of length 3
        Returns:
            True if the lists have the same length and larr1 subset larr2
        (In particular: If each list contains no element twice, returns true iff, as sets, larr1 == larr2)
    '''
    if (len(larr1) != len(larr2)):
        return False
    for arr in larr1:
        if not any(np.array_equal(arr,arr2) for arr2 in larr2):
            return False
    return True

def maxMonomsAt(value_poly,eps):
    '''
        Args:
            np.array of doubles
        Returns:
            List of exponents of max monoms
    '''
    exps = []
    themax = np.max(value_poly)
    for i in range(len(value_poly)):
        if (abs(themax-value_poly[i])<eps):
            exps.append([i // (d+1)**2,(i // (d+1))%(d+1),i % (d+1)])
    return exps

def maxIndices(vals,eps):
    '''
        Args:
            np array of doubles, double
        Returns:
            The indices of the maximal (up to eps) values of vals
    '''
    mval = np.max(vals)
    maxinds = []
    for i in range(len(vals)):
        if (abs(mval - vals[i])<eps):
            maxinds.append(i)
    return maxinds

def pointsUnEqual(x1,y1,z1,x2,y2,z2,eps):
    '''
        Args:
            x1,y1,z1,x2,y2,z2 (float): The entries of two points in 3d space
            eps (float): error tolerance
        Returns:
            boolean: True if the points have distance less than eps in the infinity-norm, else False
    '''
    return (abs(x1 - x2) > eps) or (abs(y1-y2) > eps) or (abs(z1 - z2) > eps)

        
''' Code - do stuff '''

if (a_x >= b_x) or (a_y >= b_y) or (a_z >= b_z):
    raise ValueError("Your plotting range has too few dimensions or the boundarys dont make sense (e.g. intervall [2,1]). If you want to plot something 2-dim there is an extra script for that.")

# Get the relevant entries of the matrix A. This saves us a lot of asking for entries not beeing infty
inds = []
for l in range((d+1)**3):
    i,j,k = linToExp(l)
    if not (A[i,j,k] == infty):
        inds.append(l)

# The number of non infty entries in A
n = len(inds)

if (n < 2):
    raise ValueError("Please enter a polynomial, that is, set at least two entries of A to something different than infty.")

# A list of the corners of the plotting range. We will need this later
corners = [[a_x, a_y, a_z],
           [a_x, a_y, b_z],
           [a_x, b_y, a_z],
           [a_x, b_y, b_z],
           [b_x, a_y, a_z],
           [b_x, a_y, b_z],
           [b_x, b_y, a_z],
           [b_x, b_y, b_z]]

            
''' Walk through matrix A, check all quadrupels of monomials '''

list_x_vals = [] # In here go the x values for quadruple-max points
list_y_vals = [] # In here go the y values for quadruple-max points
list_z_vals = []
list_is = [] # In here go, as lists, the x-exponents of the maximal monomials
list_js = [] # In here go, as lists, the y-exponents of the maximal monomials
list_ks = []


for i in range(n-3):
    for j in range(i+1,n-2):
        for k in range(j+1,n-1):
            for l in range(k+1,n):
                ## Solve for points where a value is attained three times
                # Recover the exponents from the linear counting
                i1, j1, k1 = linToExp(inds[i])
                i2, j2, k2 = linToExp(inds[j])
                i3, j3, k3 = linToExp(inds[k])
                i4, j4, k4 = linToExp(inds[l])

                # The coefficients of the system of linear equations we have to solve
                B = np.array([[i2-i1, j2-j1, k2-k1],[i3-i1, j3-j1, k3-k1],[i4-i1, j4-j1, k4-k1]])
                b = np.array([A[i1,j1,k1] - A[i2,j2,k2],A[i1,j1,k1] - A[i3,j3,k3],A[i1,j1,k1] - A[i4,j4,k4]])

                # Proceed only if the system of lin eq is uniquely solvable
                if (0 != np.linalg.det(B)):
                    xyz = np.linalg.solve(B, b)
                    # Check if it is actually a maximum
                    if (abs(np.max(tropPolySet(xyz[0],xyz[1],xyz[2],A,inds)) - A[i1,j1,k1] - i1*xyz[0] - j1*xyz[1] - k1*xyz[2]) < eps):
                        list_x_vals.append(xyz[0])
                        list_y_vals.append(xyz[1])
                        list_z_vals.append(xyz[2])
                        if (xyz[0] <= a_x) or (xyz[0] >= b_x) or (xyz[1] <= a_y) or (xyz[1] >= b_y) or (xyz[2] <= a_z) or (xyz[2] >= b_z):
                            print(f"The point ({xyz[0]},{xyz[1]},{xyz[2]}) is not in your plotting range, which should be a problem...")

# Get all max monomials at each point
for i in range(len(list_x_vals)):
    temp = maxIndices(tropPolySet(list_x_vals[i],list_y_vals[i],list_z_vals[i],A,inds),eps)
    temp_is = []
    temp_js = []
    temp_ks = []
    for tp in temp:
        i1,j1,k1 = linToExp(inds[tp])
        temp_is.append(i1)
        temp_js.append(j1)
        temp_ks.append(k1)
    list_is.append(temp_is)
    list_js.append(temp_js)
    list_ks.append(temp_ks)

''' Now consider the boundary-triple max points '''

list_x_vals_bdry = [] # In here go the x values for points in the trop variety on the boundary
list_y_vals_bdry = [] # In here go the y values for points in the trop variety on the boundary
list_z_vals_bdry = []
list_is_bdry = [] # In here go, as length 3 lists, the x-exponents of the maximal monomials
list_js_bdry = [] # In here go, as length 3 listss, the y-exponents of the maximal monomials
list_ks_bdry = []

# It follows a lot of almost copy-pasted code. However, as there are always small but significiant changes, it cannot be easily be turned into a function
for i in range(n-2):
    for j in range(i+1,n-1):
        for k in range(j+1,n):
            ## Solve for points on boundary where a value is attained three times
            # Recovers the exponents from the linear counting

            i1, j1, k1 = linToExp(inds[i])
            i2, j2, k2 = linToExp(inds[j])
            i3, j3, k3 = linToExp(inds[k])
            
            # a_x boundary
            B = np.array([[j2-j1, k2-k1],[j3-j1, k3-k1]])
            b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(i2-i1)*a_x,A[i1,j1,k1]-A[i3,j3,k3]-(i3-i1)*a_x])
            if (np.linalg.det(B) != 0): # If uniquely solvable
                rest_coord = np.linalg.solve(B, b)
                # if in plotting range
                if (rest_coord[0] >= a_y) and (rest_coord[0] <= b_y) and (rest_coord[1] >= a_z) and (rest_coord[1] <= b_z):
                    # if actually in the variety
                    if(abs(np.max(tropPolySet(a_x,rest_coord[0],rest_coord[1],A,inds)) - A[i1,j1,k1] - i1*a_x - j1* rest_coord[0] - k1*rest_coord[1]) < eps):
                        list_x_vals_bdry.append(a_x)
                        list_y_vals_bdry.append(rest_coord[0])
                        list_z_vals_bdry.append(rest_coord[1])

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])
                # b_x
                b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(i2-i1)*b_x,A[i1,j1,k1]-A[i3,j3,k3]-(i3-i1)*b_x])
                rest_coord = np.linalg.solve(B, b)
                #if in plotting range
                if (rest_coord[0] >= a_y) and (rest_coord[0] <= b_y) and (rest_coord[1] >= a_z) and (rest_coord[1] <= b_z):
                    # if actually in the variety
                    if(abs(np.max(tropPolySet(b_x,rest_coord[0],rest_coord[1],A,inds)) - A[i1,j1,k1] - i1*b_x - j1* rest_coord[0] - k1*rest_coord[1]) < eps):
                        list_x_vals_bdry.append(b_x)
                        list_y_vals_bdry.append(rest_coord[0])
                        list_z_vals_bdry.append(rest_coord[1])

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])
            # a_y
            B = np.array([[i2-i1, k2-k1],[i3-i1, k3-k1]])
            b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(j2-j1)*a_y, A[i1,j1,k1]-A[i3,j3,k3]-(j3-j1)*a_y])
            if (np.linalg.det(B) != 0): # If uniquely solvable
                rest_coord = np.linalg.solve(B, b)
                # if in plotting range
                if (rest_coord[0] >= a_x) and (rest_coord[0] <= b_x) and (rest_coord[1] >= a_z) and (rest_coord[1] <= b_z):
                    # if in variety
                    if(abs(np.max(tropPolySet(rest_coord[0],a_y,rest_coord[1],A,inds)) - A[i1,j1,k1] - i1*rest_coord[0] - j1* a_y - k1*rest_coord[1]) < eps):
                        list_x_vals_bdry.append(rest_coord[0])
                        list_y_vals_bdry.append(a_y)
                        list_z_vals_bdry.append(rest_coord[1])

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])
                # b_y
                b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(j2-j1)*b_y, A[i1,j1,k1]-A[i3,j3,k3]-(j3-j1)*b_y])
                rest_coord = np.linalg.solve(B, b)
                # if in range
                if (rest_coord[0] >= a_x) and (rest_coord[0] <= b_x) and (rest_coord[1] >= a_z) and (rest_coord[1] <= b_z):
                    if(abs(np.max(tropPolySet(rest_coord[0],b_y,rest_coord[1],A,inds)) - A[i1,j1,k1] - i1*rest_coord[0] - j1* b_y - k1*rest_coord[1]) < eps):
                        list_x_vals_bdry.append(rest_coord[0])
                        list_y_vals_bdry.append(b_y)
                        list_z_vals_bdry.append(rest_coord[1])

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])
            # a_z
            B = np.array([[i2-i1, j2-j1],[i3-i1, j3-j1]])
            b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(k2-k1)*a_z, A[i1,j1,k1]-A[i3,j3,k3]-(k3-k1)*a_z])
            if (np.linalg.det(B) != 0): # If uniquely solvable
                rest_coord = np.linalg.solve(B, b)
                #if in range
                if (rest_coord[0] >= a_x) and (rest_coord[0] <= b_x) and (rest_coord[1] >= a_y) and (rest_coord[1] <= b_y):
                    #if in variety
                    if(abs(np.max(tropPolySet(rest_coord[0],rest_coord[1],a_z,A,inds)) - A[i1,j1,k1] - i1*rest_coord[0] - j1* rest_coord[1] - k1*a_z) < eps):
                        list_x_vals_bdry.append(rest_coord[0])
                        list_y_vals_bdry.append(rest_coord[1])
                        list_z_vals_bdry.append(a_z)

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])
                # b_z
                b = np.array([A[i1,j1,k1]-A[i2,j2,k2]-(k2-k1)*b_z, A[i1,j1,k1]-A[i3,j3,k3]-(k3-k1)*b_z])
                rest_coord = np.linalg.solve(B, b)
                if (rest_coord[0] >= a_x) and (rest_coord[0] <= b_x) and (rest_coord[1] >= a_y) and (rest_coord[1] <= b_y):
                    if(abs(np.max(tropPolySet(rest_coord[0],rest_coord[1],b_z,A,inds)) - A[i1,j1,k1] - i1*rest_coord[0] - j1* rest_coord[1] - k1*b_z) < eps):
                        list_x_vals_bdry.append(rest_coord[0])
                        list_y_vals_bdry.append(rest_coord[1])
                        list_z_vals_bdry.append(b_z)

                        list_is_bdry.append([i1,i2,i3])
                        list_js_bdry.append([j1,j2,j3])
                        list_ks_bdry.append([k1,k2,k3])


'''Now we have all points at wich edges start and end and the exponents of the maximal monomials'''
# We check which points we should connect with a line
# Sadly just checking if the point in the middle is in the variety is not sufficient.
# The max has to be attained by the right monomials, too

m = len(list_x_vals) # We will need this number often. So let's just remember it
list_x_edge_vals = [] # In here go the x values, as length 2 arrays, of start and endpoint of edges of the trop variety
list_y_edge_vals = [] # In here go the y values, as length 2 arrays, of start and endpoint of edges of the trop variety
list_z_edge_vals = [] # In here go the y values, as length 2 arrays, of start and endpoint of edges of the trop variety


''' First: Check which 4-max points to connect to other 4-max points'''

for i in range(m-1):
    for j in range(i+1,m):
        # If the two points are actually different points
        if pointsUnEqual(list_x_vals[i],list_y_vals[i],list_z_vals[i],list_x_vals[j],list_y_vals[j],list_z_vals[j],eps):
            mid_x = (list_x_vals[i] + list_x_vals[j])/2
            mid_y = (list_y_vals[i] + list_y_vals[j])/2
            mid_z = (list_z_vals[i] + list_z_vals[j])/2
            val = np.max(tropPolySet(mid_x,mid_y,mid_z,A,inds))
            a = []
            b = []
            for k in range(0,len(list_is[i])): # Fill a with the values of the max monomials of point one at the middle point
                a.append(A[list_is[i][k],list_js[i][k],list_ks[i][k]] + list_is[i][k]*mid_x + list_js[i][k]*mid_y+list_ks[i][k]*mid_z)
            for k in range(0,len(list_is[j])): # Fill b with the values of the max monomials of point two at the middle point
                b.append(A[list_is[j][k],list_js[j][k],list_ks[j][k]] + list_is[j][k]*mid_x + list_js[j][k]*mid_y + +list_ks[j][k]*mid_z)
            if threeEqual(a,b,val,eps): # in trop var the max is attained by the right monomials at the middle point?
                list_x_edge_vals.append(np.array([list_x_vals[i],list_x_vals[j]]))
                list_y_edge_vals.append(np.array([list_y_vals[i],list_y_vals[j]]))
                list_z_edge_vals.append(np.array([list_z_vals[i],list_z_vals[j]]))


''' Second: Check which 4-max points to connect to triple max boundary points'''

for i in range(m):
    for j in range(len(list_x_vals_bdry)):
        # If the two points are actually different points
        if pointsUnEqual(list_x_vals[i],list_y_vals[i],list_z_vals[i],list_x_vals_bdry[j],list_y_vals_bdry[j],list_z_vals_bdry[j],eps):
            # If the max is attained by the right monomials
            if isSubSet3(list_is_bdry[j],list_js_bdry[j],list_ks_bdry[j],list_is[i],list_js[i],list_ks[i]):
                list_x_edge_vals.append(np.array([list_x_vals[i],list_x_vals_bdry[j]]))
                list_y_edge_vals.append(np.array([list_y_vals[i],list_y_vals_bdry[j]]))
                list_z_edge_vals.append(np.array([list_z_vals[i],list_z_vals_bdry[j]]))


''' Third: Check which boundary triple max points to connect.'''
# This is a degenerate case, but it can happen - but only if there are no 4-max points

if not list_x_vals:
    for i in range(len(list_x_vals_bdry)):
        for j in range(len(list_x_vals_bdry)):
            # If the two points are actually different points
            if pointsUnEqual(list_x_vals_bdry[i],list_y_vals_bdry[i],list_z_vals_bdry[i],list_x_vals_bdry[j],list_y_vals_bdry[j],list_z_vals_bdry[j],eps):
                if isSubSet3(list_is_bdry[j],list_js_bdry[j],list_ks_bdry[j], list_is_bdry[i],list_js_bdry[i],list_ks_bdry[i]):
                    list_x_edge_vals.append(np.array([list_x_vals_bdry[i],list_x_vals_bdry[j]]))
                    list_y_edge_vals.append(np.array([list_y_vals_bdry[i],list_y_vals_bdry[j]]))
                    list_z_edge_vals.append(np.array([list_z_vals_bdry[i],list_z_vals_bdry[j]]))


'''Compute the intersections with the edges of the plotting range'''

list_x_vals_bbdry = [] # List for the x values of said points
list_y_vals_bbdry = []
list_z_vals_bbdry = []

list_is_bbdry = [] # List for the exponents
list_js_bbdry = []
list_ks_bbdry = []

if list_is_bdry:
    # If this list is non empty, it automatically has at least length two
    list_exp_bdry = [] # List for the max exponents at boundary points. These monoms are the only ones which can be maximal at the edges of the range
    for i in range(len(list_is_bdry)):
        temp = []
        temp.append([list_is_bdry[i][0], list_js_bdry[i][0], list_ks_bdry[i][0]])
        temp.append([list_is_bdry[i][1], list_js_bdry[i][1], list_ks_bdry[i][1]])
        temp.append([list_is_bdry[i][2], list_js_bdry[i][2], list_ks_bdry[i][2]])
        list_exp_bdry.append(temp)
    # Go through pairs of boundary triple max points
    for exp1 in range(len(list_exp_bdry)):
        for exp2 in range(exp1 + 1,len(list_exp_bdry)):
            # Get the monomials that agree
            monoms = getTheEquals(list_exp_bdry[exp1],list_exp_bdry[exp2])
            # If more at least two monomials agree
            if (len(monoms) > 1):
                # Go through all pairs of elements in monoms
                for mon1 in range(len(monoms)):
                    for mon2 in range(mon1 + 1, len(monoms)):
                        i1 = monoms[mon1][0]
                        j1 = monoms[mon1][1]
                        k1 = monoms[mon1][2]

                        i2 = monoms[mon2][0]
                        j2 = monoms[mon2][1]
                        k2 = monoms[mon2][2]
                        
                        di = i1 - i2
                        dj = j1 - j2
                        dk = k1 - k2
                        # Solve for the intersection points with the edges of the range
                        # This is, again, a lot of copy-paste work with small variations
                        if (dk != 0):# ax,ay - ax,by - bx,ay - bx,by
                            # ax,ay
                            last_val = (-1/dk)*(A[i1,j1,k1]-A[i2,j2,k2]+a_x*di+a_y*dj)
                            if (abs(np.max(tropPolySet(a_x,a_y,last_val,A,inds) - tropMonomVal(a_x,a_y,last_val,monoms[mon1],A)))<eps) and (last_val >= a_z) and (last_val <= b_z):
                                list_x_vals_bbdry.append(a_x)
                                list_y_vals_bbdry.append(a_y)
                                list_z_vals_bbdry.append(last_val)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            # ax,by
                            last_val = (-1/dk)*(A[i1,j1,k1]-A[i2,j2,k2]+a_x*di+b_y*dj)
                            if (abs(np.max(tropPolySet(a_x,b_y,last_val,A,inds) - tropMonomVal(a_x,b_y,last_val,monoms[mon1],A)))<eps) and (last_val >= a_z) and (last_val <= b_z):
                                list_x_vals_bbdry.append(a_x)
                                list_y_vals_bbdry.append(b_y)
                                list_z_vals_bbdry.append(last_val)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            # bx,ay
                            last_val = (-1/dk)*(A[i1,j1,k1]-A[i2,j2,k2]+b_x*di+a_y*dj)
                            if (abs(np.max(tropPolySet(b_x,a_y,last_val,A,inds) - tropMonomVal(b_x,a_y,last_val,monoms[mon1],A)))<eps) and (last_val >= a_z) and (last_val <= b_z):
                                list_x_vals_bbdry.append(b_x)
                                list_y_vals_bbdry.append(a_y)
                                list_z_vals_bbdry.append(last_val)
                                
                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            # bx,by
                            last_val = (-1/dk)*(A[i1,j1,k1]-A[i2,j2,k2]+b_x*di+b_y*dj)
                            if (abs(np.max(tropPolySet(b_x,b_y,last_val,A,inds) - tropMonomVal(b_x,b_y,last_val,monoms[mon1],A)))<eps) and (last_val >= a_z) and (last_val <= b_z):
                                list_x_vals_bbdry.append(b_x)
                                list_y_vals_bbdry.append(b_y)
                                list_z_vals_bbdry.append(last_val)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                        if (dj != 0):# ax,az - ax,bz - bx,az - bx,bz
                            #ax,az
                            last_val = (-1/dj)*(A[i1,j1,k1]-A[i2,j2,k2]+a_x*di+a_z*dk)
                            if (abs(np.max(tropPolySet(a_x,last_val,a_z,A,inds) - tropMonomVal(a_x,last_val,a_z,monoms[mon1],A)))<eps) and (last_val >= a_y) and (last_val <= b_y):
                                list_x_vals_bbdry.append(a_x)
                                list_y_vals_bbdry.append(last_val)
                                list_z_vals_bbdry.append(a_z)
                                
                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #ax,bz
                            last_val = (-1/dj)*(A[i1,j1,k1]-A[i2,j2,k2]+a_x*di+b_z*dk)
                            if (abs(np.max(tropPolySet(a_x,last_val,b_z,A,inds) - tropMonomVal(a_x,last_val,b_z,monoms[mon1],A)))<eps) and (last_val >= a_y) and (last_val <= b_y):
                                list_x_vals_bbdry.append(a_x)
                                list_y_vals_bbdry.append(last_val)
                                list_z_vals_bbdry.append(b_z)
                                
                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #bx,az
                            last_val = (-1/dj)*(A[i1,j1,k1]-A[i2,j2,k2]+b_x*di+a_z*dk)
                            if (abs(np.max(tropPolySet(b_x,last_val,a_z,A,inds) - tropMonomVal(b_x,last_val,a_z,monoms[mon1],A)))<eps) and (last_val >= a_y) and (last_val <= b_y):
                                list_x_vals_bbdry.append(b_x)
                                list_y_vals_bbdry.append(last_val)
                                list_z_vals_bbdry.append(a_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #bx,bz
                            last_val = (-1/dj)*(A[i1,j1,k1]-A[i2,j2,k2]+b_x*di+b_z*dk)
                            if (abs(np.max(tropPolySet(b_x,last_val,b_z,A,inds) - tropMonomVal(b_x,last_val,b_z,monoms[mon1],A)))<eps) and (last_val >= a_y) and (last_val <= b_y):
                                list_x_vals_bbdry.append(b_x)
                                list_y_vals_bbdry.append(last_val)
                                list_z_vals_bbdry.append(b_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                        if (di != 0):# ayaz, aybz, byaz, bybz
                            #ayaz
                            last_val = (-1/di)*(A[i1,j1,k1]-A[i2,j2,k2]+a_y*dj+a_z*dk)
                            if (abs(np.max(tropPolySet(last_val,a_y,a_z,A,inds) - tropMonomVal(last_val,a_y,a_z,monoms[mon1],A)))<eps) and (last_val >= a_x) and (last_val <= b_x):
                                list_x_vals_bbdry.append(last_val)
                                list_y_vals_bbdry.append(a_y)
                                list_z_vals_bbdry.append(a_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #aybz
                            last_val = (-1/di)*(A[i1,j1,k1]-A[i2,j2,k2]+a_y*dj+b_z*dk)
                            if (abs(np.max(tropPolySet(last_val,a_y,b_z,A,inds) - tropMonomVal(last_val,a_y,b_z,monoms[mon1],A)))<eps) and (last_val >= a_x) and (last_val <= b_x):
                                list_x_vals_bbdry.append(last_val)
                                list_y_vals_bbdry.append(a_y)
                                list_z_vals_bbdry.append(b_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #byaz
                            last_val = (-1/di)*(A[i1,j1,k1]-A[i2,j2,k2]+b_y*dj+a_z*dk)
                            if (abs(np.max(tropPolySet(last_val,b_y,a_z,A,inds) - tropMonomVal(last_val,b_y,a_z,monoms[mon1],A)))<eps) and (last_val >= a_x) and (last_val <= b_x):
                                list_x_vals_bbdry.append(last_val)
                                list_y_vals_bbdry.append(b_y)
                                list_z_vals_bbdry.append(a_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])
                            #bybz
                            last_val = (-1/di)*(A[i1,j1,k1]-A[i2,j2,k2]+b_y*dj+b_z*dk)
                            if (abs(np.max(tropPolySet(last_val,b_y,b_z,A,inds) - tropMonomVal(last_val,b_y,b_z,monoms[mon1],A)))<eps) and (last_val >= a_x) and (last_val <= b_x):
                                list_x_vals_bbdry.append(last_val)
                                list_y_vals_bbdry.append(b_y)
                                list_z_vals_bbdry.append(b_z)

                                list_is_bbdry.append([i1,i2])
                                list_js_bbdry.append([j1,j2])
                                list_ks_bbdry.append([k1,k2])

# Initialisation for plotting the dual
x_vals_edge_dual = []
y_vals_edge_dual = []
z_vals_edge_dual = []

# Initialisation for the plots
fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12,6))
ax1.set_xlim([a_x, b_x])
ax1.set_ylim([a_y, b_y])
ax1.set_zlim([a_z, b_z])

# If there are edges. In the else case, the variety is just a plane
if list_x_edge_vals:
    '''Compute the intersections with the corners of the plotting range'''
    # There shouldn't be too many, but it can happen and f the plot up if we don't check for them
    list_x_vals_bbbdry = []
    list_y_vals_bbbdry = []
    list_z_vals_bbbdry = []

    list_is_bbbdry = []
    list_js_bbbdry = []
    list_ks_bbbdry = []

    for corner in corners:
        values_poly = tropPolySet(corner[0],corner[1],corner[2],A,inds)
        corner_monom = maxMonomsAt(values_poly,eps)
        if (len(corner_monom) > 1):
            list_x_vals_bbbdry.append(corner[0])
            list_y_vals_bbbdry.append(corner[1])
            list_z_vals_bbbdry.append(corner[2])
            temp_i = []
            temp_j = []
            temp_k = []
            for exp in corner_monom:
                temp_i.append(exp[0])
                temp_j.append(exp[1])
                temp_k.append(exp[2])
            list_is_bbbdry.append(temp_i)
            list_js_bbbdry.append(temp_j)
            list_ks_bbbdry.append(temp_k)


    '''Go for the 2-max planes'''
    # Go through all possible pairs of monomials and check which vertices contain them in their max monomials
    
    vertices_polygon = [] # In here go the coordinates of a polygon as a list

    # Get all the possible coordinates of the polygons of the variety
    list_x_vals_all = list_x_vals + list_x_vals_bdry + list_x_vals_bbdry + list_x_vals_bbbdry
    list_y_vals_all = list_y_vals + list_y_vals_bdry + list_y_vals_bbdry + list_y_vals_bbbdry
    list_z_vals_all = list_z_vals + list_z_vals_bdry + list_z_vals_bbdry + list_z_vals_bbbdry
    # And the exponents of the maximal monomials
    list_is_all = list_is + list_is_bdry + list_is_bbdry + list_is_bbbdry
    list_js_all = list_js + list_js_bdry + list_js_bbdry + list_js_bbbdry
    list_ks_all = list_ks + list_ks_bdry + list_ks_bbdry + list_ks_bbbdry

    num_polygon = 0 # Counts the polygons

    ## Compute the vertices of all polygons
    # Polygons with more than two maximal monomials occur more often than once

    max_monoms = [] # At entry i there is a list of the exponents of the maximal monomials that lead to the i th polygon

    # Go through pairs of monomials
    for i in range(n):
        for j in range(i+1,n):
            temp = [] # In here go the relevant indices
            
            # Gives us the exponents of the monomials
            i1, j1, k1 = linToExp(inds[i])
            i2, j2, k2 = linToExp(inds[j])

            for verts in range(len(list_x_vals_all)):
                # If max is attained by the two monomials
                if bothAreIn(i1,j1,k1,i2,j2,k2,list_is_all[verts],list_js_all[verts],list_ks_all[verts]):
                    temp.append(verts)
            # If it is actually a polygon
            if (len(temp) > 2):
                max_monoms.append([[i1,j1,k1],[i2,j2,k2]])
                vertices_polygon.append([])
                for k in range(len(temp)):
                    vertices_polygon[num_polygon].append(np.array([list_x_vals_all[temp[k]], list_y_vals_all[temp[k]], list_z_vals_all[temp[k]]]))
                num_polygon = num_polygon + 1

    # For each polygon, each vertex should be mentioned only once
    for i in range(len(vertices_polygon)):
        vertices_polygon[i]= noDouble(vertices_polygon[i],eps)

    # And each polygon should be mentionned only once
    # Memorize the maximal monomials at each polygon
    vertices_polygon, max_monoms  = noDoublePolygon(vertices_polygon, max_monoms)

    ''' Plotting'''

    # Draw the lines
    for i in range(len(list_x_edge_vals)):
        ax1.plot(list_x_edge_vals[i],list_y_edge_vals[i],list_z_edge_vals[i],color='black')
        
    # Draw the polygons & compute the dual
    for i in range(len(vertices_polygon)):
        orth = np.cross(vertices_polygon[i][0]-vertices_polygon[i][1],vertices_polygon[i][0]-vertices_polygon[i][2])
        # If plane not orthogonal to x_3 = 0: Take first two coords for convex hull
        if (abs(np.dot(orth,np.array([0,0,1]))) > eps):
            proj_2d = [array[:2] for array in vertices_polygon[i]]
            hull_2d = ConvexHull(proj_2d)
            verts = []
            for ind in hull_2d.vertices:
                verts.append(vertices_polygon[i][ind])
            poly = Poly3DCollection([verts], facecolors='b', linewidths=0, edgecolors='r', alpha=0.5)
            ax1.add_collection3d(poly)

            # (0,0,1) is not in direction of polygon
            max_monoms[i] = sorted(max_monoms[i], key=lambda x: x[2])
            temp_len = len(max_monoms[i])-1
            x_vals_edge_dual.append([max_monoms[i][0][0],max_monoms[i][temp_len][0]])
            y_vals_edge_dual.append([max_monoms[i][0][1],max_monoms[i][temp_len][1]])
            z_vals_edge_dual.append([max_monoms[i][0][2],max_monoms[i][temp_len][2]])
        # else if not orthogonal to x_2 = 0: Take 1st and 3rd coord
        elif (abs(np.dot(orth,np.array([0,1,0]))) > eps):
            proj_2d = [array[[0,2]] for array in vertices_polygon[i]]
            hull_2d = ConvexHull(proj_2d)
            verts = []
            for ind in hull_2d.vertices:
                verts.append(vertices_polygon[i][ind])
            poly = Poly3DCollection([verts], facecolors='b', linewidths=0, edgecolors='r', alpha=0.5)
            ax1.add_collection3d(poly)

            # (0,1,0) is not in direction of polygon
            max_monoms[i] = sorted(max_monoms[i], key=lambda x: x[1])
            temp_len = len(max_monoms[i])-1
            x_vals_edge_dual.append([max_monoms[i][0][0],max_monoms[i][temp_len][0]])
            y_vals_edge_dual.append([max_monoms[i][0][1],max_monoms[i][temp_len][1]])
            z_vals_edge_dual.append([max_monoms[i][0][2],max_monoms[i][temp_len][2]])
        # else if not orthogonal to 1_1 = 0: Take 2nd and 3rd
        elif (abs(np.dot(orth,np.array([1,0,0]))) > eps):
            proj_2d = [array[[1,2]] for array in vertices_polygon[i]]
            hull_2d = ConvexHull(proj_2d)
            verts = []
            for ind in hull_2d.vertices:
                verts.append(vertices_polygon[i][ind])
            poly = Poly3DCollection([verts], facecolors='b', linewidths=0, edgecolors='r', alpha=0.5)
            ax1.add_collection3d(poly)
            
            # (1,0,0) is not in direction of polygon
            max_monoms[i] = sorted(max_monoms[i], key=lambda x: x[0])
            temp_len = len(max_monoms[i])-1
            x_vals_edge_dual.append([max_monoms[i][0][0],max_monoms[i][temp_len][0]])
            y_vals_edge_dual.append([max_monoms[i][0][1],max_monoms[i][temp_len][1]])
            z_vals_edge_dual.append([max_monoms[i][0][2],max_monoms[i][temp_len][2]])
        else:
            print('???')
    

else:
    # This is a very degenerate case: We will get a plane.
    # In this case, however, in the plot of the dual it does not tell you which monomials ever attain a maximum. Improve this? TODO
    
    # Check the values at the corners of the plotting range. If the variety has an honest intersection with the plotting range, this gives us all the information we need
    # We do not care about the values, only about the max monoms

    inds_rel = []
    count = 0
    # Get the two relevant monomials
    while (len(inds_rel) < 2):
        if count == 8:
            raise ValueError("The surface does not intersect your plotting range")
        inds_temp = maxIndices(tropPolySet(corners[count][0],corners[count][1],corners[count][2],A,inds),eps)
        if (len(inds_temp) == 1) and (not (inds_temp[0] in inds_rel)):
            inds_rel.append(inds_temp[0])
        count += 1
    exp1 = linToExp(inds_rel[0])
    exp2 = linToExp(inds_rel[1])

    x_vals_edge_dual.append([exp1[0],exp2[0]])
    y_vals_edge_dual.append([exp1[1],exp2[1]])
    z_vals_edge_dual.append([exp1[2],exp2[2]])
    
    # Equation of the plane is given by:
    # A[exp1[0],exp1[1],exp1[2]] - A[exp2[0],exp2[1],exp2[2]] + (exp1[0] - exp2[0]) * x + (exp1[1] - exp2[1]) * y + (exp1[2] - exp2[2]) * z = 0
    if (abs(exp1[2] - exp2[2]) > eps):
        # We can solve for z
        X = np.linspace(a_x, b_x, 100)
        Y = np.linspace(a_y, b_y, 100)
        X, Y = np.meshgrid(X, Y)

        Z = (A[exp1[0],exp1[1],exp1[2]] - A[exp2[0],exp2[1],exp2[2]] + (exp1[0] - exp2[0]) * X + (exp1[1] - exp2[1]) * Y) / (exp2[2] - exp1[2])
        Z = np.where((Z >= a_z) & (Z <= b_z), Z, np.nan)
        if np.all(np.isnan(Z)):
            print("The surface does not intersect your plotting range")
        # plot
        ax1.plot_surface(X, Y, Z, alpha=0.5, rstride=100, cstride=100, cmap='viridis', facecolors = 'b')
    elif (abs(exp1[1] - exp2[1]) > eps):
        # We can solve for y and the z part is zero
        X = np.linspace(a_x, b_x, 100)
        Z = np.linspace(a_z, b_z, 100)
        X, Z = np.meshgrid(X, Z)

        Y = (A[exp1[0],exp1[1],exp1[2]] - A[exp2[0],exp2[1],exp2[2]] + (exp1[0] - exp2[0]) * X) / (exp2[1] - exp1[1])
        Y = np.where((Y >= a_y) & (Y <= b_y), Y, np.nan)
        if np.all(np.isnan(Y)):
            print("The surface does not intersect your plotting range")
        #plot
        ax1.plot_surface(X, Y, Z, alpha=0.5, rstride=100, cstride=100, cmap='viridis', facecolors = 'b')
    elif (abs(exp1[0] - exp2[0]) > eps):
        # We can solve for x and y,z part are zero
        Y = np.linspace(a_y, b_y, 100)
        Z = np.linspace(a_z, b_z, 100)
        Y, Z = np.meshgrid(Y, Z)

        X = (A[exp1[0],exp1[1],exp1[2]] - A[exp2[0],exp2[1],exp2[2]]) / (exp2[0] - exp1[0])
        X = np.where((X >= a_x) & (X <= b_x), X, np.nan)
        if np.all(np.isnan(X)):
            print("The surface does not intersect your plotting range")
        #plot
        ax1.plot_surface(X, Y, Z, alpha=0.5, rstride=100, cstride=100, cmap='viridis', facecolors = 'b')
    else:
        print("This shouldn't happen... in the sense that it should be impossible. Case with only one plane but there was a contradiction in the equations?!")


''' Draw the dual '''
# Create the coordinates for \Delta_deg
x_coords_dual = []
y_coords_dual = []
z_coords_dual = []

# Compute the degree of the polynomial
degF = 0
for i in range(n):
    i1, j1, k1 = linToExp(inds[i])
    if (i1 + j1 + k1 > degF):
        degF = i1 + j1 + k1
        
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
