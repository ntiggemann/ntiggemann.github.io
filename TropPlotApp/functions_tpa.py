import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.ticker import MaxNLocator
import itertools
import math

''' Functions '''

def evaluate_trop_monomials(x,y,A):
    '''
    Given a point (x,y) and a tropical polynomial with coefficient dictionary A, returns the values of the tropical 
    monomials.
    Args:
        x (int or float): x-coordinate
        y (int or float): y-coordinate
        A (dictionars): The coefficients of a tropical polynomial with key the exponents
    Returns:
        array of floats: The values in the max of the tropical polynomial
    '''
    values = np.zeros(len(A))
    count = 0

    for (i1,j1), coeff in A.items():
        values[count] = (i1)*x + (j1)*y + coeff
        count = count + 1
    return values

def get_values_indices(a,b,value,eps):
    '''
    Given lists of values a and b, returns lists of indices of a and b for which value is attained up to eps
    Args:
        a(list or array of int or float): First tuple
        b(list or array of int or float): Second tuple
        val (int or float): A value
        eps (int or float): error tolerance
    Returns:
        aind (list of int): list of indices for which a has value val up to eps
        bind (list of int): list of indices for which b has value val up to eps
    '''    
    a_value_indices = []
    b_value_indices = []

    for i in range(len(a)):
        if (abs(a[i] - value) < eps):
            a_value_indices.append(i)
    for i in range(len(b)):
        if (abs(b[i] - value) < eps):
            b_value_indices.append(i)

    return a_value_indices, b_value_indices
            
def write_the_polynomial(A,eps,is_title):
    '''
    Given the coefficients of a tropical polynomial, generates a string of LaTeX code for it
    Args:
        A (dictionary): Coefficients of a polynomial with key the exponents
        eps (float): Difference before two numbers are considered equal
        is_title (boolean): Decides if there will be breaks in the string 
    Returns
        string: LaTeX code that writes the polynomial
    '''

    if not is_title:
        first_monom = True
        for (i1,j1), coeff in A.items():
            if first_monom:
                first_monom = False
                if (abs(coeff-round(coeff,0)) < eps):
                    str_coeff = f"{int(round(coeff,0))}"
                else:
                    str_coeff = f"{coeff}"
                latex_polynomial = f"{str_coeff}\\odot x^{ {i1} }y^{ {j1} }" #f"{coeff}\\odot x^{{({i},{j})}}"
            else:
                if (abs(coeff-round(coeff,0)) < eps):
                    str_coeff = f"{int(round(coeff,0))}"
                else:
                    str_coeff = f"{coeff}"
                latex_polynomial = latex_polynomial + f"\\oplus {str_coeff}\\odot x^{ {i1} }y^{ {j1} }"    
    
    else:
        monom_count = 0
        for (i1,j1), coeff in A.items():
            if monom_count == 0:
                if (abs(coeff-round(coeff,0)) < eps):
                    str_coeff = f"{int(round(coeff,0))}"
                else:
                    str_coeff = f"{coeff}"
                latex_polynomial = f"${str_coeff}\\odot x^{ {i1} }y^{ {j1} }$" #f"{coeff}\\odot x^{{({i},{j})}}"
                monom_count = monom_count +1
            else:
                if (abs(coeff-round(coeff,0)) < eps):
                    str_coeff = f"{int(round(coeff,0))}"
                else:
                    str_coeff = f"{coeff}"
                if monom_count != 0 and monom_count % 5 == 0:
                    latex_polynomial = latex_polynomial + f"$\\oplus {str_coeff}\\odot x^{ {i1} }y^{ {j1} }$" + "\n"
                else:
                    latex_polynomial = latex_polynomial + f"$\\oplus {str_coeff}\\odot x^{ {i1} }y^{ {j1} }$"
                monom_count = monom_count + 1
    return latex_polynomial

def change_signs_of_polynomial(A):
    '''
    Args:
        A (dictionary): The coefficients
    Returns:
        Dictionary A with the signs of all entries indicated by non_infty_entries changed.
    '''
    for exp, coef in A.items():
        A[exp] = (-1)*coef
    return A

def tupel_is_in_list(x,L,eps):
    '''
    Args:
        x (list of length two of numbers): 
        L (list of lists of length two of numbers)
        eps : A (small) number
    Returns:
        i if x is the i-th entry of L up to eps in infty norm, -1 if x not in L
    '''
    for i in range(len(L)):
        if (abs(L[i][0] - x[0]) < eps) and (abs(L[i][1] - x[1]) < eps):
            return i
    
    return -1

def parse_polynomial_trop(poly_str):
    """
    Parses a polynomial string in two variables (x and y) and creates a dictionary.

    Args:
        poly_str (str): The polynomial string, e.g., "3.2x2y + yx + 9".

    Returns:
        dict: A dictionary representing the coefficients of the polynomial.
        int: The size of the array
    """
    # Remove spaces and handle '+' and '-' for separation
    poly_str = poly_str.replace(' ', '')
    poly_str = poly_str.replace('*', '')
    poly_str = poly_str.replace('^', '')
    poly_str = poly_str.replace('(', '')
    poly_str = poly_str.replace(')', '')
    poly_str = poly_str.replace('+-', '-')

    terms = re.findall(r'[+-]?\d*\.?\d*[xy]?\d*[xy]?\d*', poly_str)

    
    # Initialize a dictionary to store coefficients
    coefficients = {}

    for term in terms:
        if not term:  # Skip empty terms
            continue

        # Extract coefficient
        match = re.match(r'([+-]?\d*\.?\d*)([xy]\d*)?([xy]\d*)?', term)
        if not match:
            raise ValueError(f"Invalid term: {term}")

        coeff_str, var1, var2 = match.groups()

        # Coefficient handling
        coefficient = float(coeff_str) if coeff_str not in ('', '+', '-') else float(coeff_str + '0')

        # Initialize exponents
        x_exp, y_exp = 0, 0

        # Exponent handling for each variable
        for var in (var1, var2):
            if var:
                if var[0] == 'x':
                    x_exp = int(var[1:]) if len(var) > 1 else 1
                elif var[0] == 'y':
                    y_exp = int(var[1:]) if len(var) > 1 else 1

        # Store coefficient (add if multiple terms affect the same coefficient)
        if (x_exp, y_exp) in coefficients:
            coefficients[(x_exp, y_exp)] = max(coefficient,coefficients[(x_exp, y_exp)])
        else:
            coefficients[(x_exp, y_exp)] = coefficient

    # Neues Dictionary mit sortierten Schlüsseln
    sorted_coefficients = dict(sorted(coefficients.items(), key=lambda x: x[0], reverse=True))

    return sorted_coefficients

# Function to plot the zero set
def plot_trop_polynomial(poly_str,a_x,b_x,a_y,b_y,weight_text_size,weight_color,
                        color_curve,color_DS,color_NP,
                        eps):
    try:
        A = parse_polynomial_trop(poly_str)

        # The number of monomials
        number_of_monomials = len(A)

        if (number_of_monomials < 2):
            raise ValueError("Please enter a polynomial with at least two different monomials.")


        # Get ready for linear algebra!
        '''Compute the vertices of the curve, aka triple-max points'''

        x_vals_triple = [] # List of the x values for triple-max points
        y_vals_triple = [] # List of the y values for triple-max points
        is_triple = [] # List of length 3 arrays, the x-exponents of the maximal monomials
        js_triple = [] # List of length 3 arrays, the y-exponents of the maximal monomials

        missing_points = []

        # Walk through dictionary A, check all triples of monomials
        for ((i1,j1),(i2,j2),(i3,j3)) in itertools.combinations(A,3):
            ## Solve for points where a value is attained three times

            # The coefficients of the system of linear equations we have to solve
            B = np.array([[i1-i2, j1-j2],[i1-i3, j1-j3]])
            b = np.array([A[(i2,j2)] - A[(i1,j1)],A[(i3,j3)] - A[(i1,j1)]])
            # Proceed only if  the syst of lin eq is uniquely solvable
            if (0 != np.linalg.det(B)):
                x = np.linalg.solve(B, b)
                # Check if it is actually a maximum
                if (abs(np.max(evaluate_trop_monomials(x[0],x[1],A)) - A[(i1,j1)] - i1*x[0] - j1*x[1]) < eps):
                    x_vals_triple.append(x[0])
                    y_vals_triple.append(x[1])
                    is_triple.append(np.array([i1,i2,i3]))
                    js_triple.append(np.array([j1,j2,j3]))
                    if (x[0] <= a_x) or (x[0] >= b_x) or (x[1] <= a_y) or (x[1] >= b_y):
                        # If a vertex of the curve is not in the plot, it will look weird, and you will probably be unhappy with the plot
                        missing_points.append([x[0],x[1]])
        if missing_points:
            error_message = "The points\n"
            for vert in missing_points:
                error_message = error_message + f"({vert[0]},{vert[1]})\n"
            error_message = error_message + "are vertices of the curve and not in the plotting range, which will most likely ruin the plot"
#            messagebox.showerror("Missing vertices", error_message)


        ''' Compute the intersection points of the curve with the boundary of the plot'''

        x_vals_bdry = [] # List of the x values for points in the trop variety on the boundary
        y_vals_bdry = [] # List of the y values for points in the trop variety on the boundary
        is_bdry = [] # List of length 2 arrays, the x-exponents of the maximal monomials
        js_bdry = [] # List of length 2 arrays, the y-exponents of the maximal monomials

        # Check all pairs of monomials
        for ((i1,j1),(i2,j2)) in itertools.combinations(A,2):
            ## Solve for points on boundary where a value is attained twice

            rhs_a_x = A[(i2,j2)] + (i2 - i1)* a_x - A[(i1,j1)]
            rhs_b_x = A[(i2,j2)] + (i2 - i1)* b_x - A[(i1,j1)]
            rhs_a_y = A[(i2,j2)] + (j2 - j1)* a_y - A[(i1,j1)]
            rhs_b_y = A[(i2,j2)] + (j2 - j1)* b_y - A[(i1,j1)]

            if (j1-j2 != 0): # if equation solvable
                # if max is actually attained
                if (abs(np.max(evaluate_trop_monomials(a_x,(1/(j1-j2))*rhs_a_x,A))- A[(i1,j1)] - i1*a_x - j1*(1/(j1-j2))*rhs_a_x) < eps):
                    val_bdry = (1/(j1-j2))*rhs_a_x
                    # if point in plotting range
                    if (val_bdry >= a_y) and (val_bdry <= b_y):
                        x_vals_bdry.append(a_x)
                        y_vals_bdry.append(val_bdry)
                        is_bdry.append(np.array([i1,i2]))
                        js_bdry.append(np.array([j1,j2]))
                # if max is actually attained
                if (abs(np.max(evaluate_trop_monomials(b_x,(1/(j1-j2))*rhs_b_x,A)) - A[(i1,j1)] - i1*b_x - j1*(1/(j1-j2))*rhs_b_x) < eps):
                    val_bdry = (1/(j1-j2))*rhs_b_x
                    # if point in plotting range
                    if (val_bdry >= a_y) and (val_bdry <= b_y):
                        x_vals_bdry.append(b_x)
                        y_vals_bdry.append(val_bdry)
                        is_bdry.append(np.array([i1,i2]))
                        js_bdry.append(np.array([j1,j2]))
            if (i1-i2 != 0):# if equation solvable
                # if max is actually attained
                if (abs(np.max(evaluate_trop_monomials((1/(i1-i2))*rhs_a_y,a_y,A))- A[(i1,j1)] - i1*(1/(i1-i2))*rhs_a_y - j1*a_y) < eps):
                    val_bdry = (1/(i1-i2))*rhs_a_y
                    # if point in plotting range
                    if (val_bdry >= a_x) and (val_bdry <= b_x):
                        x_vals_bdry.append(val_bdry)
                        y_vals_bdry.append(a_y)
                        is_bdry.append(np.array([i1,i2]))
                        js_bdry.append(np.array([j1,j2]))
                # if max is actually attained
                if (abs(np.max(evaluate_trop_monomials((1/(i1-i2))*rhs_b_y,b_y,A))- A[(i1,j1)] - i1*(1/(i1-i2))*rhs_b_y - j1*b_y) < eps):
                    val_bdry = (1/(i1-i2))*rhs_b_y
                    # if point in plotting range
                    if (val_bdry >= a_x) and (val_bdry <= b_x):
                        x_vals_bdry.append(val_bdry)
                        y_vals_bdry.append(b_y)
                        is_bdry.append(np.array([i1,i2]))
                        js_bdry.append(np.array([j1,j2]))

        '''Now we have all points at wich edges start and end'''
        # We check which points we should connect with a line
        # Sadly just checking if the point in the middle is in the curve is not sufficient.
        # The max has to be attained by the right monomials, too

        # Identify points that occur more often than once
        all_vertices = [] # List of all vertices: Tuples [x,y] of coordinates of vertices of the curve (including boundarys of the plot)
        all_vertices_exp = [] # And a list of lists of their exponents [i,j], all listed only one time

        for i in range(len(x_vals_triple)):
            index_of_triple_point = tupel_is_in_list([x_vals_triple[i], y_vals_triple[i]],all_vertices,eps)
            if index_of_triple_point == -1:
                all_vertices.append([x_vals_triple[i], y_vals_triple[i]])
                all_vertices_exp.append([[is_triple[i][0],js_triple[i][0]],[is_triple[i][1],js_triple[i][1]],[is_triple[i][2],js_triple[i][2]]])
            else:
                if not([is_triple[i][0],js_triple[i][0]] in all_vertices_exp[index_of_triple_point]):
                    all_vertices_exp[index_of_triple_point].append([is_triple[i][0],js_triple[i][0]])
                if not([is_triple[i][1],js_triple[i][1]] in all_vertices_exp[index_of_triple_point]):
                    all_vertices_exp[index_of_triple_point].append([is_triple[i][1],js_triple[i][1]])
                if not([is_triple[i][2],js_triple[i][2]] in all_vertices_exp[index_of_triple_point]):
                    all_vertices_exp[index_of_triple_point].append([is_triple[i][2],js_triple[i][2]])

        for i in range(len(x_vals_bdry)):
            index_of_triple_point = tupel_is_in_list([x_vals_bdry[i], y_vals_bdry[i]],all_vertices,eps)
            if index_of_triple_point == -1:
                all_vertices.append([x_vals_bdry[i], y_vals_bdry[i]])
                all_vertices_exp.append([[is_bdry[i][0],js_bdry[i][0]],[is_bdry[i][1],js_bdry[i][1]]])
            else:
                if not([is_bdry[i][0],js_bdry[i][0]] in all_vertices_exp[index_of_triple_point]):
                    all_vertices_exp[index_of_triple_point].append([is_bdry[i][0],js_bdry[i][0]])
                if not([is_bdry[i][1],js_bdry[i][1]] in all_vertices_exp[index_of_triple_point]):
                    all_vertices_exp[index_of_triple_point].append([is_bdry[i][1],js_bdry[i][1]])


        ''' Check which vertices to connect'''
        number_of_vertices = len(all_vertices) # The number of vertices
        x_edge_vals = [] # List of the x values, as length 2 arrays, of start and endpoint of edges of the trop variety
        y_edge_vals = [] # List of the y values, as length 2 arrays, of start and endpoint of edges of the trop variety
        edge_weights = [] # List of the weights of the edges

        x_vals_dual = [] # List of the x values of start and endpoint of edges in the dual
        y_vals_dual = [] # List of the y values of start and endpoint of edges in the dual

        # Go through pairs of vertices
        edge_count = 0
        for i in range(number_of_vertices-1):
            for j in range(i+1,number_of_vertices):
                # Compute the coordinates of the point between the vertices
                mid_x = (all_vertices[i][0] + all_vertices[j][0])/2
                mid_y = (all_vertices[i][1] + all_vertices[j][1])/2
                val = np.max(evaluate_trop_monomials(mid_x,mid_y,A))
                a = []
                b = []
                for exp in all_vertices_exp[i]: # Fill a with the values of the monomials of vertex 0 at the middle point
                    a.append(A[(exp[0],exp[1])] + exp[0] * mid_x + exp[1] * mid_y)
                for exp in all_vertices_exp[j]: # Fill b with the values of the monomials of vertex 1 at the middle point
                    b.append(A[(exp[0],exp[1])] + exp[0] * mid_x + exp[1] * mid_y)

                a_ind, b_ind = get_values_indices(a,b,val,eps)

                if (abs(max(a)-val) < eps) and (abs(max(b)-val) < eps) and (len(a_ind) > 1) and (len(b_ind) > 1): # If middle point in variety
                    # We found an edge!
                    x_edge_vals.append(np.array([all_vertices[i][0], all_vertices[j][0]]))
                    y_edge_vals.append(np.array([all_vertices[i][1], all_vertices[j][1]]))

                    # And its dual + weight
                    # The exponents which give the weight are the ones that define the dual edge. So compute all possible weights.
                    edge_weights.append(0)
                    x_vals_dual.append([])
                    y_vals_dual.append([])
                    num_max_exps = len(a_ind)

                    for k in range(0,num_max_exps-1):
                        for l in range(k+1,num_max_exps):
                            w = math.gcd(
                                abs(all_vertices_exp[i][a_ind[k]][0] - all_vertices_exp[i][a_ind[l]][0]),
                                abs(all_vertices_exp[i][a_ind[k]][1] - all_vertices_exp[i][a_ind[l]][1]))
                            if (edge_weights[edge_count] < w):
                                edge_weights[edge_count] = w
                                x_vals_dual[edge_count] = [all_vertices_exp[i][a_ind[k]][0], all_vertices_exp[i][a_ind[l]][0]]
                                y_vals_dual[edge_count] = [all_vertices_exp[i][a_ind[k]][1], all_vertices_exp[i][a_ind[l]][1]]
                    edge_count = edge_count + 1

        '''Computations for the dual'''
        # First, we need to compute the degree of the polynomial
        degF = 0
        for (i1,j1) in A:
            if (i1 + j1 > degF):
                degF = i1 + j1
                
        # Then we can create the coordinates for \Delta_degF
        x_coords_dual = []
        y_coords_dual = []

        for i in range(degF+1):
            for j in range(degF-i+1):
                x_coords_dual.append(i)
                y_coords_dual.append(j)
                
        ''' Plotting'''

        ''' Plotting the curve'''
        plt.figure(1)

        # Draw the lines
        for i in range(len(x_edge_vals)):
            plt.plot(x_edge_vals[i],y_edge_vals[i],color=color_curve)
            if (edge_weights[i] > 1):
                plt.text((x_edge_vals[i][0] + x_edge_vals[i][1]) / 2, (y_edge_vals[i][0] + y_edge_vals[i][1]) / 2, edge_weights[i], fontsize=weight_text_size, color=weight_color)

        plt.axis([a_x, b_x, a_y, b_y])  # Defines the axes

        # Rescales the plot such that the image is not distorted
        ratio = (b_x - a_x) / (b_y - a_y)
        plt.gca().set_aspect(ratio, adjustable='box')
        #plt.gcf().set_size_inches(9, 9)  # Windowsize of the plots

        # Write the polynomial as title of the plot
        plt.title(f"{write_the_polynomial(A,eps,True)}")

        # What are the axes called?

        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')

        plt.grid(True) # Makes a grid. Or doesn't. You decide(d)!

        ''' Plotting the dual'''
        plt.figure(2)

        plt.title("Dual subdivision")

        # Draw all points of \Delta_deg
        plt.scatter(x_coords_dual, y_coords_dual, color=color_NP, marker='o')

        # Draw the lines
        for i in range(len(x_vals_dual)):
            plt.plot(x_vals_dual[i],y_vals_dual[i],color=color_DS, marker = 'o')

        # Only integers at the axes. Decimals do not really make sense
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().set_aspect(1, adjustable='box')
    except Exception as e:
        print(f"Error: {e}")