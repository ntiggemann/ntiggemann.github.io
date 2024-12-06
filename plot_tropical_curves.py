###############################
# This script plots the tropical curve defined by a bivariate tropical polynomial
# and its dual subdivision
# v2.0
###############################
''' Preamble '''
# Don't change anything
import numpy as np
import matplotlib.pyplot as plt
import math # pre-installed
from matplotlib.ticker import MaxNLocator
import tkinter as tk # often pre-installed
import re # pre-installed

# Default values:
default_polynomial = "-2x2-2y2-2+x+y+xy"
default_use_max = True
default_a_x = -5 # Lower bound x-values
default_b_x = 5 # Upper bound x-values
default_a_y = -5 # Lower bound y-values
default_b_y = 5 # Upper bound y-values
default_color_curve = 'blue'
default_color_DS = 'blue'
default_color_NP = 'black'
default_make_grid = True
default_print_vertices = False
default_show_weights = True
default_size_weights = 12
default_color_weights = 'black'
default_print_tikz_curve = False
default_print_tikz_DS = False


''' If something weird happens, read this:'''
## Non-symbolic computing problems: One does not simply check floats for equality
# If you don't do sth weird, this part should not concern you. Some words:
# This eps can usually be chosen very small (~10^{-14}). If the polynomial is of very high degree you might need a bigger eps
# So if the sanity check of a plot fails, this is the first thing to change and check again.

eps = 0.00000000000001  # Epsilon - this number defines under which distance two floats are interpreted as the same number

# Auxiliar-minus-infty - just some very negative number
infty = -9999999999999999 #DO NOT CHANGE, even if you want to use the min convention!

''' Functions '''
def get_exponent(index,degree):
    '''
    Calculate the exponent of the monomial corresponding to the given index in a coefficient matrix.

    Args:
        index (int): The linear index in the coefficient matrix.
        degree (int): The maximum degree (number of rows/columns - 1) in the matrix.

    Returns:
        tuple: The row (x exponent) and column (y exponent) corresponding to the monomial.
    '''
    return index // (degree+1),(index % (degree+1))

def evaluate_trop_monomials(x,y,A,non_infty_entries,max_exponent):
    '''
    Given a point (x,y) and a tropical polynomial with coefficient matrix A and not-infty entries non_infty_entries, returns the values of the tropical 
    monomials.
    Args:
        x (int or float): x-coordinate
        y (int or float): y-coordinate
        A (2-dim Array of int or float): The coefficient matrix of a tropical polynomial
        non_infty_entries (list of int): List of the non-infty entries of A
    Returns:
        array of floats: The values in the max of the tropical polynomial
    '''
    values = np.zeros(len(non_infty_entries))

    for idx, linear_index in enumerate(non_infty_entries):
        i1, j1 = get_exponent(linear_index,max_exponent)
        values[idx] = (i1)*x + (j1)*y + A[i1,j1]
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
            
def write_the_polynomial(A,non_infty_entries,max_exponent,eps,is_title):
    '''
    Given the coefficients of a tropical polynomial, generates a string of LaTeX code for it
    Args:
        A (2-dim array of int): Coefficients of a polynomial
        non_infty_entries (list of int): list of relevant entries in the polynomial
    Returns
        string: LaTeX code that writes the polynomial
    '''
    i, j = get_exponent(non_infty_entries[0],max_exponent)

    if not is_title:
        if (abs(A[i,j]-round(A[i,j],0)) < eps):
            coeff = f"{int(round(A[i,j],0))}"
        else:
            coeff = f"{A[i,j]}"
        latex_polynomial = f"{coeff}\\odot x^{ {i} }y^{ {j} }" #f"{coeff}\\odot x^{{({i},{j})}}"
        for k in range(1,len(non_infty_entries)):
            i,j = get_exponent(non_infty_entries[k],max_exponent)
            if (abs(A[i,j]-round(A[i,j],0)) < eps):
                coeff = f"{int(round(A[i,j],0))}"
            else:
                coeff = f"{A[i,j]}"
            latex_polynomial = latex_polynomial + f"\\oplus {coeff}\\odot x^{ {i} }y^{ {j} }"
    else:
        monom_count = 0
        if (abs(A[i,j]-round(A[i,j],0)) < eps):
            coeff = f"{int(round(A[i,j],0))}"
        else:
            coeff = f"{A[i,j]}"
        latex_polynomial = f"${coeff}\\odot x^{ {i} }y^{ {j} }$" #f"{coeff}\\odot x^{{({i},{j})}}"

        for k in range(1,len(non_infty_entries)):
            i,j = get_exponent(non_infty_entries[k],max_exponent)
            if (abs(A[i,j]-round(A[i,j],0)) < eps):
                coeff = f"{int(round(A[i,j],0))}"
            else:
                coeff = f"{A[i,j]}"
            if monom_count != 0 and monom_count % 5 == 0:
                latex_polynomial = latex_polynomial + f"$\\oplus {coeff}\\odot x^{ {i} }y^{ {j} }$" + "\n"
            else:
                latex_polynomial = latex_polynomial + f"$\\oplus {coeff}\\odot x^{ {i} }y^{ {j} }$"
            monom_count = monom_count + 1
    return latex_polynomial

def change_signs_of_polynomial(A,non_infty_entries,max_exponent):
    '''
    Args:
        Matrix A: The coefficient matrix
        non_infty_entries (list of int): Some indices
    Returns:
        Matrix A with the signs of all entries indicated by non_infty_entries changed.
    '''
    for k in non_infty_entries:
        i,j = get_exponent(k,max_exponent)
        A[i,j] = (-1)*A[i,j]
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

def parse_polynomial_trop(poly_str,use_min):
    """
    Parses a polynomial string in two variables (x and y) and creates a 2D array of coefficients.

    Args:
        poly_str (str): The polynomial string, e.g., "3.2x2y + yx + 9".

    Returns:
        np.array: A 2D array representing the coefficients of the polynomial.
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
            if use_min:
                coefficients[(x_exp, y_exp)] = min(coefficient,coefficients[(x_exp, y_exp)])
            else:
                coefficients[(x_exp, y_exp)] = max(coefficient,coefficients[(x_exp, y_exp)])
        else:
            coefficients[(x_exp, y_exp)] = coefficient

    # Determine the size of the array
    max_x = max((key[0] for key in coefficients), default=0)
    max_y = max((key[1] for key in coefficients), default=0)

    max_exponent = max(max_x,max_y)

    # Create a 2D array
    result = infty * np.ones((max_exponent + 1, max_exponent + 1), dtype=float)

    # Fill the array with coefficients
    for (x_exp, y_exp), coeff in coefficients.items():
        result[x_exp, y_exp] = coeff

    return result, max_exponent

# Function to plot the zero set
def plot_trop_polynomial(poly_str,a_x,b_x,a_y,b_y,use_min,
                        show_weights,weight_text_size,weight_color,
                        color_curve,color_DS,color_NP,make_grid,write_polynomial,axes_xy,
                        print_vertices,print_tikz_code_curve,print_tikz_code_dual,tikz_color_curve,tikz_color_DS,tikz_color_NP,eps):
    try:
        A, max_exponent = parse_polynomial_trop(poly_str,use_min)

        # Get the relevant entries of the matrix.
        non_infty_entries = []
        for k in range((max_exponent+1)**2):
            i,j = get_exponent(k,max_exponent)
            if not (A[i,j] == infty):
                non_infty_entries.append(k)

        # The number of monomials
        number_of_monomials = len(non_infty_entries)

        if (number_of_monomials < 2):
            raise ValueError("Please enter a polynomial with at least two different monomials.")

        if use_min:
            A = change_signs_of_polynomial(A,non_infty_entries,max_exponent)

            temp = a_x
            a_x = (-1) * b_x
            b_x = (-1) * temp

            temp = a_y
            a_y = (-1) * b_y
            b_y = (-1) * temp


        # Get ready for linear algebra!
        '''Compute the vertices of the curve, aka triple-max points'''

        x_vals_triple = [] # List of the x values for triple-max points
        y_vals_triple = [] # List of the y values for triple-max points
        is_triple = [] # List of length 3 arrays, the x-exponents of the maximal monomials
        js_triple = [] # List of length 3 arrays, the y-exponents of the maximal monomials

        # Walk through matrix A, check all triples of monomials

        for i in range(number_of_monomials-2):
            for j in range(i+1,number_of_monomials-1):
                for k in range(j+1,number_of_monomials):
                    ## Solve for points where a value is attained three times
                    # Recovers the exponents from the linear counting

                    i1,j1 = get_exponent(non_infty_entries[i],max_exponent)
                    i2,j2 = get_exponent(non_infty_entries[j],max_exponent)
                    i3,j3 = get_exponent(non_infty_entries[k],max_exponent)

                    # The coefficients of the system of linear equations we have to solve
                    B = np.array([[i1-i2, j1-j2],[i1-i3, j1-j3]])
                    b = np.array([A[i2,j2] - A[i1,j1],A[i3,j3] - A[i1,j1]])
                    # Proceed only if  the syst of lin eq is uniquely solvable
                    if (0 != np.linalg.det(B)):
                        x = np.linalg.solve(B, b)
                        # Check if it is actually a maximum
                        if (abs(np.max(evaluate_trop_monomials(x[0],x[1],A,non_infty_entries,max_exponent)) - A[i1,j1] - i1*x[0] - j1*x[1]) < eps):
                            x_vals_triple.append(x[0])
                            y_vals_triple.append(x[1])
                            is_triple.append(np.array([i1,i2,i3]))
                            js_triple.append(np.array([j1,j2,j3]))
                            if (x[0] <= a_x) or (x[0] >= b_x) or (x[1] <= a_y) or (x[1] >= b_y):
                                # If a vertex of the curve is not in the plot, it will look weird, and you will probably be unhappy with the plot
                                print(f"The point ({x[0]},{x[1]}) is not in your plotting range, which could be a problem...")

        ''' Compute the intersection points of the curve with the boundary of the plot'''

        x_vals_bdry = [] # List of the x values for points in the trop variety on the boundary
        y_vals_bdry = [] # List of the y values for points in the trop variety on the boundary
        is_bdry = [] # List of length 2 arrays, the x-exponents of the maximal monomials
        js_bdry = [] # List of length 2 arrays, the y-exponents of the maximal monomials

        # Walk through matrix A, check all pairs of monomials
        for i in range(number_of_monomials-1):
            for j in range(i+1,number_of_monomials):
                ## Solve for points on boundary where a value is attained twice
                # Recovers the exponents from the linear counting
                i1,j1 = get_exponent(non_infty_entries[i],max_exponent)
                i2,j2 = get_exponent(non_infty_entries[j],max_exponent)

                rhs_a_x = A[i2,j2] + (i2 - i1)* a_x - A[i1,j1]
                rhs_b_x = A[i2,j2] + (i2 - i1)* b_x - A[i1,j1]
                rhs_a_y = A[i2,j2] + (j2 - j1)* a_y - A[i1,j1]
                rhs_b_y = A[i2,j2] + (j2 - j1)* b_y - A[i1,j1]
                if (j1-j2 != 0): # if equation solvable
                    # if max is actually attained
                    if (abs(np.max(evaluate_trop_monomials(a_x,(1/(j1-j2))*rhs_a_x,A,non_infty_entries,max_exponent))- A[i1,j1] - i1*a_x - j1*(1/(j1-j2))*rhs_a_x) < eps):
                        val_bdry = (1/(j1-j2))*rhs_a_x
                        # if point in plotting range
                        if (val_bdry >= a_y) and (val_bdry <= b_y):
                            x_vals_bdry.append(a_x)
                            y_vals_bdry.append(val_bdry)
                            is_bdry.append(np.array([i1,i2]))
                            js_bdry.append(np.array([j1,j2]))
                    # if max is actually attained
                    if (abs(np.max(evaluate_trop_monomials(b_x,(1/(j1-j2))*rhs_b_x,A,non_infty_entries,max_exponent)) - A[i1,j1] - i1*b_x - j1*(1/(j1-j2))*rhs_b_x) < eps):
                        val_bdry = (1/(j1-j2))*rhs_b_x
                        # if point in plotting range
                        if (val_bdry >= a_y) and (val_bdry <= b_y):
                            x_vals_bdry.append(b_x)
                            y_vals_bdry.append(val_bdry)
                            is_bdry.append(np.array([i1,i2]))
                            js_bdry.append(np.array([j1,j2]))
                if (i1-i2 != 0):# if equation solvable
                    # if max is actually attained
                    if (abs(np.max(evaluate_trop_monomials((1/(i1-i2))*rhs_a_y,a_y,A,non_infty_entries,max_exponent))- A[i1,j1] - i1*(1/(i1-i2))*rhs_a_y - j1*a_y) < eps):
                        val_bdry = (1/(i1-i2))*rhs_a_y
                        # if point in plotting range
                        if (val_bdry >= a_x) and (val_bdry <= b_x):
                            x_vals_bdry.append(val_bdry)
                            y_vals_bdry.append(a_y)
                            is_bdry.append(np.array([i1,i2]))
                            js_bdry.append(np.array([j1,j2]))
                    # if max is actually attained
                    if (abs(np.max(evaluate_trop_monomials((1/(i1-i2))*rhs_b_y,b_y,A,non_infty_entries,max_exponent))- A[i1,j1] - i1*(1/(i1-i2))*rhs_b_y - j1*b_y) < eps):
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
                val = np.max(evaluate_trop_monomials(mid_x,mid_y,A,non_infty_entries,max_exponent))
                a = []
                b = []
                for exp in all_vertices_exp[i]: # Fill a with the values of the monomials of vertex 0 at the middle point
                    a.append(A[exp[0],exp[1]] + exp[0] * mid_x + exp[1] * mid_y)
                for exp in all_vertices_exp[j]: # Fill b with the values of the monomials of vertex 1 at the middle point
                    b.append(A[exp[0],exp[1]] + exp[0] * mid_x + exp[1] * mid_y)

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
        for i in range(number_of_monomials):
            i1,j1 = get_exponent(non_infty_entries[i],max_exponent)
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

        if use_min:
            for i in range(len(x_edge_vals)):
                x_edge_vals[i] = (-1) * x_edge_vals[i]
                y_edge_vals[i] = (-1) * y_edge_vals[i]

        ''' Plotting the curve'''
        plt.figure(1)

        # Draw the lines
        if print_tikz_code_curve:
            print("% --- tikz code for the curve ---")
            print("%" + "$F = " + f"{write_the_polynomial(A,non_infty_entries,max_exponent,eps,False)}" + "$")
            for i in range(len(x_edge_vals)):
                plt.plot(x_edge_vals[i],y_edge_vals[i],color=color_curve)
                if show_weights and (edge_weights[i] > 1):
                    plt.text((x_edge_vals[i][0] + x_edge_vals[i][1]) / 2, (y_edge_vals[i][0] + y_edge_vals[i][1]) / 2, edge_weights[i], fontsize=weight_text_size, color=weight_color)
                    print(f"\\node[{weight_color}] at ({(x_edge_vals[i][0] + x_edge_vals[i][1]) / 2},{(y_edge_vals[i][0] + y_edge_vals[i][1]) / 2}) {{{edge_weights[i]}}};")
                print(f"\\draw[{tikz_color_curve}] ({x_edge_vals[i][0]},{y_edge_vals[i][0]}) -- ({x_edge_vals[i][1]},{y_edge_vals[i][1]});")
        else:
            for i in range(len(x_edge_vals)):
                plt.plot(x_edge_vals[i],y_edge_vals[i],color=color_curve)
                if show_weights and (edge_weights[i] > 1):
                    plt.text((x_edge_vals[i][0] + x_edge_vals[i][1]) / 2, (y_edge_vals[i][0] + y_edge_vals[i][1]) / 2, edge_weights[i], fontsize=weight_text_size, color=weight_color)

        if use_min: 
            temp = a_x
            a_x = (-1) * b_x
            b_x = (-1) * temp

            temp = a_y
            a_y = (-1) * b_y
            b_y = (-1) * temp

        plt.axis([a_x, b_x, a_y, b_y])  # Defines the axes

        # Rescales the plot such that the image is not distorted
        ratio = (b_x - a_x) / (b_y - a_y)
        plt.gca().set_aspect(ratio, adjustable='box')
        #plt.gcf().set_size_inches(9, 9)  # Windowsize of the plots

        # Write the polynomial?
        if write_polynomial:
            if use_min:
                A = change_signs_of_polynomial(A,non_infty_entries,max_exponent)
            plt.title("$F = $" + f"{write_the_polynomial(A,non_infty_entries,max_exponent,eps,True)}")

        # What are the axes called?
        if axes_xy:
            plt.xlabel('$x$')
            plt.ylabel('$y$')
        else:
            plt.xlabel(r'$x_1$')
            plt.ylabel(r'$x_2$')

        plt.grid(make_grid) # Makes a grid. Or doesn't. You decide(d)!

        # Should we give the interesting points?
        if print_vertices:
            # If yes, then here we go!
            relevant_pts = []
            redundant_pts = []
            for i in range(len(x_edge_vals)):
                redundant_pts.append([x_edge_vals[i][0],y_edge_vals[i][0]])
                redundant_pts.append([x_edge_vals[i][1],y_edge_vals[i][1]])
            for i in range(len(redundant_pts)):
                # This eliminates the points that occure more often than once
                if not any(np.array_equal(redundant_pts[i], bi) for bi in relevant_pts):
                    relevant_pts.append(redundant_pts[i])
            # Sort by x and y coordinate, in this order
            relevant_pts = sorted(relevant_pts, key=lambda x: (x[0], x[1]))
            print("Vertices & intersection points with boundary:")
            for i in range(len(relevant_pts)):
                print(f"({relevant_pts[i][0]}, " + f"{relevant_pts[i][1]})")

        ''' Plotting the dual'''
        plt.figure(2)

        plt.title("Dual subdivision")

        # Draw all points of \Delta_deg
        plt.scatter(x_coords_dual, y_coords_dual, color=color_NP, marker='o')


        if print_tikz_code_dual:
            print("%--- tikz code for the dual ---")
            for i in range(len(x_coords_dual)):
                print(f"\\filldraw[color = {tikz_color_NP}, fill={tikz_color_NP}] ({x_coords_dual[i]},{y_coords_dual[i]}) circle (1.5pt);")
            # Draw the lines & print the code
            for i in range(len(x_vals_dual)):
                plt.plot(x_vals_dual[i],y_vals_dual[i],color=color_DS, marker = 'o')    
                print(f"\\draw[{tikz_color_DS}] ({x_vals_dual[i][0]},{y_vals_dual[i][0]}) -- ({x_vals_dual[i][1]},{y_vals_dual[i][1]});")
        else:
            # Draw the lines
            for i in range(len(x_vals_dual)):
                plt.plot(x_vals_dual[i],y_vals_dual[i],color=color_DS, marker = 'o')

        # Only integers at the axes. Decimals do not really make sense
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

        # Tadaaaaaaa!
        plt.show()
    except Exception as e:
        print(f"Error: {e}")

# Tkinter GUI
def main():
    def on_plot():
        plt.close('all')
        poly_str = entry.get() # Get polynomial
        # Get boundaries
        try:
            a = float(entry_a.get())
            b = float(entry_b.get())
            c = float(entry_c.get())
            d = float(entry_d.get())
            if a >= b or c >= d:
                raise ValueError("Invalid boundaries: Ensure a < b and c < d.")
        except ValueError as e:
            print(f"Error: {e}")
            return

        # Get use_min
        use_min = not use_max_var.get()
        # Get show_weights
        show_weights = show_weights_var.get()
        # Get weight_size
        weight_text_size = entry_weight_size.get()
        # Get weight_color
        weight_color = entry_weight_color.get()
        # Get color_curve
        color_curve = entry_color_curve.get()
        # Get color_DS
        color_DS = entry_color_DS.get()
        # Get color_NP
        color_NP = entry_color_NP.get()
        # Get make_grid
        make_grid = make_grid_var.get()
        # Get print_vertices
        print_vertices = print_vertices_var.get()

        write_polynomial = True # This stays true, as it is probably good for checking the input
        axes_xy = True # If one would want to set this to true, one would have to change the write_the_polnomial function

        # Get print_tikz_code_curve
        print_tikz_code_curve = print_tikz_code_curve_var.get()
        # Get print_tikz_code_dual
        print_tikz_code_dual = print_tikz_code_dual_var.get()
        
        tikz_color_curve, tikz_color_DS, tikz_color_NP = color_curve, color_DS, color_NP # This was an unnecessary feature.

        plot_trop_polynomial(poly_str,a,b,c,d,use_min,
                                show_weights,weight_text_size,weight_color,
                                color_curve,color_DS,color_NP,make_grid,write_polynomial,axes_xy,
                                print_vertices,
                                print_tikz_code_curve,print_tikz_code_dual,tikz_color_curve,tikz_color_DS,tikz_color_NP,
                                eps)

    # Create main window
    window = tk.Tk()
    window.title("Plotter for tropical curves")

    # Ensure proper closure of the application
    def on_close():
        window.destroy()
        tk.Tk().destroy()
        print("Application closed.")

    window.protocol("WM_DELETE_WINDOW", on_close)

    # Input field for the polynomial
    frame = tk.Frame(window)
    frame.pack(pady=10)

    tk.Label(frame, text="Enter Polynomial:").pack(side=tk.LEFT)
    entry = tk.Entry(frame, width=50)
    entry.insert(0, default_polynomial)  # Standard value for polynomial
    entry.pack(side=tk.LEFT, padx=5)

    # Button to plot
    plot_button = tk.Button(window, text="Plot", command=on_plot)
    plot_button.pack(pady=10)

    # Checkbox for use_min
    frame_use_max = tk.Frame(window)
    frame_use_max.pack(pady=10)

    use_max_var = tk.BooleanVar(value=default_use_max)
    use_max_checkbox = tk.Checkbutton(
        frame_use_max, text="Use max-convention", variable=use_max_var
    )
    use_max_checkbox.pack(side=tk.LEFT)

    # Input fields for plot boundaries
    frame_bounds = tk.Frame(window)
    frame_bounds.pack(pady=10)

    tk.Label(frame_bounds, text="x-range:").pack(side=tk.LEFT)
    entry_a = tk.Entry(frame_bounds, width=5)
    entry_a.insert(0, default_a_x)  # Standard value for a
    entry_a.pack(side=tk.LEFT, padx=2)
    entry_b = tk.Entry(frame_bounds, width=5)
    entry_b.insert(0, default_b_x)  # Standard value for b
    entry_b.pack(side=tk.LEFT, padx=2)

    tk.Label(frame_bounds, text="y-range:").pack(side=tk.LEFT, padx=10)
    entry_c = tk.Entry(frame_bounds, width=5)
    entry_c.insert(0, default_a_y)  # Standard value for c
    entry_c.pack(side=tk.LEFT, padx=2)
    entry_d = tk.Entry(frame_bounds, width=5)
    entry_d.insert(0, default_b_y)  # Standard value for d
    entry_d.pack(side=tk.LEFT, padx=2)

    # Options for plot
    instruction = tk.Label(
            window,
            text=(
                "Options for the plot\n"
                "_______________"
            ),
            font=("Arial", 12),
            justify="center",  # Text zentrieren
            wraplength=500     # Max. Breite für Zeilenumbruch (in Pixeln)
        )

    instruction.pack(pady=10)  # Abstände zum Rand hinzufügen

    # Input field for the color of the curve
    frame_color_curve = tk.Frame(window)
    frame_color_curve.pack(pady=10)

    tk.Label(frame_color_curve, text="Color of curve:").pack(side=tk.LEFT)
    entry_color_curve = tk.Entry(frame_color_curve, width=50)
    entry_color_curve.insert(0, default_color_curve)  # Standard value for color of the curve
    entry_color_curve.pack(side=tk.LEFT, padx=5)

    # Input field for the color of the DS
    frame_color_DS = tk.Frame(window)
    frame_color_DS.pack(pady=10)

    tk.Label(frame_color_DS, text="Color of dual subdivision:").pack(side=tk.LEFT)
    entry_color_DS = tk.Entry(frame_color_DS, width=50)
    entry_color_DS.insert(0, default_color_DS)  # Standard value for color of the curve
    entry_color_DS.pack(side=tk.LEFT, padx=5)


    # Input field for the color of the Newton polygon
    frame_color_NP = tk.Frame(window)
    frame_color_NP.pack(pady=10)

    tk.Label(frame_color_NP, text="Color of Newton polygon:").pack(side=tk.LEFT)
    entry_color_NP = tk.Entry(frame_color_NP, width=50)
    entry_color_NP.insert(0, default_color_NP)  # Standard value for color of the curve
    entry_color_NP.pack(side=tk.LEFT, padx=5)

    # Checkbox for make_grid
    frame_make_grid = tk.Frame(window)
    frame_make_grid.pack(pady=10)

    make_grid_var = tk.BooleanVar(value=default_make_grid)
    make_grid_checkbox = tk.Checkbutton(
        frame_make_grid, text="Draw a grid", variable=make_grid_var
    )
    make_grid_checkbox.pack(side=tk.LEFT)

    # Checkbox for print_vertices
    frame_print_vertices = tk.Frame(window)
    frame_print_vertices.pack(pady=10)

    print_vertices_var = tk.BooleanVar(value=default_print_vertices)
    print_vertices_checkbox = tk.Checkbutton(
        frame_print_vertices, text="Print ends of edges", variable=print_vertices_var
    )
    print_vertices_checkbox.pack(side=tk.LEFT)

    # Options for weights
    instruction = tk.Label(
            window,
            text=(
                "Options for weights\n"
                "_______________"
            ),
            font=("Arial", 12),
            justify="center",  # Text zentrieren
            wraplength=500     # Max. Breite für Zeilenumbruch (in Pixeln)
        )

    instruction.pack(pady=10)  # Abstände zum Rand hinzufügen

    # Checkbox for show_weights
    frame_show_weights = tk.Frame(window)
    frame_show_weights.pack(pady=10)

    show_weights_var = tk.BooleanVar(value=True)
    show_weights_checkbox = tk.Checkbutton(
        frame_show_weights, text="Show weights of edges", variable=show_weights_var
    )
    show_weights_checkbox.pack(side=tk.LEFT)

    # Input field for the size of the weights
    frame_weight_size = tk.Frame(window)
    frame_weight_size.pack(pady=10)

    tk.Label(frame_weight_size, text="Size of weights:").pack(side=tk.LEFT)
    entry_weight_size = tk.Entry(frame_weight_size, width=50)
    entry_weight_size.insert(0, default_size_weights)  # Standard value for weight_size
    entry_weight_size.pack(side=tk.LEFT, padx=5)

    # Input field for the color of the weights
    frame_weight_color = tk.Frame(window)
    frame_weight_color.pack(pady=10)

    tk.Label(frame_weight_color, text="Color of weights:").pack(side=tk.LEFT)
    entry_weight_color = tk.Entry(frame_weight_color, width=50)
    entry_weight_color.insert(0, default_color_weights)  # Standard value for weight_color
    entry_weight_color.pack(side=tk.LEFT, padx=5)
    
    # Options for tikz
    instruction = tk.Label(
            window,
            text=(
                "Options for tikz code\n"
                "_______________"
            ),
            font=("Arial", 12),
            justify="center",  # Text zentrieren
            wraplength=500     # Max. Breite für Zeilenumbruch (in Pixeln)
        )

    instruction.pack(pady=10)  # Abstände zum Rand hinzufügen

    # Checkbox for print_tikz_code_curve
    frame_print_tikz_code_curve = tk.Frame(window)
    frame_print_tikz_code_curve.pack(pady=10)

    print_tikz_code_curve_var = tk.BooleanVar(value=default_print_tikz_curve)
    print_tikz_code_curve_checkbox = tk.Checkbutton(
        frame_print_tikz_code_curve, text="Print tikz code curve", variable=print_tikz_code_curve_var
    )
    print_tikz_code_curve_checkbox.pack(side=tk.LEFT)

    # Checkbox for print_tikz_code_dual
    frame_print_tikz_code_dual = tk.Frame(window)
    frame_print_tikz_code_dual.pack(pady=10)

    print_tikz_code_dual_var = tk.BooleanVar(value=default_print_tikz_DS)
    print_tikz_code_dual_checkbox = tk.Checkbutton(
        frame_print_tikz_code_dual, text="Print tikz code dual", variable=print_tikz_code_dual_var
    )
    print_tikz_code_dual_checkbox.pack(side=tk.LEFT)

    # Start GUI loop
    window.mainloop()
if __name__ == "__main__":
    main()