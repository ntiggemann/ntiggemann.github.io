###############################
# This script plots the tropical curve defined by a bivariate tropical polynomial
# and its dual subdivision
# v2.3.1
###############################


''' Preamble '''
# Don't change anything
import numpy as np
import matplotlib.pyplot as plt
import math # pre-installed
from matplotlib.ticker import MaxNLocator
import tkinter as tk # often pre-installed
import re # pre-installed
import itertools
from tkinter import messagebox
from tkinter import ttk
import webbrowser

# Default values:
default_polynomial = "-2x2-2y2-2+x+y+xy"
default_use_max = True
default_auto_adjust = True
default_auto_square = True
default_auto_boundary_dist = 1
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

def parse_polynomial_trop(poly_str,use_min):
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
            messagebox.showerror("Invalid term",f"Invalid term: {term} in the polynomial input field")
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

    # Neues Dictionary mit sortierten Schlüsseln
    sorted_coefficients = dict(sorted(coefficients.items(), key=lambda x: x[0], reverse=True))

    return sorted_coefficients

def common_sublists(list1, list2):
    """
    Finds all sublists that are present in both input lists.
    
    :param list1: List of lists, where each sublist has two integers
    :param list2: List of lists, where each sublist has two integers
    :return: A list of sublists that are present in both list1 and list2
    """
    set1 = {tuple(sublist) for sublist in list1}  # Convert to set of tuples
    set2 = {tuple(sublist) for sublist in list2}
    
    common = set1.intersection(set2)  # Find common tuples
    
    return [list(sublist) for sublist in common]  # Convert back to list of lists

# Function to plot the zero set
def plot_trop_polynomial(poly_str,a_x,b_x,a_y,b_y,use_min,
                        show_weights,weight_text_size,weight_color,
                        color_curve,color_DS,color_NP,make_grid,
                        print_vertices,print_tikz_code_curve,print_tikz_code_dual,auto_adjust,auto_boundary_dist,auto_square,eps):
    try:
        A = parse_polynomial_trop(poly_str,use_min)

        # The number of monomials
        number_of_monomials = len(A)

        if (number_of_monomials < 2):
            messagebox.showerror("Missing vertices", "Please enter a polynomial with at least two different monomials.")
            raise ValueError("Please enter a polynomial with at least two different monomials.")

        if use_min:
            A = change_signs_of_polynomial(A)

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

        if x_vals_triple:
            if auto_adjust:
                # Automatic adjustment of the boundary
                a_x = min(x_vals_triple)
                b_x = max(x_vals_triple)
                a_y = min(y_vals_triple)
                b_y = max(y_vals_triple)
            else:
                if (a_x >= b_x) or (a_y >= b_y):
                    messagebox.showerror("Invalid boundaries", "Your boundaries do not make sense.")
                    raise ValueError("Your boundaries do not make sense.")

                missing_points = []
                for i in range(len(x_vals_triple)):
                    x = [x_vals_triple[i],y_vals_triple[i]]
                    if (x[0] <= a_x) or (x[0] >= b_x) or (x[1] <= a_y) or (x[1] >= b_y):
                        # If a vertex of the curve is not in the plot, it will look weird, and you will probably be unhappy with the plot
                        missing_points.append([x[0],x[1]])
                if missing_points:
                    error_message = "The points\n"
                    for vert in missing_points:
                        error_message = error_message + f"({vert[0]},{vert[1]})\n"
                    error_message = error_message + "are vertices of the curve and not in the plotting range, which will most likely ruin the plot. Adjust the plotting range or opt for auto-adjustment."
                    messagebox.showerror("Missing vertices", error_message)
        else:
            if auto_adjust:
                pseudo_vertices = [] # In here go the points that we compute below
                for ((i1,j1),(i2,j2)) in itertools.combinations(A,2):
                    # The equation where they agree is given by (i1-i2)x+(j1-j2)y+A[(i1,j1)]-A[(i2,j2)] = 0
                    
                    # If i1-i2 not= 0: one point is: (...,0)
                    if (i1-i2 != 0):
                        #Check if max point
                        x_value = (A[(i2,j2)] - A[(i1,j1)]) / (i1-i2)
                        y_value = 0
                        # If its a point of the "curve"
                        if abs(max(evaluate_trop_monomials(x_value,y_value,A)) - A[(i1,j1)] - i1*x_value) < eps:
                            pseudo_vertices.append([x_value,y_value])
                    # Else (0,...)
                    else:
                        #Check if max point
                        x_value = 0
                        y_value = (A[(i2,j2)] - A[(i1,j1)]) / (j1-j2)
                        # If its a point of the "curve"
                        if abs(max(evaluate_trop_monomials(x_value,y_value,A)) - A[(i1,j1)] - i1*x_value) < eps:
                            pseudo_vertices.append([x_value,y_value])

                # Min/Max for x values
                a_x = min(pseudo_vertices, key=lambda x: x[0])[0]
                b_x = max(pseudo_vertices, key=lambda x: x[0])[0]

                # Min/Max for y values
                a_y = min(pseudo_vertices, key=lambda x: x[1])[1]
                b_y = max(pseudo_vertices, key=lambda x: x[1])[1]

        if auto_adjust and auto_square:
            a_x_new = ((a_x+b_x)/2) - (max(b_x - a_x, b_y - a_y)/2) - auto_boundary_dist
            b_x_new = ((a_x+b_x)/2) + (max(b_x - a_x, b_y - a_y)/2) + auto_boundary_dist
            a_y_new = ((a_y+b_y)/2) - (max(b_x - a_x, b_y - a_y)/2) - auto_boundary_dist
            b_y_new = ((a_y+b_y)/2) + (max(b_x - a_x, b_y - a_y)/2) + auto_boundary_dist

            a_x = a_x_new
            b_x = b_x_new
            a_y = a_y_new
            b_y = b_y_new
        elif auto_adjust:
            a_x = a_x - auto_boundary_dist
            b_x = b_x + auto_boundary_dist
            a_y = a_y - auto_boundary_dist
            b_y = b_y + auto_boundary_dist

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
                # If the points should get connected
                common_exps = common_sublists(all_vertices_exp[i],all_vertices_exp[j])
                if len(common_exps) >= 2:
                    # Add the start and end points of the edge
                    x_edge_vals.append(np.array([all_vertices[i][0], all_vertices[j][0]]))
                    y_edge_vals.append(np.array([all_vertices[i][1], all_vertices[j][1]]))

                    # And its dual + weight
                    # The exponents which give the weight are the ones that define the dual edge. So compute all possible weights.
                    edge_weights.append(0)
                    x_vals_dual.append([])
                    y_vals_dual.append([])

                    for ((i1,j1),(i2,j2)) in itertools.combinations(common_exps,2):
                        w = math.gcd(abs(i1-i2),abs(j1-j2))
                        if (edge_weights[edge_count] < w):
                            edge_weights[edge_count] = w
                            x_vals_dual[edge_count] = [i1,i2]
                            y_vals_dual[edge_count] = [j1,j2]
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

        if use_min:
            for i in range(len(x_edge_vals)):
                x_edge_vals[i] = (-1) * x_edge_vals[i]
                y_edge_vals[i] = (-1) * y_edge_vals[i]

        ''' Plotting the curve'''
        plt.figure(1)

        # Draw the lines
        if print_tikz_code_curve:
            print("% --- tikz code for the curve ---")
            print("%" + "$F = " + f"{write_the_polynomial(A,eps,False)}" + "$")
            for i in range(len(x_edge_vals)):
                plt.plot(x_edge_vals[i],y_edge_vals[i],color=color_curve)
                if show_weights and (edge_weights[i] > 1):
                    plt.text((x_edge_vals[i][0] + x_edge_vals[i][1]) / 2, (y_edge_vals[i][0] + y_edge_vals[i][1]) / 2, edge_weights[i], fontsize=weight_text_size, color=weight_color)
                    print(f"\\node[{weight_color}] at ({(x_edge_vals[i][0] + x_edge_vals[i][1]) / 2},{(y_edge_vals[i][0] + y_edge_vals[i][1]) / 2}) {{{edge_weights[i]}}};")
                print(f"\\draw[{color_curve}] ({x_edge_vals[i][0]},{y_edge_vals[i][0]}) -- ({x_edge_vals[i][1]},{y_edge_vals[i][1]});")
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

        # Makes that the plot is not distorted
        plt.gca().set_aspect('equal')

        #plt.gcf().set_size_inches(9, 9)  # Windowsize of the plots

        # Write the polynomial as title of the plot
        if use_min:
            A = change_signs_of_polynomial(A)
        plt.title(f"{write_the_polynomial(A,eps,True)}")

        # What are the axes called?
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')

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
            # Draw the lines & print the code
            for i in range(len(x_vals_dual)):
                plt.plot(x_vals_dual[i],y_vals_dual[i],color=color_DS, marker = 'o')    
                print(f"\\draw[{color_DS}] ({x_vals_dual[i][0]},{y_vals_dual[i][0]}) -- ({x_vals_dual[i][1]},{y_vals_dual[i][1]});")
            for i in range(len(x_coords_dual)):
                print(f"\\filldraw[color = {color_NP}, fill={color_NP}] ({x_coords_dual[i]},{y_coords_dual[i]}) circle (1.5pt);")
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
    def on_plot(event=None):
        plt.close('all')
        poly_str = entry.get() # Get polynomial
        # Get boundaries
        try:
            a = float(entry_a.get())
            b = float(entry_b.get())
            c = float(entry_c.get())
            d = float(entry_d.get())
        except ValueError as e:
            messagebox.showerror("Invalid input", "Please enter numbers for the plotting range (e.g. 1 or 1.5).")
            raise ValueError("Please enter numbers for the plotting range (e.g. 1 or 1.5).")

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

        auto_adjust = auto_adjust_var.get()
        try:
            auto_boundary_dist = float(auto_boundary_dist_entry.get())
        except:
            messagebox.showerror("Invalid input", "Please enter a number for the minimal distance to the boundary of the plot (e.g. 1 or 1.5).")
            raise ValueError("Please enter a number for the minimal distance to the boundary of the plot (e.g. 1 or 1.5).")
        auto_square = auto_square_var.get()

        # Get print_tikz_code_curve
        print_tikz_code_curve = print_tikz_code_curve_var.get()
        # Get print_tikz_code_dual
        print_tikz_code_dual = print_tikz_code_dual_var.get()
        
        print("Started plotting curve defined by " + poly_str)
        plot_trop_polynomial(poly_str,a,b,c,d,use_min,
                                show_weights,weight_text_size,weight_color,
                                color_curve,color_DS,color_NP,make_grid,
                                print_vertices,
                                print_tikz_code_curve,print_tikz_code_dual,auto_adjust,auto_boundary_dist,auto_square,
                                eps)

    def on_configure(event):
        """Matches scrolling region to size of frame."""
        canvas.configure(scrollregion=canvas.bbox("all"))

    def on_mouse_wheel(event):
        """Enables scrolling with mouse wheel."""
        canvas.yview_scroll(-1 * (event.delta // 120), "units")
    
    # Create main window
    window = tk.Tk()
    window.title("Plotter for tropical curves")
    window.geometry("450x1000")

    # Create canvas and scrollbar
    canvas = tk.Canvas(window)
    scrollbar = tk.Scrollbar(window, orient="vertical",command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    # Place scroll frame in canvas
    scrollable_window = canvas.create_window((0,0), window=scrollable_frame, anchor="nw")
    
    # Update scroll region if frame changes
    scrollable_frame.bind("<Configure>", on_configure)

    # Connect scrollbar with canvas
    canvas.configure(yscrollcommand=scrollbar.set)

    # Place widgets
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    # Add scrolling with mouse wheel
    canvas.bind_all("<MouseWheel>", on_mouse_wheel)

    # Ensure proper closure of the application
    def on_close():
        plt.close('all')
        window.destroy()
        tk.Tk().destroy()
        print("Application closed.")

    window.protocol("WM_DELETE_WINDOW", on_close)

    # Input field for the polynomial
    frame = tk.Frame(scrollable_frame, borderwidth=2, relief="ridge", padx=10, pady=5)
    frame.pack(pady=10)

    tk.Label(frame, text="Enter Polynomial:").pack(side=tk.LEFT)
    entry = tk.Entry(frame, width=50)
    entry.insert(0, default_polynomial)  # Standard value for polynomial
    entry.pack(side=tk.LEFT, padx=5)
    entry.bind("<Return>", on_plot)

    # Button to plot
    plot_button = tk.Button(scrollable_frame, text="Plot", width= 8, command=on_plot)
    plot_button.pack(pady=10)

    # Checkbox for use_min
    frame_use_max = tk.Frame(scrollable_frame)
    frame_use_max.pack(pady=10)

    use_max_var = tk.BooleanVar(value=default_use_max)
    use_max_checkbox = tk.Checkbutton(
        frame_use_max, text="Use max-convention", variable=use_max_var
    )
    use_max_checkbox.pack(side=tk.LEFT)

    # Checkboxes for auto adjust
    frame_use_auto_adjust = tk.Frame(scrollable_frame)
    frame_use_auto_adjust.pack(pady=10)

    auto_adjust_var = tk.BooleanVar(value=default_auto_adjust)
    use_auto_adjust_checkbox = tk.Checkbutton(
        frame_use_auto_adjust, text="Auto-adjust range", variable=auto_adjust_var
    )
    use_auto_adjust_checkbox.pack(side=tk.LEFT)

    auto_square_var = tk.BooleanVar(value=default_auto_square)
    use_auto_square_checkbox = tk.Checkbutton(
        frame_use_auto_adjust, text="make it a square", variable=auto_square_var
    )
    use_auto_square_checkbox.pack(side=tk.LEFT)

    frame_auto_boundary_dist = tk.Frame(scrollable_frame)
    frame_auto_boundary_dist.pack(pady=10)
    tk.Label(frame_auto_boundary_dist, text="Min distance vertices to boundary of plot:").pack(side=tk.LEFT, padx=10)
    auto_boundary_dist_entry = tk.Entry(frame_auto_boundary_dist, width=5)
    auto_boundary_dist_entry.insert(0, default_auto_boundary_dist)
    auto_boundary_dist_entry.pack(side=tk.LEFT, padx=5)

    # Separator
    separator = ttk.Separator(scrollable_frame, orient="horizontal")
    separator.pack(fill="x", padx=10, pady=5)

    # Options for plot
    instruction = tk.Label(
            scrollable_frame,
            text=(
                "Manual plotting range"
            ),
            font=("Arial", 12),
            justify="center",  # Center text
            wraplength=500     # Max width for new line (in pixels)
        )

    instruction.pack(pady=5)  # Add distances to boundary

    # Input fields for plot boundaries
    frame_bounds = tk.Frame(scrollable_frame)
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

    # Separator
    separator = ttk.Separator(scrollable_frame, orient="horizontal")
    separator.pack(fill="x", padx=10, pady=5)

    # Options for plot
    instruction = tk.Label(
            scrollable_frame,
            text=(
                "Options for the plot"
            ),
            font=("Arial", 12),
            justify="center",  # Center text
            wraplength=500     # Max width for new line (in pixels)
        )

    instruction.pack(pady=5)  # Add distances to boundary

    # Input field for the color of the curve
    frame_color_curve = tk.Frame(scrollable_frame)
    frame_color_curve.pack(pady=10)

    tk.Label(frame_color_curve, text="Color of curve:").pack(side=tk.LEFT)
    entry_color_curve = tk.Entry(frame_color_curve, width=10)
    entry_color_curve.insert(0, default_color_curve)  # Standard value for color of the curve
    entry_color_curve.pack(side=tk.LEFT, padx=5)

    # Input field for the color of the DS
    frame_color_DS = tk.Frame(scrollable_frame)
    frame_color_DS.pack(pady=10)

    tk.Label(frame_color_DS, text="Color of dual subdivision:").pack(side=tk.LEFT)
    entry_color_DS = tk.Entry(frame_color_DS, width=10)
    entry_color_DS.insert(0, default_color_DS)  # Standard value for color of the curve
    entry_color_DS.pack(side=tk.LEFT, padx=5)


    # Input field for the color of the Newton polygon
    frame_color_NP = tk.Frame(scrollable_frame)
    frame_color_NP.pack(pady=10)

    tk.Label(frame_color_NP, text="Color of Newton polygon:").pack(side=tk.LEFT)
    entry_color_NP = tk.Entry(frame_color_NP, width=10)
    entry_color_NP.insert(0, default_color_NP)  # Standard value for color of the curve
    entry_color_NP.pack(side=tk.LEFT, padx=5)

    # Checkbox for make_grid
    frame_make_grid = tk.Frame(scrollable_frame)
    frame_make_grid.pack(pady=10)

    make_grid_var = tk.BooleanVar(value=default_make_grid)
    make_grid_checkbox = tk.Checkbutton(
        frame_make_grid, text="Draw a grid", variable=make_grid_var
    )
    make_grid_checkbox.pack(side=tk.LEFT)

    # Checkbox for print_vertices
    frame_print_vertices = tk.Frame(scrollable_frame)
    frame_print_vertices.pack(pady=10)

    print_vertices_var = tk.BooleanVar(value=default_print_vertices)
    print_vertices_checkbox = tk.Checkbutton(
        frame_print_vertices, text="Print ends of edges", variable=print_vertices_var
    )
    print_vertices_checkbox.pack(side=tk.LEFT)

    # Separator
    separator = ttk.Separator(scrollable_frame, orient="horizontal")
    separator.pack(fill="x", padx=10, pady=5)

    # Options for weights
    instruction = tk.Label(
            scrollable_frame,
            text=(
                "Options for weights"
            ),
            font=("Arial", 12),
            justify="center",
            wraplength=500
        )

    instruction.pack(pady=5)

    # Checkbox for show_weights
    frame_show_weights = tk.Frame(scrollable_frame)
    frame_show_weights.pack(pady=10)

    show_weights_var = tk.BooleanVar(value=True)
    show_weights_checkbox = tk.Checkbutton(
        frame_show_weights, text="Show weights of edges", variable=show_weights_var
    )
    show_weights_checkbox.pack(side=tk.LEFT)

    # Input field for the size of the weights
    frame_weight_size = tk.Frame(scrollable_frame)
    frame_weight_size.pack(pady=10)

    tk.Label(frame_weight_size, text="Size of weights:").pack(side=tk.LEFT)
    entry_weight_size = tk.Entry(frame_weight_size, width=10)
    entry_weight_size.insert(0, default_size_weights)  # Standard value for weight_size
    entry_weight_size.pack(side=tk.LEFT, padx=5)

    # Input field for the color of the weights
    frame_weight_color = tk.Frame(scrollable_frame)
    frame_weight_color.pack(pady=10)

    tk.Label(frame_weight_color, text="Color of weights:").pack(side=tk.LEFT)
    entry_weight_color = tk.Entry(frame_weight_color, width=10)
    entry_weight_color.insert(0, default_color_weights)  # Standard value for weight_color
    entry_weight_color.pack(side=tk.LEFT, padx=5)
    
    # Separator
    separator = ttk.Separator(scrollable_frame, orient="horizontal")
    separator.pack(fill="x", padx=10, pady=5)

    # Options for tikz
    instruction = tk.Label(
            scrollable_frame,
            text=(
                "Options for tikz code"
            ),
            font=("Arial", 12),
            justify="center",
            wraplength=500
        )

    instruction.pack(pady=5)

    # Checkbox for print_tikz_code_curve
    frame_print_tikz_code_curve = tk.Frame(scrollable_frame)
    frame_print_tikz_code_curve.pack(pady=10)

    print_tikz_code_curve_var = tk.BooleanVar(value=default_print_tikz_curve)
    print_tikz_code_curve_checkbox = tk.Checkbutton(
        frame_print_tikz_code_curve, text="Print tikz code curve", variable=print_tikz_code_curve_var
    )
    print_tikz_code_curve_checkbox.pack(side=tk.LEFT)

    # Checkbox for print_tikz_code_dual
    frame_print_tikz_code_dual = tk.Frame(scrollable_frame)
    frame_print_tikz_code_dual.pack(pady=10)

    print_tikz_code_dual_var = tk.BooleanVar(value=default_print_tikz_DS)
    print_tikz_code_dual_checkbox = tk.Checkbutton(
        frame_print_tikz_code_dual, text="Print tikz code dual", variable=print_tikz_code_dual_var
    )
    print_tikz_code_dual_checkbox.pack(side=tk.LEFT)

    # Separator
    separator = ttk.Separator(scrollable_frame, orient="horizontal")
    separator.pack(fill="x", padx=10, pady=5)

    instruction = tk.Label(
            scrollable_frame,
            text=(
                "Version 2.3.0 \n Author: Nathan Tiggemann"
            ),
            font=("Arial", 8),
            justify="center",  # Center text
            wraplength=500     # Max width for new line (in pixels)
        )

    instruction.pack(pady=0)  # Add distances to boundary

    def open_link():
        webbrowser.open("https://ntiggemann.github.io/coding.html")  # opens link in browser

    link_label = tk.Label(scrollable_frame, text="More information at https://ntiggemann.github.io/coding.html", font=("Arial", 8), justify="center", fg="blue", cursor="hand2")
    link_label.pack(pady=0)

    # Add open link
    link_label.bind("<Button-1>", lambda e: open_link())

    # Start GUI loop
    window.mainloop()
if __name__ == "__main__":
    main()