###############################
# This script prints the product of bivariate tropical polynomials, so that it can be
# used in my program to plot tropical curves
###############################

''' Preamble '''
# Don't change anything
import numpy as np
import re as re
import itertools
# Auxiliar-minus-infty - just some very negative number
infty = -9999999999999999 #DO NOT CHANGE, even if you want to use the min convention!

''' Define the polynomials '''
# Should we use min? (False or True)
useMin = False

# The two polynomials you want to multiply. The input works just like in my other program
poly_a = "-4x3+-1xy+x-2y2+y-2"
poly_b = "-1x+1y+0"

''' Functions '''
# Dont change stuff after here
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

    # New dictionary with sorted keys
    sorted_coefficients = dict(sorted(coefficients.items(), key=lambda x: x[0], reverse=True))

    return sorted_coefficients

coeff_a = parse_polynomial_trop(poly_a,useMin)
coeff_b = parse_polynomial_trop(poly_b,useMin)

coeff = {}

# Multiply the polynomials
for exp_a, exp_b in itertools.product(coeff_a.keys(), coeff_b.keys()):
    exp = (exp_a[0] + exp_b[0],exp_a[1] + exp_b[1])
    if exp in coeff:
        if useMin:
            coeff[exp] = min(coeff[exp],coeff_a[exp_a] + coeff_b[exp_b])
        else:
            coeff[exp] = max(coeff[exp],coeff_a[exp_a] + coeff_b[exp_b])
    else:
        coeff[exp] = coeff_a[exp_a] + coeff_b[exp_b]

# The future output
poly_c = ""
# Write the output
for exp in coeff:
    poly_c = poly_c + "+" + str(coeff[exp]) + "x" + str(exp[0]) + "y" + str(exp[1])

# Print output
print(poly_c)