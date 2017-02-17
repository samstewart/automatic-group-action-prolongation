from itertools import combinations_with_replacement, chain, combinations, imap
import pdb
import functools
import operator
from sympy import *



# all derivatives up to and including max_order
def all_derivatives(symbol, indep_vars, max_order):
	# concatenate the list of each order of derivative
	return reduce(chain, 
		map(lambda k: all_derivatives_of_order_k(symbol, indep_vars, k), 
			range(max_order + 1)))

# generating list of all derivatives at to order k
def all_derivatives_of_order_k(symbol, indep_vars, k):
	return imap(lambda v: symbol + ''.join(v), combinations_with_replacement(indep_vars, k))

def generate_vector_field_coefficients(indep_vars, max_derivative_order, base_letter):
	# we add +1 because we need a coefficient for the independent variable
	variable_names = [ base_letter + str(i) for i in range(len(indep_vars) + 1) ]

	# now we generate all possible derivatives of each variable
	all_ds = reduce(chain,
		map(lambda d: all_derivatives(d, indep_vars, max_derivative_order), variable_names))

	return list(all_ds)

indep_vars = 'xt'
vector_field_coeffs = generate_vector_field_coefficients('uxt', 2, 'A') + generate_vector_field_coefficients('uxt', 2, 'B')
domain = 'ZZ[%s]' % (','.join(vector_field_coeffs), )
vector_derivatives = symbols(vector_field_coeffs)

u_derivatives = symbols(list(all_derivatives('u', indep_vars, 2)))

def differentiate_var(u, var):
	"""Derivative of function with direct dependence on var (don't need a total derivative)"""
	name = str(u)

	new_symbol_name = (name + var)

	# make sure the derivatives are in a canonical order
	new_symbol_name = name[0] + ''.join(sorted(new_symbol_name[1:]))

	# should already be in the symbols bank
	return sympify(new_symbol_name)

def dxi_dj(xi, var):
	"""Differentiates the coefficient of the generator vector that depends on both 
	the independent and dependent variables.
	The formula is really just the chain rule

	D(A, t) = A_t + A_u u_t
	"""

	# we compute the A_t part
	At = differentiate_var(xi, var)

	# compute the A_u
	ut = differentiate_var('u', var)

	# compute the u_t part
	Au = differentiate_var(xi, 'u')

	return At + ut * Au

# differentiates a single factor depending on its type
# constant -> 0
# vector coefficient dxi_dj
# independent function -> differentiate_var
def diff_factor(factor, var):
	if factor.is_integer:
		return 0

	elif factor in u_derivatives:
		# independent variable
		return differentiate_var(factor, var)

	elif factor in vector_derivatives:
		# vector coefficients that depend on the independent solution as well
		return dxi_dj(factor, var)
	else:
		# do nothing if we don't know what it is
		return factor 

def product_rule(monomial, f):
	"""Higher order function that applies a Leibniz rule and applies
	the function to each term (like a derivative)"""
	
	factors = monomial.args

	# edge case of just one factor
	if len(factors) == 0:
		return f(monomial)

	# generate the indices for the nondifferentiated terms
	non_diff = list(combinations(factors, len(factors) - 1))

	# differentiate the other terms
	differentiated_factors = map(f, factors)

	# reverse the list to agree with the ordering of the nondifferentiated elements
	differentiated_factors.reverse()

	# now recombine the lists by multiplying out the factors
	return sum(map(lambda a: operator.mul(reduce(operator.mul, a[1], 1), a[0]), zip(differentiated_factors, non_diff)))

def diff_term(monomial, var):
	return product_rule(monomial, lambda term: diff_factor(term, var))

def addition_rule(polynomial, f):
	"""Higher order function that applies a given function to each term in a sum"""

	if len(polynomial.args[0].args) == 0:
		return f(polynomial.args[0])

	return sum(map(f, polynomial.args[0].args))

def diff_entire_poly(p, var):
	return Poly(addition_rule(p, lambda term: diff_term(term, var)), domain=domain)

# compute multivariable derivative (i.e. D^j p(x))
def multi_variable_derivative(p, derivatives):
	if len(derivatives) == 0:
		return p

	if len(derivatives) == 1:
		return diff_entire_poly(p, derivatives[0])

	return diff_entire_poly(multi_variable_derivative(p, derivatives[:-1]), derivatives[-1])

# compute the coefficient of a prolongation (so one coordinate in the jet bundle).
# Symbolically, the prolongation formula is given by (with one dependent variable)
# \phi^J = D^J(\phi - \sum_{i = 1}^p \xi_i u_i) + \sum_{i = 1}^p \xi_i u_{J, i}.
def prolongation(indep_coeff, derivatives, indep_vars):
	# the first sum
	first_sum = map(lambda indexAndVariable: sympify('A%d * u%s' % (indexAndVariable[0], indexAndVariable[1])), zip(range(len(indep_vars)), indep_vars))
	first_sum = Poly(sympify(indep_coeff) - sum(first_sum), domain=domain)
	print first_sum
	# now we compute the total derivative
	first_sum = multi_variable_derivative(first_sum, derivatives)

	# the second sum
	second_sum = map(lambda indexAndVariable: sympify('A%d * %s' % (indexAndVariable[0], differentiate_var('u', derivatives + indexAndVariable[1]))), zip(range(len(indep_vars)), indep_vars))
	second_sum = sum(second_sum)


	return first_sum


pdb.set_trace()
print prolongation('B0', 'x', 'xt')





