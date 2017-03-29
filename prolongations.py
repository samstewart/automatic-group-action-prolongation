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
	return imap(lambda v: symbol + ''.join(sorted(v)), combinations_with_replacement(indep_vars, k))

def generate_vector_field_coefficients(indep_vars, dep_vars, max_derivative_order, base_letter):
	variable_names = [ base_letter + str(i) for i in range(len(indep_vars)) ]

	# now we generate all possible derivatives of each variable
	all_ds = reduce(chain,
		map(lambda d: all_derivatives(d, dep_vars + indep_vars, max_derivative_order), variable_names))

	return list(all_ds)

indep_vars = 'xt'
vector_field_coeffs = generate_vector_field_coefficients('xt', 'u', 2, 'A') + list(all_derivatives('B0', 'uxt', 2))
domain = 'ZZ[%s]' % (','.join(vector_field_coeffs), ) # adjoin them as elements to the polynomial ring
vector_derivatives = symbols(vector_field_coeffs)

u_derivatives = symbols(list(all_derivatives('u', indep_vars, 3)))


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
		return 0 # differentiating a constant gives zero

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
	"""Higher order function that applies a given function to each monomial in a polynomial"""
	pdb.set_trace()

	if degree(polynomial) == 0:
		return f(polynomial.args[0])

	return sum(map(f, polynomial.args[0].args))

def diff_entire_poly(p, var):
	# note: we must expand the polynomial into sum of monomials
	return Poly(addition_rule(expand(p), lambda term: diff_term(term, var)), *u_derivatives, domain=domain)

# compute multivariable derivative (i.e. D^j p(x))
def multi_variable_derivative(poly, derivatives):
	pdb.set_trace()

	if len(derivatives) == 0:
		return poly

	if len(derivatives) == 1:
		return diff_entire_poly(poly, derivatives[0])

	return diff_entire_poly(multi_variable_derivative(poly, derivatives[:-1]), derivatives[-1])

# compute the coefficient of a prolongation (so one coordinate in the jet bundle).
# Symbolically, the prolongation formula is given by (with one dependent variable)
# \phi^J = D^J(\phi - \sum_{i = 1}^p \xi_i u_i) + \sum_{i = 1}^p \xi_i u_{J, i}.
def prolongation(indep_coeff, derivatives, indep_vars):
	# the first sum
	first_sum = map(lambda indexAndVariable: sympify('A%d * u%s' % (indexAndVariable[0], indexAndVariable[1])), zip(range(len(indep_vars)), indep_vars))

	first_sum = Poly(sympify(indep_coeff) - sum(first_sum), domain=domain)

	# now we compute the total derivative
	first_sum = multi_variable_derivative(first_sum, derivatives)

	# the second sum
	second_sum = map(lambda indexAndVariable: sympify('A%d * %s' % (indexAndVariable[0], differentiate_var('u', derivatives + indexAndVariable[1]))), zip(range(len(indep_vars)), indep_vars))
	second_sum = Poly(sum(second_sum), domain=domain)

	return first_sum

# TODO: mod out by relation u_t = u_xx for example (the original PDE needs to be baked into the computation).
# TODO: extract the coefficients from the polynomial.
# pdb.set_trace()
# find \phi_{xx}
phixx = Poly(
	"""B0xx + (2 * B0ux - A0xx)*ux - A1xx * ut 
	+ (B0uu - 2 * A0ux) * ux**2 - 2*A1ux * ux * ut 
	- A0uu * ux**3 - A1uu * ux**2 * ut 
	+ (B0u - 2 * A0x) * uxx - 2 * A1x * utx 
	- 3 * A0u * ux * uxx - A1u * ut * uxx 
	- 2 * A1u * ux * utx""", *u_derivatives, domain=domain);

phix = Poly('B0x + (B0u - A0x) * ux - A1x * ut - A0u * ux**2 - A1u * ux * ut', *u_derivatives, domain=domain);

#print multi_variable_derivative(Poly('-3*A0', sympify('ux'), domain=domain), 'xx').args[0] + sympify('3*A0u*uxx')


# testing down to the level of total derivatives works.
# p1 = multi_variable_derivative(Poly('B0', *u_derivatives, domain=domain), 'xx')
# p2 = Poly('-ux', *u_derivatives, domain=domain) * multi_variable_derivative(Poly('A0', *u_derivatives, domain=domain), 'xx')
# p3 = Poly('-ut', *u_derivatives, domain=domain) * multi_variable_derivative(Poly('A1', *u_derivatives, domain=domain), 'xx')
# p4 = Poly('-2*uxx', *u_derivatives, domain=domain) * multi_variable_derivative(Poly('A0', *u_derivatives, domain=domain), 'x')
# p5 = Poly('-2*utx', *u_derivatives, domain=domain) * multi_variable_derivative(Poly('A1', *u_derivatives, domain=domain), 'x')

# test to see if the same answer
#print ((p1 + p2 + p3 + p4 + p5) - phixx).args[0]

# now test prolongation formula
# p = multi_variable_derivative(Poly('B0 - A0 * ux - A1 * ut', *u_derivatives, domain=domain), 'xx') + Poly('A0 * uxxx + A1 * utxx', *u_derivatives, domain=domain)

print multi_variable_derivative(Poly('B0 - A0', *u_derivatives, domain=domain), 'x').args[0]

