
from collections import MutableMapping

import sympy

from ode_system import ODESystem
from sympy_helper import unique_array_stable

class ChemicalSpecies(object):
    ''' Chemical species, A_i from Harrington paper '''
    def __init__(self, id_):
        self.id_ = id_

    def _key(self):
        ''' Return the key of the object for comparison '''
        return str(self.id_)

    def __eq__(self, other):
        if isinstance(other, ChemicalSpecies):
            return self._key() == other._key()
        return False

    def __hash__(self):
        return hash(self._key())

    def __repr__(self):
        return str(self.id_)

class Complex(MutableMapping):
    ''' A complex of ChemicalSpecies '''
    def __init__(self, *args, **kwargs):
        self._species_dict = dict(*args, **kwargs)
        for key in self._species_dict:
            if not isinstance(key, ChemicalSpecies):
                raise ValueError('{} must be a ChemicalSpecies'.format(key))

    def __iter__(self):
        return self._species_dict.__iter__()

    def __len__(self):
        return self._species_dict.__len__()

    def __getitem__(self, item):
        return self._species_dict.__getitem__(item)

    def __setitem__(self, key, value):
        if not isinstance(key, ChemicalSpecies):
            raise ValueError('{} must be a ChemicalSpecies'.format(key))
        self._species_dict[key] = value

    def __repr__(self):
        return ' + '.join(['{}.{}'.format(val, key) for key, val in self.iteritems()])

    def __delitem__(self, key):
        self._species_dict.__delitem__(key)

    def as_vector(self, variable_order):
        return tuple([self._species_dict.get(variable, 0) for variable in variable_order])

class Reaction(object):
    ''' Represents a reaction between complexes '''
    def __init__(self, complex1, complex2):
        self.complex1 = complex1
        self.complex2 = complex2

    def __repr__(self):
        return '{} -> {}'.format(self.complex1, self.complex2)

class ChemicalReactionNetwork(object):
    ''' Chemical reaction network '''
    def __init__(self, chemical_species, complexes, reactions):
        assert all([isinstance(c_s, ChemicalSpecies) for c_s in chemical_species])

        self.chemical_species = tuple(chemical_species)
        self.complexes = tuple(complexes)
        self.reactions = tuple(reactions)

        # Check the complexes only involve our species
        set_c_s = set(self.chemical_species)
        for complex in self.complexes:
            assert set(complex.keys()).issubset(set_c_s)

        # Check our reactions only involve our complexes
        for reaction in self.reactions:
            assert reaction.complex1 in self.complexes
            assert reaction.complex2 in self.complexes

    @property
    def p(self):
        return len(self.complexes)

    @property
    def n(self):
        return len(self.chemical_species)

    @property
    def r(self):
        return len(self.reactions)

    def ode_equations(self):
        ''' Return a tuple of differential equations for each species '''
        sympy_chem_spec = map(lambda x: sympy.var(str(x)), self.chemical_species)

        rate_reaction_function = []
        stochiometric_matrix = []

        for reaction in self.reactions:
            # Calculate the rate reaction function
            # Create a rate parameter
            i, j = self.complexes.index(reaction.complex1), self.complexes.index(reaction.complex2)
            rate_constant = sympy.var('k_{}_{}'.format(i, j))

            concentration = [spec ** exponent for spec, exponent in
                             zip(sympy_chem_spec, reaction.complex1.as_vector(self.chemical_species))]

            rate_reaction_function.append(rate_constant * sympy.prod(concentration))

            # Calculate the stochiometric matrix
            col = (sympy.Matrix(reaction.complex2.as_vector(self.chemical_species)) -
                   sympy.Matrix(reaction.complex1.as_vector(self.chemical_species)))
            stochiometric_matrix.append(col.T)

        rate_reaction_function = sympy.Matrix(rate_reaction_function)
        assert rate_reaction_function.shape == (self.r, 1)
        stochiometric_matrix = sympy.Matrix(stochiometric_matrix).T
        assert stochiometric_matrix.shape == (self.n, self.r)

        equations = stochiometric_matrix * rate_reaction_function
        assert len(sympy_chem_spec) == len(equations)
        return tuple(equations)

    def to_ode_system(self):
        ''' Generate a system of ODEs based on the current network '''
        sympy_chem_spec = map(lambda x: sympy.var(str(x)), self.chemical_species)
        equations = self.ode_equations()
        deriv_dict = dict(zip(sympy_chem_spec, equations))
        return ODESystem.from_dict(deriv_dict=deriv_dict)

    def __repr__(self):
        return '\n'.join([reaction.__repr__() for reaction in self.reactions])

    @classmethod
    def from_diagram(cls, diagram):
        ''' Given a text diagram, return an interpreted chemical reaction network '''
        species = []
        complexes = []
        reactions = []
        complex0 = Complex({})
        for reaction in diagram.strip().split('\n'):
            _complexes = reaction.strip().split('->')
            if len(_complexes) != 2:
                raise ValueError('Invalid reaction: {}'.format(reaction))

            # Process the left hand side
            complex_left = _complexes[0].strip()
            if complex_left:
                complex_left = sympy.sympify(complex_left)
                complex_left = {ChemicalSpecies(k): v for k, v in complex_left.as_coefficients_dict().iteritems()}
                species.extend(complex_left.keys())
                complex_left = Complex(complex_left)
            else:
                complex_left = complex0

            complex_right = _complexes[1].strip()
            if complex_right:
                complex_right = sympy.sympify(complex_right)
                complex_right = {ChemicalSpecies(k): v for k, v in complex_right.as_coefficients_dict().iteritems()}
                species.extend(complex_right.keys())
                complex_right = Complex(complex_right)
            else:
                complex_right = complex0

            if complex_left not in complexes:
                complexes.append(complex_left)
            if complex_right not in complexes:
                complexes.append(complex_right)

            reactions.append(Reaction(complex_left, complex_right))

        species = unique_array_stable(species)
        return cls(species, complexes, reactions)



if __name__ == '__main__':
    import doctest
    doctest.testmod()