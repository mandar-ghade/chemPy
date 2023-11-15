from .compound import Compound
from .data import H_rxn
from .element import Element
from .utils import split_str
from collections import Counter
from fractions import Fraction
import numpy as np
from string import ascii_lowercase
import sympy as smp
from typing import Optional, Self


def solution(U):
    # find the eigenvalues and eigenvector of U(transpose).U
    e_vals, e_vecs = np.linalg.eig(np.dot(U.T, U))  
    # extract the eigenvector (column) associated with the minimum eigenvalue
    return e_vecs[:, np.argmin(e_vals)] 


class Equation:
    def __init__(self, reactants: list[Compound], products: list[Compound]):
        if not isinstance(reactants, list):
            raise TypeError('Expected `list[Compound]` for `reactants` argument but received '
                            f'`{reactants.__class__.__name__}` instead.')
        if not isinstance(products, list):
            raise TypeError('Expected `list[Compound]` for `products` argument but received '
                            f'`{products.__class__.__name__}` instead.')
        self.reactants = reactants
        self.products = products
        self.coefficients = [1 for _ in range(len(self.reactants) + len(self.products))]
        self.equation = self._equation()
        self.is_balanced = False

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}([{self.reactants}], [{self.products}])"

    def _center_str(self, inp: int, latex: bool) -> str:
        """Encases the string in \\text{} for use when displayed in latex. Typically used to format coefficients."""
        return str(inp) if latex in [None, False] else f'\\text{{{inp}}}'

    def _latex(self, inp: Compound, latex: bool) -> str:
        """Returns a latexified or normal representation of the equation, denoted by whether latex is specified to be True or False."""
        return inp.comp_str.strip() if latex in [None, False] else inp.latexify()  

    def _equation(self, latex: bool = False) -> str:
        """Returns the equation as a string in either latex or non-latex format."""
        return ' + '.join(
            [
                self._center_str(coef, latex) + f'({self._latex(reactant, latex)})' 
                if int(str(coef).replace(',', '')) > 1 
                else self._latex(reactant, latex)
                for reactant, coef in 
                zip(self.reactants, self.coefficients[:len(self.reactants)])
            ]) \
        + ' → ' + \
            ' + '.join(
                [
                    self._center_str(coef, latex) + f'({self._latex(product, latex)})' 
                    if int(str(coef).replace(',', '')) > 1 
                    else self._latex(product, latex)
                    for product, coef in zip(self.products, self.coefficients[len(self.reactants):])
                ])

    def _repr_latex_(self) -> str:
        """Returns latex representation of an equation."""
        return f'${self._equation(latex=True)}$'
    
    def __str__(self) -> str:
        """Returns equation."""
        return self.equation
    
    def total_left(self) -> Counter[Element]:
        """Returns a Counter of elements to the left of the equation (in the reactants)."""
        reactants = Counter()
        for reactant in self.reactants:
            compound = reactant.elements.most_common()
            for (element, count) in compound:
                reactants[element] += count 
        return reactants

    def total_right(self) -> Counter[Element]:
        """Returns a Counter of elements to the right of the equation (in the products)."""
        products = Counter()
        for product in self.products:
            compound = product.elements.most_common()
            for (element, count) in compound:
                products[element] += count 
        return products

    def count_left(self, element: Element) -> list[int]:
        """Counts occurances of an element out of all the reactants."""
        return [compound.count(element) for compound in self.reactants]

    def count_right(self, element: Element) -> list[int]:
        """Counts occurances of an element out of all the products."""
        return [compound.count(element) for compound in self.products]

    def _get_coefficients(self, data: list[list[int]]) -> list[str]:
        """
        Returns a list of coefficients, 
        whose indices correspond to the positions of compounds in the given equation.
        """
        matrix = smp.Matrix(data)
        length = matrix.shape[1]
        symbols = [smp.Symbol(f'{ascii_lowercase[n]}') for n in range(length)]
        solutions = smp.linsolve(matrix, symbols)
        solutions = solutions.subs([(symbol, 1) for symbol in symbols])
        if solutions.args == ():
            print('Impossible equation')
            return [f'{arg:,}' for arg in self.coefficients]
        args = solutions.args[0]
        args: list[Fraction] = [Fraction(abs(arg)).limit_denominator() for arg in args]
        denom_list = [frac.denominator for frac in args if frac.denominator != 1]
        if denom_list != []:
            while not all(arg.denominator == 1 for arg in args):
                denom_list = [frac.denominator for frac in args if frac.denominator != 1]
                args = [arg * max(denom_list) for arg in args]
        return [f'{arg.numerator:,}' for arg in args]

    def balance(self) -> None:
        """Balances the chemical equation."""
        reactant_elements = self.total_left()
        product_elements = self.total_right()
        if reactant_elements == product_elements:
            print('Equation already balanced')
        matrix = []
        for element in reactant_elements:
            row = [0 for _ in range(len(self.reactants) + len(self.products))]
            left_count = self.count_left(element)
            right_count = self.count_right(element)
            for count, i in zip(left_count, range(len(self.reactants))):
                row[i] = count
            for count, i in zip(right_count, range(len(self.reactants), len(self.products)+len(self.reactants)+1)):
                row[i] = -count
            matrix.append(row)
        self.coefficients = self._get_coefficients(matrix)
        self.equation = self._equation()
        self.is_balanced = True

    def extend(self, equation: Self, 
               h_rxn: H_rxn = None, 
               desired_eq: Optional[tuple[Self, H_rxn]] = None
              ) -> None:
        """
        Performs Hess's law problems.
    
        ## Example:
        ```
        >>> eq = Equation.parse_from_string('H2 + O2 -> H2O')
        >>> assert eq.H_rxn == 0  # (kJ/mol)
        >>> eq.extend(Equation.parse_from_string('H2O -> H2 + O'), -234)
        ```
        (Equation.parse_from_string('H2O -> H2 + O'), -234) represents (Equation, H_rxn in kJ/mol)
        """
        # make input a list of tuples
        assert isinstance(desired_eq, tuple)
        assert isinstance(equation, self.__class__)
        assert isinstance(H_rxn, Optional[float | int])
        assert self.is_balanced
        equation.balance()
        intermediates = set(self.products).intersection(equation.reactants)



    @classmethod
    def parse_from_string(cls, line: str):
        reactants, products = line.split(split_str(line))
        reactants = [Compound(reactant) 
                     for reactant in reactants.split('+')]
        products = [Compound(product) 
                    for product in products.split('+')]
        return cls(reactants, products)
