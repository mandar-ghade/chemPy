from .compound import Compound
from .compound_counter import CompoundCounter
import sympy as smp
from string import ascii_lowercase
from fractions import Fraction
from typing import Any, Self, Optional
from .utils import split_str

class Meta(type):
    def __call__(cls, *args: Any, **kwargs: Any) -> Any:
        if len(args) == 0:
            return 'No arguments specified.'
        reactants = None
        products = None
        for arg in args:
            if type(arg) is str:
                if reactants and products:
                    break
                reactants, products = arg.split(split_str(arg))
                reactants = CompoundCounter([Compound(reactant) 
                            for reactant in reactants.split('+')])
                products = CompoundCounter([Compound(product) 
                            for product in products.split('+')])
                break
            elif type(arg) is CompoundCounter:
                if reactants and products:
                    break
                elif not reactants and not products:
                    reactants = arg
                elif reactants and not products:
                    products = arg
            else:
                raise TypeError('Expected `CompoundCounter` or `str` for `reactants` or `products` argument but received '
                                f'`{arg.__class__.__name__}` instead.')
        cls.reactants = reactants
        cls.products = products
        return super().__call__(*args, **kwargs)


class Equation2(metaclass=Meta):
    def __new__(cls, *args, **kwargs) -> Self:
        instance = super().__new__(cls)
        return instance
    
    def __init__(self, *args, **kwargs) -> None:
        self.reactants: CompoundCounter = self.reactants
        self.products: CompoundCounter = self.products
        self.compounds: list[Compound] = self.reactants.flatten() + self.products.flatten()
        self.equation = self._equation()
        self.h_rxn: Optional[int | float] = None

    def _center_str(self, inp: int | float, latex: bool) -> str:
        """Encases the string in \\text{} for use when displayed in latex. Typically used to format coefficients."""
        return f'{inp:,}' if latex in [None, False] else f'\\text{{{inp:,}}}'
    
    def _latex(self, inp: Compound, latex: bool) -> str:
        """Returns a latexified or normal representation of the equation, denoted by whether latex is specified to be True or False."""
        return inp.comp_str.strip() if latex in [None, False] else inp.latexify()  

    def _equation(self, latex: bool = False) -> str:
        return ' + '.join(self._center_str(reactant.coefficient, latex) + \
                          f'({self._latex(reactant, latex)})'
                          if float(str(reactant.coefficient).replace(',', '')) != 1
                          else self._latex(reactant, latex)
                          for reactant in self.reactants.keys()) + ' → ' + \
               ' + '.join(self._center_str(product.coefficient, latex) + \
                          f'({self._latex(product, latex)})'
                          if float(str(product.coefficient).replace(',', '')) != 1
                          else self._latex(product, latex)
                          for product in self.products.keys())
    
    def _get_coefficients(self, data: list[list[int]]) -> None:
        matrix = smp.Matrix(data)
        length = matrix.shape[1]
        symbols = [smp.Symbol(f'{ascii_lowercase[n]}') for n in range(length)]
        solutions = smp.linsolve(matrix, symbols)
        solutions = solutions.subs([(symbol, 1) for symbol in symbols])
        if solutions.args == ():
            print('Impossible equation')
            return None
        args = solutions.args[0]
        args: list[Fraction] = [Fraction(abs(arg)).limit_denominator() for arg in args]
        denom_list = [frac.denominator for frac in args if frac.denominator != 1]
        if denom_list != []:
            while not all(arg.denominator == 1 for arg in args):
                denom_list = [frac.denominator for frac in args if frac.denominator != 1]
                args = [arg * max(denom_list) for arg in args]
        coefficients = [arg.numerator for arg in args]
        for i, compound in enumerate(self.reactants.flatten() + self.products.flatten()):
            matching_coefficient = coefficients[i]
            if i <= len(self.reactants) - 1:
                self.reactants[compound] = matching_coefficient
            else:
                self.products[compound] = matching_coefficient
            compound.coefficient = matching_coefficient
        self.compounds = self.reactants.flatten() + self.products.flatten()

    def balance(self) -> None:
        """
        Balances the Chemical Equation.
        """
        if self.reactants.is_equal(self.products):
            self.equation = self._equation()
            return
        matrix = []
        for element in self.reactants.elements:
            row: list[int] = []
            row.extend(self.reactants.search_element_occurances(element).values())
            row.extend(-val for val in self.products.search_element_occurances(element).values())
            matrix.append(row)
        self._get_coefficients(matrix)
        if not self.reactants.is_equal(self.products):
            print('An error occured while balancing the chemical equation.')
        self.equation = self._equation()

    def _repr_latex_(self) -> str:
        """Returns latex representation of an equation."""
        return f'${self._equation(latex=True)}$'

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reactants}, {self.products})"
    
    def __str__(self) -> str:
        """Returns equation."""
        return self.equation