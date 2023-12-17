from .utils import tokenize, strip_coefficients, extract_coefficient
from collections import Counter
from .element import Element
from .subscript import Subscript
from typing import Optional, Self

class Compound:
    def __init__(self, comp_str: str, tokens: Optional[list[Element]] = None,
                 subscripts: Optional[list[Subscript]] = None,
                 mass: Optional[float] = None):
        if not isinstance(comp_str, str):
            raise TypeError('Expected `str` for `comp_str` argument but received '
                            f'`{comp_str.__class__.__name__}` instead.')
        size, no_subs = extract_coefficient(comp_str.strip())
        comp_str = strip_coefficients(comp_str)
        if tokens is None:
            tokens, subscripts = tokenize(comp_str)
        self.tokens = tokens
        self.elements = Counter(tokens)
        self.comp_str = comp_str
        self.subscripts: Optional[list[Subscript]] = subscripts
        self.molar_mass = self._get_molar_mass()
        self.valence_electrons = self._get_valence_electrons()
        self.electrons = self._get_total_electrons()
        self.coefficient: int | float = 1
        if not no_subs:
            if size.is_integer():
                self.coefficient = int(size)
            else:
                self.coefficient = size
        self.mass = mass
        self.moles: Optional[float] = None
        if self.mass is not None:
            self.moles = self.mass / self.molar_mass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.comp_str}')"

    def __add__(self, other: Self):
        if isinstance(other, Element):
            return Compound(self.comp_str + other.symbol)
        elif isinstance(other, self.__class__):
            if self.comp_str == other.comp_str:
                coef = self.coefficient + other.coefficient
                return Compound(f'{coef}{self.comp_str}')
            return Compound(self.comp_str + other.comp_str)
        
    def __mul__(self, other: int):
        assert isinstance(other, int)
        return Compound(f'[{self.coefficient if self.coefficient != 1 else ''}{self.comp_str}]{other}')

    def _get_total_electrons(self) -> int:
        """Returns the sum of the total electrons in the compound per number of elements. (Assigned to self.electrons)"""
        return sum(element.electrons * count 
                   for element, count in self.elements.items())
    
    def _get_molar_mass(self) -> float:
        """Returns the total molar mass of all the elements in the Compound. (Assigned to self.molar_mass)"""
        return sum(element.molar_mass  * count 
                   for element, count in self.elements.items())
    
    def _get_valence_electrons(self) -> int:
        """Returns the sum of the valence electrons in the Compound (calculated using the QM model)."""
        return sum(element.valence_electrons * count 
                   for element, count in self.elements.items())

    def count(self, element) -> int:
        """Returns the number of occurances an Element has in a compound."""
        return sum(count for e, count in Counter(self.tokens).items() if e == element)

    def subscript(self, index: int) -> Optional[int]:
        """Returns the count, or subscript of an element at a particular index in the compound string."""
        if self.subscripts is None:
            return None
        subs: list[int] = [
            sub.subs for sub in self.subscripts
            if sub.start_index - 1 == index and sub.subs is not None]
        if len(subs) == 0:
            return None
        elif len(subs) > 1:
            raise Exception(f'Multiple subscripts show under comp_str={self.comp_str}, index={index}')
        return int(subs[0])

    def latexify(self) -> str:
        """Returns a latex representation of the compound string."""
        if self.subscripts == []:
            return self.comp_str
        latex_str = ''
        for i, c in enumerate(self.comp_str):
            if c.isdigit() or c == ' ' or c == '.':
                continue
            subs = self.subscript(i)
            prefix = '\\text{'
            if subs is None: # check delims to eliminate \text redundancy.
                latex_str += prefix + c + '}'
            else:
                latex_str += prefix + c + '}' + '_{' + str(subs) + '}'
        return latex_str

    def __hash__(self) -> int:
        return hash(tuple(item for item in self.elements.items()))

    def __eq__(self, other: Self) -> bool:
        return self.elements == other.elements
    