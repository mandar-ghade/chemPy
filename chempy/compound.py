from .utils import tokenize
from collections import Counter
from .element import Element
from .subscript import Subscript
from typing import Optional, Self

class Compound:
    def __init__(self, comp_str: str, tokens: Optional[list[Element]] = None,
                 subscripts: Optional[list[Subscript]] = None):
        if not isinstance(comp_str, str):
            raise TypeError('Expected `str` for `comp_str` argument but received '
                            f'`{comp_str.__class__.__name__}` instead.')
        if tokens is None:
            tokens, subscripts = tokenize(comp_str)
        self.tokens = tokens
        self.elements = Counter(tokens)
        self.comp_str = comp_str
        self.subscripts: list[Subscript] = subscripts
        self.molar_mass = self._get_molar_mass()
        self.valence_electrons = self._get_valence_electrons()
        self.electrons = self._get_total_electrons()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.comp_str}', {self.tokens}, {self.subscripts})"

    def __add__(self, other: Self):
        if isinstance(other, Element):
            return Compound(self.comp_str + other.symbol)
        elif isinstance(other, self.__class__):
            return Compound(self.comp_str + other.comp_str)
        
    def __mul__(self, other: int):
        assert isinstance(other, int)
        return Compound(f'[{self.comp_str}]{other}')

    def _get_total_electrons(self) -> int:
        return sum(token.electrons * count 
                   for token, count in self.elements.items())
    
    def _get_molar_mass(self) -> float:
        return sum(token.molar_mass  * count 
                   for token, count in self.elements.items())
    
    def _get_valence_electrons(self) -> int:
        return sum(token.valence_electrons * count 
                   for token, count in self.elements.items())

    def count(self, element) -> int:
        count = [count for e, count in Counter(self.tokens).items() if e == element]
        return sum(count)

    def subscript(self, index: int) -> Optional[int]:
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
        if self.subscripts == []:
            return self.comp_str
        latex_str = ''
        for i, c in enumerate(self.comp_str):
            if c.isdigit() or c == ' ':
                continue
            subs = self.subscript(i)
            prefix = '\\text{'
            if subs is None: # check delims to eliminate \text redundancy.
                latex_str += prefix + c + '}'
            else:
                latex_str += prefix + c + '}' + '_{' + str(subs) + '}'
        return latex_str


    def __eq__(self, other: Self) -> bool:
        return sorted(self.tokens, key=str) == sorted(other.tokens, key=str)
