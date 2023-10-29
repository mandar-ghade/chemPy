from utils import tokenize
from .element import Element
from Collections import Counter
from .subscript import Subscript
from typing import Optional, Self

class Compound:
    def __init__(self, comp_str: str, tokens: list[Element] = None, subscripts: list[Subscript] = None):
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
        return f"{self.__class__.__name__}({self.tokens})"
    
    def _get_total_electrons(self) -> int:
        return sum(token.electrons * count 
                   for token, count in self.elements.items())
    
    def _get_molar_mass(self) -> int:
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
        return self.tokens == other.tokens