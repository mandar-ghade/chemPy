from collections import Counter
import fractions
import json
from IPython.display import display, Markdown
import os
import pandas as pd
from string import ascii_lowercase
import sympy as smp
from typing import List, Self, Optional

ELEMENTS_PATH = os.path.dirname(os.path.realpath(__file__)) + '\\elements.json'
LEFT_DELIMS = ['[', '(', '{']
RIGHT_DELIMS = [']', ')', '}']

with open(ELEMENTS_PATH, 'r') as fp:
    ELEMENT_ITEMS: json = json.load(fp)

ELEMENTS = sorted(ELEMENT_ITEMS.keys(), key=len, reverse=True)


class Subscript:
    def __init__(self, comp_str: str, element_index = None, start_index = None, size = None):
        assert isinstance(comp_str, str)
        self.comp_str = comp_str
        self.element_index = element_index
        if start_index is None and size is None:
            raise Exception(f'Cannot derive subscript size from {comp_str}')
        self.start_index = start_index
        self.element_str = self.comp_str[element_index:start_index]
        self.size, self.subs = self._get_sizes()

    def _get_sizes(self) -> list[int, Optional[int]]:
        size = ''
        for n in range(self.start_index, len(self.comp_str)+1):
            digit = self.comp_str[n:n+1]
            if not digit.isdigit():
                break
            size += digit
        subs = int(size) if size != '' else None
        if size == '':
            size = 1
        size = int(size)
        return [size, subs]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}([element='{self.element_str}', subs={self.subs}, index={self.element_index}])"

    def __hash__(self) -> int:
        return hash((self.element_str, self.subs, self.element_index))
    
    def __eq__(self, other: Self) -> bool:
        return (self.element_str, self.subs, self.element_index) == (other.element_str, other.subs, other.element_index)


class Element:
    def __init__(self, symbol: str, amu = None):
        assert isinstance(symbol, str)
        self.symbol = symbol
        self.molar_mass = float(ELEMENT_ITEMS[self.symbol])
        
    def __str__(self) -> str:
        return self.symbol
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.symbol}')"
    
    def __hash__(self) -> int:
        return hash(self.symbol)
    
    def __eq__(self, other: Self) -> bool:
        return self.symbol == other.symbol


class Token(Element):
    def __init__(self, element_str: str):
        super().__init__(element_str)


#work on __eq__
class Tokenize:
    def __init__(self, comp_str: str):
        assert isinstance(comp_str, str)
        self.comp_str = comp_str
        self.index = 0
        self.multiplier_list = self._multipliers(0, len(self.comp_str))
        self.multipliers = [self._get_multiplier(n) for n in range(len(self.comp_str))]
        self.elements: Optional[Counter[Token]] = None
        self.compound: tuple[list[Token], list[Subscript]] = self._parse_elements()
        self.molar_mass = self._get_molar_mass()

    def _nests_parenthesis(self, s, e) -> bool:
        return any(c in LEFT_DELIMS + RIGHT_DELIMS 
                    for c in self.comp_str[s+1:e-1])

    def _get_count(self, r: int) -> list[int, bool]:
        subs = ''
        no_subs = False
        for n in range(r, len(self.comp_str)+1):
            digit = self.comp_str[n:n+1]
            if not digit.isdigit():
                break
            subs += digit
        if subs == '':
            subs = 1
            no_subs = True
        return [int(subs), no_subs]

    def _matching_delim(self, delim_list: list, delim: str) -> Optional[str]:
        for i, ldm in enumerate(delim_list):
            if ldm != delim:
                continue
            return i

    def _retrieve_delims(self, ldi = None, rdi = None) -> None:
        right_count = 0
        right_delim = None
        if ldi is None and rdi is None:
            return [None, None]
        ld, rd = None, None
        for i in range(ldi, rdi):
            c = self.comp_str[i]
            if c in LEFT_DELIMS and ld is None and rd is None:
                ld = i
                right_delim = RIGHT_DELIMS[self._matching_delim(LEFT_DELIMS, c)]
                continue
            if c in LEFT_DELIMS:
                right_count -= 1
                continue
            if c in RIGHT_DELIMS:
                right_count += 1
            if right_count > 0 and c == right_delim:
                rd = i
                break
        return [ld, rd]

    def _multipliers(self, ldi: int, rdi: int, subs = None, p_list: list = None) -> list[tuple[int, int, int]]:
        if subs is None:
            subs = 1
        if p_list is None:
            p_list = []
        ldi, rdi = self._retrieve_delims(ldi, rdi)
        if ldi is not None and rdi is None:
            raise Exception('Mismatched delims error. Compound={}'.format(self.comp_str))
        if ldi is None or rdi is None:
            return p_list
        multiplier, no_subs = self._get_count(rdi+1)
        nested = self._nests_parenthesis(ldi+1, rdi-1)
        original_subs = subs
        subs *= multiplier
        p_list.append([ldi, rdi, subs])
        if nested:
            ldi += 1
            rdi -= 1
        else:
            ldi = rdi + 1 + (len(str(multiplier)) if no_subs is False else 0)
            rdi = len(self.comp_str)
            subs = original_subs
        return self._multipliers(ldi, rdi, subs, p_list)
        
    def _get_multiplier(self, n: int) -> int:
        multiplier_list = [multiplier for ldi, rdi, multiplier in self.multiplier_list if ldi <= n <= rdi]
        return multiplier_list[-1] if len(multiplier_list) > 0 else 1
    
    def _parse_elements(self) -> tuple[list[Token], list[Subscript]]:
        elements = []
        subs_list: list[Subscript] = []
        comp_str = self.comp_str
        for i, char in enumerate(comp_str):
            multiplier = self.multipliers[i]
            if char in LEFT_DELIMS + RIGHT_DELIMS:
                subscript = Subscript(comp_str, i, i + 1)
                subs_list.append(subscript)
            if char not in ELEMENTS and comp_str[i:i+2] not in ELEMENTS:
                continue
            two_letter = comp_str[i:i+2] in ELEMENTS
            start = i + 2 if two_letter else i + 1
            subscript = Subscript(comp_str, i, start)
            count = subscript.size
            element = comp_str[i:i+2] if two_letter else char
            subs_list.append(subscript)
            for _ in range(count * multiplier):
                elements.append(Token(element))
        self.elements = Counter(elements)
        self.subs_list = subs_list
        return (elements, subs_list)
    
    def _get_molar_mass(self) -> float:
        return sum(token.molar_mass * count
                   for token, count in self.elements.items())

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.compound})'
    
    def __iter__(self) -> Self:
        self.index = 0
        return self
    
    def __next__(self) -> int:
        tk = [token for i, token in enumerate(self.compound.tokens) if i == self.index]
        if tk == []:
            raise StopIteration
        tk = tk[0]
        self.index += 1
        return tk


class Compound:
    def __init__(self, comp_str: str, tokens: list[Token] = None, subscripts: list[Subscript] = None):
        self.tokenized_comp = Tokenize(comp_str)
        if tokens is None:
            tokens, subscripts = self.tokenized_comp.compound
        self.tokens = tokens
        self.elements = Counter(tokens)
        self.comp_str = comp_str
        self.subscripts = subscripts
        self.molar_mass = self.tokenized_comp.molar_mass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.tokens})"
    
    def count(self, element) -> int:
        count = [count for e, count in Counter(self.tokens).items() if e == element]
        if len(count) == 0:
            return 0
        return sum(count)

    def subscript(self, index: int) -> Optional[int]:
        if self.subscripts is None:
            return None
        subs = [sub.subs for sub in self.subscripts if sub.element_index == index and sub.subs is not None]
        if len(subs) == 0:
            return None
        elif len(subs) > 1:
            raise Exception(f'Multiple subscripts show under comp_str={self.comp_str}, index={index}')
        return int(subs[0])

    def latexify(self) -> str:
        if self.subscripts == []:
            raise Exception('Subscripts not inputted correctly.')
        latex_str = ''
        for i, c in enumerate(self.comp_str):
            if c.isdigit():
                continue
            subs = self.subscript(i)
            if subs is None:
                latex_str += c
            else:
                latex_str += c + '_{' + str(subs) + '}'
        return latex_str


    def __eq__(self, other: Self) -> bool:
        return self.tokens == other.tokens
    



class Equation:
    def __init__(self, reactants: list[Compound], products: list[Compound]):
        assert isinstance(reactants, list)
        assert isinstance(products, list)
        self.reactants = reactants
        self.products = products
        self.coefficients = [1 for _ in range(len(self.reactants) + len(self.products))]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}([{self.reactants}], [{self.products}])"

    def _equation(self, latex = None) -> str:
        return ' + '.join(
            [
                coef + f'({reactant.latexify() if latex not in [None, False] else reactant.comp_str.strip()})' 
                if int(str(coef).replace(',','')) > 1 
                    else (
                        reactant.latexify() if latex not in [None, False]
                        else reactant.comp_str.strip()
                    )
                for reactant, coef in 
                zip(self.reactants, self.coefficients[:len(self.reactants)])
            ]) \
        + ' → ' + \
            ' + '.join(
                [
                    coef + f'({product.latexify() if latex not in [None, False] else product.comp_str.strip()})' 
                    if int(str(coef).replace(',','')) > 1 
                        else (
                            product.latexify() if latex not in [None, False]
                            else product.comp_str.strip()
                        )
                    for product, coef in zip(self.products, self.coefficients[len(self.reactants):])
                ])


    def _repr_mimebundle_(self, **kwargs) -> dict:
        return {
            "text/plain": f'${self._equation()}$',
            "text/latex": f'${self._equation(latex=True)}$'
        }
    

    def __str__(self) -> str:
        return self._equation()
    
    def total_left(self) -> Counter[Token]:
        '''Returns total elements to the left of the equation (reactants).'''
        reactants = Counter()
        for reactant in self.reactants:
            compound = reactant.elements.most_common()
            for (element, count) in compound:
                reactants[element] += count 
        return reactants

    def total_right(self) -> Counter[Token]:
        '''Returns total elements to the right of the equation (products).'''
        products = Counter()
        for product in self.products:
            compound = product.elements.most_common()
            for (element, count) in compound:
                products[element] += count 
        return products

    def count_left(self, element: Token) -> list[int]:
        '''Counts occurances of an element out of all the reactants.'''
        return [compound.count(element) for compound in self.reactants]

    def count_right(self, element: Token) -> list[int]:
        '''Counts occurances of an element out of all the products.'''
        return [compound.count(element) for compound in self.products]
    
    def _get_coefficients(self, matrix: list) -> list:
        '''Returns a list of coefficients, whose indices correspond to the positions of compounds in the given equation.'''
        matrix = smp.Matrix(matrix)
        length = matrix.shape[1]
        solutions = smp.linsolve(matrix, [smp.Symbol(f'{ascii_lowercase[n]}') for n in range(length)])
        for n in range(length):
            solutions = solutions.subs(smp.Symbol(f'{ascii_lowercase[n]}'), 1)
        if solutions.args == ():
            print('Impossible equation')
            return self.coefficients
        args = solutions.args[0]
        args: list = [fractions.Fraction(abs(arg)).limit_denominator() for arg in args]
        denom_list = [frac.denominator for frac in args if frac.denominator != 1]
        if denom_list != []:
            while not all(arg.denominator == 1 for arg in args):
                denom_list = [frac.denominator for frac in args if frac.denominator != 1]
                args = [arg * max(denom_list) for arg in args]
        args = [f'{arg.numerator:,}' for arg in args]
        return args

    def balance(self) -> None:
        '''Balances the chemical equation.'''
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

def main():
    pass


if __name__ == '__main__':
    main()