from .element import Element
from .compound import Compound
from typing import Self, TypeVar, Any
from _collections_abc import Mapping


def add_elements_laterally(tuple_object) -> dict[Element, float | int]:
    element_counter: dict[Element, float | int] = dict()
    for element, count in tuple_object:
        element_get = element_counter.get
        element_counter[element] = element_get(element, 0) + count
    return element_counter


def remove_duplicate_compounds(compounds: list[Compound]) -> list[Compound]:
    new_comp_list: list[Compound] = []
    for i, compound in enumerate(compounds):
        comp = compound
        if comp in new_comp_list:
            continue
        for comp2 in compounds[i+1:]:
            if comp2 != comp:
                continue
            comp: Compound = comp + comp2
        new_comp_list.append(comp)
    return new_comp_list


def check_is_compound(obj: Any) -> None:
    if not isinstance(obj, Compound):
        raise TypeError(f'Expected `Compound` in `iterable` argument but received 'f'`{obj.__class__.__name__}` instead.')


def _count_elements(mapping, iterable):
    if isinstance(iterable, list):
        iterable = remove_duplicate_compounds(iterable)
    for item in iterable:
        check_is_compound(item)
        item: Compound = item
        mapping[item] = item.coefficient


class CompoundCounter(dict):
    def __init__(self, iterable, **kwds):
        super().__init__()
        self.elements = dict()
        self.update(iterable)
        self.elements = self._retrieve_elements()
        self.coefficients = [val for val in self.values()]

    def update(self, iterable, **kwds):
        if iterable is not None:
            if isinstance(iterable, Mapping):
                if self:
                    for comp, count in iterable.items():
                        check_is_compound(comp)
                        comp: Compound = comp
                        comp.coefficient *= count
                    iterable = remove_duplicate_compounds(list(iterable.keys()) + list(self.keys()))
                    self.clear()
                    for item in iterable:
                        check_is_compound(item) 
                        item: Compound = item
                        self[item] = item.coefficient
                    self.elements.update(self._retrieve_elements())
                    self.coefficients = [val for val in self.values()]
                else:
                    for comp, count in iterable.items():
                        check_is_compound(comp)
                        comp: Compound = comp
                        comp.coefficient *= count
                    iterable = remove_duplicate_compounds(list(iterable.keys()))
                    self.clear()
                    for item in iterable:
                        check_is_compound(item) 
                        item: Compound = item
                        self[item] = item.coefficient
                    self.elements.update(self._retrieve_elements())
                    self.coefficients = [val for val in self.values()]
            else:
                _count_elements(self, iterable)
        if kwds:
            self.update(kwds)
    
    def _retrieve_elements(self) -> dict[Element, int]:
        elements_dict: dict[Element, int] = dict()
        for cmp in self.keys():
            cmp: Compound = cmp
            elements_get = elements_dict.get
            for element, count in cmp.elements.items():
                elements_dict[element] = elements_get(element, 0) + count 
        return elements_dict 
    
    def search_element_occurances(self, element: Element) -> dict[Compound, int]:
        return dict((cmp, cmp.count(element)) for cmp in self)
    
    def is_equal(self, other: Self) -> bool:
        if not isinstance(other, CompoundCounter):
            raise TypeError(f'Expected `CompoundCounter` for `other` argument but received '
                            f'`{other.__class__.__name__}` instead.')
        return add_elements_laterally(tuple((element, count * compound.coefficient)
                       for compound, _ in self.items()
                       for element, count in compound.elements.items())) == \
               add_elements_laterally(tuple((element, count * compound.coefficient)
                       for compound, _ in other.items()
                       for element, count in compound.elements.items()))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({dict(self)})'
    
    def flatten(self) -> list[Compound]:
        return list(self.keys())

    def __add__(self, other: Self) -> Self:
        if not isinstance(other, CompoundCounter):
            raise TypeError(f'Expected `CompoundCounter` for `other` argument but received '
                            f'`{other.__class__.__name__}` instead.')
        for comp, count in other.items():
            check_is_compound(comp)
            comp: Compound = comp
            comp.coefficient = count
        iterable = remove_duplicate_compounds(list(other.keys()) + list(self.keys()))
        self.clear()
        for item in iterable:
            check_is_compound(item) 
            item: Compound = item
            self[item] = item.coefficient
        self.elements.update(self._retrieve_elements())
        self.coefficients = [val for val in self.values()]
        return self
    
    def __eq__(self, other: Self) -> bool:
        if not isinstance(other, CompoundCounter):
            raise TypeError(f'Expected `CompoundCounter` for `other` argument but received '
                            f'`{other.__class__.__name__}` instead.')
        return self.elements == other.elements
        