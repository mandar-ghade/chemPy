from .data import *
from .qm_model import get_valence_electrons, get_electron_config
from typing import Optional, Self


class Element:
    def __init__(self, symbol: str, mass: Optional[float] = None):
        if not isinstance(symbol, str):
            raise TypeError('Expected `str` for `symbol` argument but received '
                            f'`{symbol.__class__.__name__}` instead.')
        if symbol not in ELEMENTS:
            raise ValueError(f'"{symbol}" is not a valid `Element`.')
        assert isinstance(symbol, str)
        self.symbol = symbol
        self.molar_mass = float(ELEMENT_ITEMS[self.symbol])
        self.protons = ELEMENT_PROTON_DATA[self.symbol]
        self.e_cfg = get_electron_config(self.protons)
        self.electron_configuration = ' '.join(str(o) for o in self.e_cfg)
        self.valence_electrons = get_valence_electrons(self.e_cfg)
        self.electronegativity: Optional[float] = electronegativity_data(self.protons)     # type: ignore        
        self.electrons = self.protons
        self.mass = mass
        self.moles = None
        if self.mass is not None:
            self.moles = self.mass / self.molar_mass
        
    def __str__(self) -> str:
        return self.symbol

    def __add__(self, other: Self):
        from .compound import Compound
        if not isinstance(other, self.__class__) and \
            not isinstance(other, Compound):
            raise TypeError(f'Expected `Compound` or `Element` for `other` argument but received {other.__class__.__name__}')
        if isinstance(other, self.__class__):
            return Compound(self.symbol + other.symbol)
        elif isinstance(other, Compound):
            return Compound(self.symbol + other.comp_str)

    def __mul__(self, other: int):
        assert isinstance(other, int)
        from .compound import Compound
        return Compound(f'{self.symbol}{other}')
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.symbol}')"
    
    def __hash__(self) -> int:
        return hash(self.symbol)
    
    def __eq__(self, other: Self) -> bool:
        return self.symbol == other.symbol