from typing import Optional, Self
from .data import SUBSHELL_MAP


class Orbital:
    def __init__(self, n: int, l: int, electrons: Optional[int]) -> None:
        self.n = n
        self.l = l
        self.electrons = electrons
        self.shape = SUBSHELL_MAP[l]
        self.is_special = self.l >= 2

    def __str__(self) -> str:
        return f'{self.n}{self.shape}^{self.electrons}'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.n}, {self.l}, {self.electrons})'

    def __hash__(self):
        return hash((self.n, self.l, self.electrons))

    def __eq__(self, other: Self) -> bool:
        assert isinstance(other, self.__class__)
        return hash(self) == hash(other)