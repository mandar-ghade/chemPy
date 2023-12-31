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
    
    def exception(self, electrons: int) -> Optional[str]:
        if self.is_lr_exception(electrons):
            return 'lr'
        elif self.is_th_exception(electrons):
            return 'th'
        elif self.is_five_f_exception(electrons):
            return '5f'
        elif self.is_five_d_exception(electrons):
            return '5d'
        elif self.is_unstable_transition_metal(electrons):
            return 'utm'
        elif self.is_lanthanide_exception(electrons):
            return 'lanthanide'

    def is_lr_exception(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` represents `Lawrencium` in its unstable electron configuration."""
        return self.n == 6 and self.l == 2 and electrons == 1

    def is_th_exception(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` represents `Thorium` in its unstable electron configuration."""
        return self.n == 5 and self.l == 3 and electrons == 2

    def is_five_f_exception(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` has an unstable electron configuration in `5f`."""
        return self.n == 5 and self.l == 3 and electrons in (1, 3, 4, 5, 8)

    def is_five_d_exception(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` has an unstable electron configuration in `5d` (`n`=`5`, `l`=`2`)."""
        return self.n == 5 and self.l == 2 and electrons in (8, 9)

    def is_unstable_transition_metal(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` has an unstable electron configuration in its `d` block."""
        return self.n == 3 and self.l == 2 and electrons in (4, 9)

    def is_lanthanide_exception(self, electrons: int) -> bool:
        """Returns `True` if the `Orbital` has an unstable lanthanide electron configuration."""
        return self.n == 4 and self.l == 3 and electrons in (1, 2, 8)