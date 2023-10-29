from typing import Optional, Literal

type Shape = Literal['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o']
MAX_SUBSHELL = {
    's': 2, 
    'p': 6, 
    'd': 10, 
    'f': 14
}
SUBSHELL_MAP: list[Shape] = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o']
SPECIAL_SUBSHELLS: list[Shape] = ['d', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o']


class Orbital:
    def __init__(self, n: int, l: int, electrons: int) -> None:
        self.n = n
        self.l = l
        self.electrons = electrons
        self.shape = SUBSHELL_MAP[l]

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.n}, {self.l}, {self.electrons})'
    
    def __str__(self) -> str:
        return f'{self.n}{self.shape}^{self.electrons}'
    

def total_electrons(l: int, all_electrons: int, proton_count: int) -> int:
    electrons = 4*abs(l) + 2
    excess = all_electrons + electrons - proton_count
    if excess > 0:
        electrons -= excess
    if electrons <= 0:
        return 0
    return electrons


def get_special_orbitals(
        n: int, l: int,
        left_over: list[tuple[int, str]],
    ) -> tuple[bool, list[tuple[int, str]]]:
    matches = [
        (n_, s_) for n_, s_ in reversed(left_over)
        if (n_ == n - 1 and s_ == SUBSHELL_MAP[l + 1]) or
        (n_ == n - 2 and s_ == SUBSHELL_MAP[l + 2]) or
        (n_ == n - 3 and s_ == SUBSHELL_MAP[l + 3])
    ]
    return len(matches) > 0, matches


def get_electron_config(
        protons: int,
        pqn: int = 1,
        left_over: Optional[list[tuple[int, str]]] = None,  # tuple[n, s]
        e_config: Optional[list[str]] = None,
        count: int = 0,
    ) -> list[str]:
    e_config = e_config or []

    if protons == count:
        return e_config

    original_left_over = left_over.copy() if left_over is not None else None
    shapes = [SUBSHELL_MAP[l] for l in range(pqn)]
    left_over = [s for s in shapes if s in SPECIAL_SUBSHELLS]
    shapes = sorted(set(shapes) - set(left_over), reverse=True)
    left_over = [(pqn, l_) for l_ in left_over]
    l_list = sorted(set(range(pqn)) - {2, 3, 4})

    lanthanide_exceptions: list[tuple[int, str, int]] = []
    five_f_exception = None
    th_exception = None
    lr_exception = None
    for shape, aqn in zip(shapes, l_list):
        if left_over != original_left_over and original_left_over is not None:
            matches_exist, matches = get_special_orbitals(pqn, aqn, original_left_over)
            if matches_exist:
                for n, s in matches:
                    l = [i for i, s_val in enumerate(SUBSHELL_MAP) if s_val == s][0]
                    e = total_electrons(l, count, protons)
                    count += e
                    if n in (3, 4) and s in ('d', 'f') and e in (1, 2, 4, 8, 9):
                        if s == 'd' and e in (4, 9):
                            e_config[-1] = f'{pqn}s^1'
                            e += 1
                        elif n == 4 and s == 'f' and e in (1, 2, 8):
                            lanthanide_exceptions.append((n, s, e))
                            e -= 1
                    if n == 5 and s in ('d', 'f') and e <= 9:
                        if s == 'd' and e in (8, 9):
                            e_config[-2] = '6s^1'
                            e += 1
                        elif n == 5 and s == 'f' and e in (1, 2, 3, 4, 5, 8):
                            if e == 2:
                                th_exception = True
                                e -= 2
                            else:
                                five_f_exception = True
                                e -= 1
                    if n == 6 and s == 'd' and e == 1:
                        lr_exception = True
                        e = 0
                    original_left_over.remove((n, s))
                    if e:
                        e_config.append(f'{n}{s}^{e}')
        electrons = total_electrons(aqn, count, protons)
        count += electrons
        if electrons:
            e_config.append(f'{pqn}{shape}^{electrons}')

    if five_f_exception is not None:
        e_config.append('6d^1')

    if th_exception is not None:
        e_config.append('6d^2')

    if lr_exception is not None:
        e_config.append('7p^1')

    e_config.extend([f'{n_+1}d^1' for n_, _, _ in lanthanide_exceptions])

    left_over += original_left_over or []

    return get_electron_config(protons, pqn+1, left_over, e_config, count)


def calculate_valence(segment: list[str]) -> int:
    count = 0
    for orbital in segment:
        shape = orbital[1]
        electrons = int(orbital[3:])
        if (shape in SPECIAL_SUBSHELLS 
            and electrons == MAX_SUBSHELL[shape]):
            continue
        count += electrons
    return count


def get_valence_electrons(e_cfg: list[str]) -> int:
    pqns = [o for o in e_cfg if o[1] == 's']
    max_pqn_s_index = e_cfg.index(pqns[-1])
    segment = e_cfg[max_pqn_s_index:]
    valence_es = calculate_valence(segment)
    return valence_es


def main():
    pass


if __name__ == '__main__':
    main()