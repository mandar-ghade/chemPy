from typing import Optional
from .data import MAX_SUBSHELL, SUBSHELL_MAP, SPECIAL_SUBSHELLS
from .orbital import Orbital


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
        left_over: list[Orbital],
    ) -> list[Orbital]:
    matches = [
        o for o in reversed(left_over)
        if (o.n == n - 1 and o.shape == SUBSHELL_MAP[l + 1]) or
        (o.n == n - 2 and o.shape == SUBSHELL_MAP[l + 2]) or
        (o.n == n - 3 and o.shape == SUBSHELL_MAP[l + 3])
    ]
    return matches


def get_electron_config(
        protons: int,
        pqn: int = 1,
        left_over: Optional[list[Orbital]] = None,
        e_config: Optional[list[Orbital]] = None,
        count: int = 0,
    ) -> list[Orbital]:

    e_config = e_config or []

    if protons == count: 
        return e_config

    _left_over: Optional[list[Orbital]] = left_over.copy() \
        if left_over is not None else None

    orbitals: list[Orbital] = [Orbital(pqn, l, None) for l in range(pqn)]
    non_special_orbitals: list[Orbital] = [orbital for orbital in orbitals if not orbital.is_special]
    new_left_overs: list[Orbital] = list(set(orbitals).difference(non_special_orbitals))

    lanthanide_exceptions: list[Orbital] = []
    five_f_exception = None
    th_exception = None
    lr_exception = None

    for orbital in non_special_orbitals:
        if new_left_overs != _left_over and _left_over is not None:
            matches = get_special_orbitals(orbital.n, orbital.l, _left_over)
            for s_o in matches:
                e = total_electrons(s_o.l, count, protons)
                count += e
                match s_o.exception(e):
                    case 'utm':
                        e_config[-1].n = pqn
                        e_config[-1].electrons = 1
                        e += 1
                    case 'lanthanide':
                        lanthanide_exceptions.append(s_o)
                        e -= 1
                    case '5d':
                        e_config[-2].n = 6
                        e_config[-2].shape = 's'
                        e_config[-2].electrons = 1
                    case 'th':
                        th_exception = True
                        e -= 2
                    case '5f':
                        five_f_exception = True
                        e -= 1
                    case 'lr':
                        lr_exception = True
                        e = 0
                _left_over.remove(s_o)
                if e:
                    s_o.electrons = e
                    e_config.append(s_o)
        electrons = total_electrons(orbital.l, count, protons)
        count += electrons
        if electrons:
            orbital.electrons = electrons
            e_config.append(orbital)
    
    if five_f_exception is not None:
        e_config.append(Orbital(6, 2, 1))

    if th_exception is not None:
        e_config.append(Orbital(6, 2, 2))

    if lr_exception is not None:
        e_config.append(Orbital(7, 1, 1))

    e_config.extend(Orbital(le.n+1, 2, 1) for le in lanthanide_exceptions)

    new_left_overs += _left_over or []

    return get_electron_config(protons, pqn+1, new_left_overs, e_config, count)


def calculate_valence(segment: list[Orbital]) -> int:
    count = 0
    for orbital in segment:
        shape = orbital.shape
        electrons = orbital.electrons
        if electrons is None:
            continue
        if (shape in SPECIAL_SUBSHELLS 
            and electrons == MAX_SUBSHELL[shape]):
            continue
        count += electrons
    return count


def get_valence_electrons(e_cfg: list[Orbital]) -> int:
    pqns = [o for o in e_cfg if o.shape == 's']
    max_pqn_s_index = e_cfg.index(pqns[-1])
    segment = e_cfg[max_pqn_s_index:]
    valence_electrons = calculate_valence(segment)
    return valence_electrons
