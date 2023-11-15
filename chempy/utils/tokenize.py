from ..data import *
from ..element import Element
from ..subscript import Subscript
from typing import Optional, Self
from collections import Counter


def nests_delims(comp_str: str, start: int, end: int) -> bool:
    return any(char in ALL_DELIMS for char in comp_str[start+1:end-1])


def get_element_count(comp_str: str, start: int) -> tuple[int, bool]:
    subs = ''
    no_subscript = False
    for n in range(start, len(comp_str)+1):
        digit = comp_str[n:n+1]
        if not digit.isdigit():
            break
        subs += digit
    if subs == '':
        subs = 1
        no_subscript = True
    return int(subs), no_subscript


def matching_delim_index(delim_list: list, delim: str) -> Optional[int]:
    for i, ldm in enumerate(delim_list):
        if ldm != delim:
            continue
        return i
    return None


def retrieve_delims(
    comp_str: str, 
    ldi: Optional[int] = None,
    rdi: Optional[int] = None,
) -> tuple[Optional[int], Optional[int]]:
    right_count = 0
    right_delim: Optional[str] = None
    if ldi is None and rdi is not None:
        return None, None
    assert (isinstance(ldi, int) and isinstance(rdi, int))
    ld, rd = None, None
    for i in range(ldi, rdi):
        c = comp_str[i]
        if c in LEFT_DELIMS and ld is None and rd is None:
            ld = i
            matching_rdi = matching_delim_index(LEFT_DELIMS, c)
            if matching_rdi is None:
                raise TypeError(f'Right delim indices for {c} in {comp_str} could not be found.')
            right_delim = RIGHT_DELIMS[matching_rdi]
            continue
        if c in LEFT_DELIMS:
            right_count -= 1
            continue
        if c in RIGHT_DELIMS:
            right_count += 1
        if right_count > 0 and c == right_delim:
            rd = i
            break
    return (ld, rd)


def fetch_multiplier_list(
        comp_str: str,
        ldi: Optional[int],
        rdi: Optional[int], 
        subs = None,
        p_list: Optional[list] = None,
) -> list[tuple[int, int, int]]:
    if subs is None:
        subs = 1
    if p_list is None:
        p_list = []
    ldi, rdi = retrieve_delims(comp_str, ldi, rdi)
    if ldi is not None and rdi is None:
        raise Exception('Mismatched delims error. Compound={}'.format(comp_str))
    if ldi is None or rdi is None:
        return p_list
    multiplier, no_subs = get_element_count(comp_str, rdi+1)
    nested = nests_delims(comp_str, ldi+1, rdi-1)
    original_subs = subs
    subs *= multiplier
    p_list.append((ldi, rdi, subs))
    if nested:
        ldi += 1
        rdi -= 1
    else:
        ldi = rdi + 1 + (len(str(multiplier)) if no_subs is False else 0)
        rdi = len(comp_str)
        subs = original_subs
    return fetch_multiplier_list(comp_str, ldi, rdi, subs, p_list)


def get_multiplier(multiplier_list: list[tuple[int, int, int]], n: int) -> int:
    mp_list: list[int] = [multiplier
                       for ldi, rdi, multiplier in multiplier_list
                       if ldi <= n <= rdi]
    return mp_list[-1] if len(mp_list) > 0 else 1


def tokenize(comp_str: str) -> tuple[list[Element], list[Subscript]]:
    elements: list[Element] = []
    subs_list: list[Subscript] = []
    multiplier_list = fetch_multiplier_list(comp_str, 0, len(comp_str))
    multipliers = [get_multiplier(multiplier_list, n) for n in range(len(comp_str))]
    for i, char in enumerate(comp_str):
        if char in ALL_DELIMS:
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
        for _ in range(count * multipliers[i]):
            elements.append(Element(element))
    return (elements, subs_list)