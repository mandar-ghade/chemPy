from ..data import ALL_DELIMS

def strip_coefficients(comp_str: str) -> str:
    comp_str = comp_str.strip()
    start = 0
    end = len(comp_str)
    for i, c in enumerate(comp_str):
        if not c.isalpha() and c not in ALL_DELIMS:
            continue
        start = i
        break
    for n, c in enumerate(comp_str):
        if n <= start:
            continue
        if c == '(':
            if not (comp_str[n+1:n+2].isalnum() 
                and comp_str[n+1:n+2].islower() 
                and comp_str[n+2:n+3] == ')'):
                continue
            end = n
            break
    return comp_str[start:end]


def extract_coefficient(comp_str: str) -> tuple[float, bool]:
    subs = ''
    no_subscript = False
    for n in range(len(comp_str)+1):
        digit = comp_str[n:n+1]
        if digit.isalpha() or digit in ALL_DELIMS:
            break
        subs += digit
    if subs == '':
        subs = '1'
        no_subscript = True
    return eval(subs), no_subscript
    