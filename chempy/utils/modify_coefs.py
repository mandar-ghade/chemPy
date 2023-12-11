def strip_coefficients(comp_str: str) -> str:
    comp_str = comp_str.strip()
    for i, c in enumerate(comp_str):
        if not c.isalpha():
            continue
        return comp_str[i:]


def extract_coefficient(comp_str: str) -> tuple[float, bool]:
    subs = ''
    no_subscript = False
    for n in range(len(comp_str)+1):
        digit = comp_str[n:n+1]
        if digit.isalpha():
            break
        subs += digit
    if subs == '':
        subs = '1'
        no_subscript = True
    return eval(subs), no_subscript
    