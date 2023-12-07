def strip_coefficients(comp_str: str) -> str:
    comp_str = comp_str.strip()
    for i, c in enumerate(comp_str):
        if c.isdigit():
            continue
        return comp_str[i:]