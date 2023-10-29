from typing import Optional


def split_str(line: str) -> Optional[str]:
    split_chars = ['=', '->', '→', '=>']
    for split_char in split_chars:
        if split_char not in line:
            continue
        return split_char
    return None