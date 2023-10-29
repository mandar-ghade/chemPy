from typing import Optional, Self
class Subscript:
    def __init__(self, comp_str: str, 
                 element_index: Optional[int] = None, 
                 start_index: Optional[int] = None, 
                 size: Optional[int] = None) -> None:
        assert isinstance(comp_str, str)
        self.comp_str = comp_str
        self.element_index = element_index
        if start_index is None and size is None:
            raise Exception(f'Cannot derive subscript size from {comp_str}')
        self.start_index = start_index
        self.element_str = self.comp_str[element_index:start_index]
        self.size, self.subs = self._get_sizes()

    def _get_sizes(self) -> list[int, Optional[int]]:
        size = ''
        for n in range(self.start_index, len(self.comp_str)+1):
            digit = self.comp_str[n:n+1]
            if not digit.isdigit():
                break
            size += digit
        subs = int(size) if size != '' else None
        if size == '':
            size = 1
        size = int(size)
        return [size, subs]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(comp_str='{self.comp_str}', element_index={self.element_index}, start_index={self.start_index}, size={self.size})"

    def __hash__(self) -> int:
        return hash((self.element_str, self.subs, self.element_index))
    
    def __eq__(self, other: Self) -> bool:
        return (self.element_str, self.subs, self.element_index) == (other.element_str, other.subs, other.element_index)