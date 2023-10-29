import os
import json


ELEMENTS_PATH = os.path.dirname(os.path.realpath(__file__)) + '\\elements.json'
with open(ELEMENTS_PATH, 'r') as fp:
    ELEMENT_ITEMS: json = json.load(fp)
ELEMENT_PROTON_DATA = {element: i+1 for i, [element, _] in enumerate(ELEMENT_ITEMS.items())}
ELEMENTS = sorted(ELEMENT_ITEMS.keys(), key=len, reverse=True)