import json
import os
import pandas as pd
from typing import Optional, Union


ELEMENTS_PATH = os.path.dirname(os.path.realpath(__file__)) + '/elements.json'
ELEMENT_INFO_PATH = os.path.dirname(os.path.realpath(__file__)) + '/element_info.csv'
ELEMENT_CSV_DATA = pd.read_csv(ELEMENT_INFO_PATH)
df: pd.DataFrame = pd.DataFrame(ELEMENT_CSV_DATA)
electronegativity_data: Optional[float] = lambda protons : (
    dt.item() 
    if not (
        dt := (df.loc[df['AtomicNumber'] == protons, 'Electronegativity'])
    ).isnull().item()
    else None
)
with open(ELEMENTS_PATH, 'r') as fp:
    ELEMENT_ITEMS: dict[str, str] = dict(json.load(fp))
ELEMENT_PROTON_DATA = {element: i+1 for i, [element, _] in enumerate(ELEMENT_ITEMS.items())}
ELEMENTS = sorted(ELEMENT_ITEMS.keys(), key=len, reverse=True)