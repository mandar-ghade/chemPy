from .compound import Compound
from .equation import Equation
import numpy as np
from typing import Optional


# eq = Equation.parse_from_string('NO + O2 -> NO2')
# eq.h_rxn = -116
# eqns = []
# eqns.append((Equation.parse_from_string('N2 + O2 + H2O -> HNO3'), -256))
# eqns.append((Equation.parse_from_string('N2 + O2 -> 2NO'), 183))
# eq.extend(eqns, Equation.parse_from_string('NO2 + H2O -> HNO3 + NO'))

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


class Hess_Law:
    def __init__(self, initial_eq: tuple[Equation, Optional[float | int]],
                intermediate_equations: list[tuple[Equation, Optional[float | int]]],
                desired_equation: Equation):
        self.initial_eq, self.initial_h_rxn = initial_eq
        self.initial_eq.h_rxn = self.initial_h_rxn
        self.desired_equation = desired_equation
        self.desired_equation.balance()
        self.initial_eq.balance()
        for intermediate_eq, h_rxn in intermediate_equations:
            # make capture coefficient method
            intermediate_eq.h_rxn = h_rxn
            # if any(comp.coefficient != 1 for comp in intermediate_eq.compounds):
            #     continue
            # subs_specified = False
            # for i, comp in enumerate(intermediate_eq.compounds):
            #     size, no_subs = extract_coefficient(comp.original_comp_str)
            #     if no_subs:
            #         continue
            #     print(size, comp.original_comp_str)
            #     subs_specified = True
            #     if size.is_integer():
            #         comp.coefficient = int(size)
            #     else:
            #         comp.coefficient = size
            #     intermediate_eq.coefficients[i] = comp.coefficient
            # if not subs_specified:
            intermediate_eq.balance()
        self.intermediate_equations = [eq for eq, _ in intermediate_equations]
        self.unique_compounds = set(compound 
                                    for eq in self.intermediate_equations + [self.initial_eq] + [self.desired_equation]
                                    for compound in eq.compounds)
        self.new_h_rxn = self._get_new_enthalpy()


    def _get_new_enthalpy(self) -> Optional[int | float]:
        self.A = []
        self.B = []
        self.h_rxns = np.array([equation.h_rxn for equation in [self.initial_eq] + self.intermediate_equations])
        for compound in self.unique_compounds:
            row = []
            for equation in [self.initial_eq] + self.intermediate_equations:
                row.append(equation._matching_coefficient(compound))
            self.A.append(row)
            self.B.append([self.desired_equation._matching_coefficient(compound)])
        self.solution = np.linalg.lstsq(self.A, self.B, rcond=None)[0].T[0]
        return np.sum(self.solution * self.h_rxns)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(({self.initial_eq}, {self.initial_h_rxn}), {",".join((f"[{eq.__repr__()}, {eq.h_rxn}]" for eq in self.intermediate_equations))}, {self.desired_equation.__repr__()})'

