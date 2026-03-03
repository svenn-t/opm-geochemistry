#%%
import sys
import unittest

from pyGeoChem.io_cpp_solver import Species

class SolutionTest(unittest.TestCase):
    '''
    Unit tests to test solution class
    '''

    def test_basis(self):
        specie = Species()
        i1 = specie.get_specie('Ca')
        i2 = specie.get_specie('Ca+2')
        self.assertEqual(i1,i2)

        i1 = specie.get_specie('Mg')
        i2 = specie.get_specie('Mg+2')
        self.assertEqual(i1,i2)

        i1 = specie.get_specie('Cl')
        i2 = specie.get_specie('Cl-')
        self.assertEqual(i1,i2)
    
    def test_parse_reaction(self):
        specie = Species()
        names, stoich = specie.parse_reaction('Mg+2 + 2Cl-')
        self.assertEqual(names, ['Mg+2', 'Cl-'])
        self.assertEqual(stoich, [1, 2])

        names, stoich=specie.parse_reaction(' Mg+2  +   2Cl-')
        self.assertEqual(names, ['Mg+2', 'Cl-'])
        self.assertEqual(stoich, [1, 2])

        names, stoich=specie.parse_reaction(' Ca  - H + HCO3')
        self.assertEqual(names,['Ca','H','HCO3'])
        self.assertEqual(stoich,[1,-1,1])

    def test_solution_composition(self):
        seawater_comp={'Na+':0.4500870487,'Cl-':0.52512319,
                       'SO4-2':0.02400696, 'HCO3-':0.0020236,
                       'K+':0.010060186, 'Mg+2':0.0445150,'Ca+2':0.0129917}
                       
        specie = Species()
        specie.add_salt('NaCl', reaction='Na + Cl', mass=23.38)
        specie.add_salt('Na2SO4',reaction='2Na+ + SO4-2',mass=3.41)
        specie.add_salt('NaHCO3',reaction='Na+ + HCO3-',mass=0.17)
        specie.add_salt('KCl', reaction='K+ + Cl-',mass=0.75)
        specie.add_salt('MgCl2x6H2O',reaction='Mg+2 + 2Cl- + 6H2O',mass=9.05)
        specie.add_salt('CaCl2x2H2O',reaction='Ca+2 + 2Cl- + 2H2O',mass=1.91)
        specie.calculate_solution_composition()
        for spec,conc in seawater_comp.items():
            self.assertAlmostEqual(specie.solution_species[spec][specie.conc_mol_str],conc)



if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)    
# %%
