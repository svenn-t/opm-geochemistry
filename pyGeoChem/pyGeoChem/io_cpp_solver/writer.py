# %%
import os
import pathlib as pt
import pandas as pd
import numpy as np


class Input:
    '''
    For writing input files to the geochemical solver
    '''

    def __init__(self, name, dir=''):
        '''
        name of file to run
        '''
        self.phase_dict = {}  # water properties
        self.aq_dict = {}  # ion composition, minerals etc
        self.reset_input()
        self.relative_dir = dir

        self.full_file_name = dir + str(name) + '.dat'

    def reset_input(self):
        '''
        initialize dictionaries
        '''
        self.phase_dict = {'GEOCHEM': {'Temp': 25, 'Pres': 1.01325e5, 'debug': 0}}  # water properties

        self.aq_dict = {'SOLUTION': {}, 'EQUILIBRIUM_PHASES': {}, 'GAS_PHASE': {}}  # ion composition  # minerals

    def set_temp(self, T):
        self.phase_dict['GEOCHEM']['Temp'] = T

    def set_pres(self, P):
        self.phase_dict['GEOCHEM']['Pres'] = P

    def set_specie_modification(self, full_path_to_specie_file=''):
        self.phase_dict['GEOCHEM']['add_species'] = '\"' + full_path_to_specie_file + '\"'

    def set_jacobian(self, numerical=1):
        self.phase_dict['GEOCHEM']['NumericalJacobian'] = numerical

    def set_keyword(self, key, val):
        '''
        key: SOLUTION, EQUILIBRIUM_PHASES
        val: Ion concentration, mineral wt% saturation index
        '''
        self.phase_dict[key] = val

    def add_aq(self, key, val):
        '''
        key: Ion name
        val: conc molar
        '''
        self.aq_dict['SOLUTION'][key] = val

    def remove_aq(self, key):
        '''
        key: Ion name
        removes regardless if specie is present or not
        '''
        self.aq_dict['SOLUTION'].pop(key, None)

    def add_mineral(self, key, wtp=1, si=0):
        '''
        key: mineral name
        wtp: wt% of mineral
        si: saturation index or partial pressure
        '''
        self.aq_dict['EQUILIBRIUM_PHASES'][key] = str(wtp) + ' ' + str(si)

    def add_gas_component(self, key, pressure, vol_frac=1):
        '''
        key: gas phase name
        pressure: partial pressure of gas in Pascal
        vol_frac: gas volume/water volume [dimensionless]
        '''
        self.aq_dict['GAS_PHASE'][key] = str(pressure)
        if vol_frac != 1:
            self.aq_dict['GAS_PHASE']['volume_frac'] = str(vol_frac)

    def write_block(self, f, out_dict, end='/end\n'):
        '''
        f: file object
        out_dict: a dictionary to write to file,
                  if dictionary contains another dictinary, function
                  is called recursively
        end =  token to indicate end of write block
        '''
        for key in out_dict:
            if type(out_dict[key]) == dict:
                f.write(key + '\n')
                self.write_block(f, out_dict[key], end=end)
            else:
                f.write(key + '\t' + str(out_dict[key]))
                f.write('\n')
        f.write(end)

    def write_run_file(self):
        end = '/end'
        with open(self.full_file_name, 'w') as f:
            self.write_block(f, self.phase_dict, end='')
            for key in self.aq_dict:
                if self.aq_dict[key]:  # Dont write if empty
                    f.write(key + '\n')
                    self.write_block(f, self.aq_dict[key])
            f.write(end)

    def init_seawater(self):
        self.add_aq('pH', '7 charge')
        self.add_aq('SO4', '0.024')
        self.add_aq('Na', '0.45')
        self.add_aq('Mg', 0.0445)
        self.add_aq('Cl', 0.525)
        self.add_aq('Ca', '0.013')
        self.add_aq('K', 0.01)
        self.add_aq('HCO3', '1e-8 as CO3')

    def init_cc_co2(self):
        self.add_mineral('calcite')
        self.add_mineral('CO2(G)', wtp=1, si=-3.5)


if __name__ == '__main__':
    print(os.getcwd())
    try:
        os.chdir('../../../debug/')
    except:
        pass
    input_file = 'ind_py'
    GC = Input(input_file)
    GC.init_seawater()
    GC.init_cc_co2()
    GC.write_run_file()

#    df_aq=pd.read_csv(input_file+'_aq.out',sep='\t')
#    df_buf=pd.read_csv(input_file+'_buffer.out',sep='\t')
#    print(df_aq)
#    print(df_buf)
# %%


# %%
