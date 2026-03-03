#%%
import pandas as pd 
import pathlib as pt 
from pyGeoChem.io_cpp_solver.util import  path_exists
import re
import os

class Species:
    '''
    A class for calculating concentrations of basis species in a solution
    User needs to specify disassociation in terms of basis species given in
    ../data/basis_species.csv. Example of reaction are
    MgCl2x6H20 = Mg + 2Cl + 6H2O
    '''
    def __init__(self,volume=1,ppm=False):
        self.volume=volume # unit: Liter
        self.ppm = ppm # mass entered in ppm unit
        absolute_path = os.path.dirname(__file__)
        self.data_basis_spec=os.path.join(absolute_path,'../data/basis_species.csv')
        self.elements=os.path.join(absolute_path,'../data/elements_reduced.csv')
        self.basis_str='Basis' #col name
        self.nickn_str='NickName'
        self.molw_str='MolWeight'# col name
        self.reac_str='reaction'
        self.conc_str='concentration [g/l]'
        self.conc_mol_str= 'concentration [mol/l]'
        self.stoi_str='stoichiometric_coefficient'
        self.mass_str='mass [g]'
        self.tds=0. # total dissolved solid
        assert path_exists(self.data_basis_spec)
        assert path_exists(self.elements)
        self.df_basis=pd.read_csv(self.data_basis_spec,sep='\t')
        assert len(self.df_basis[self.basis_str].unique()) == len(self.df_basis)        
        self.df_basis[self.nickn_str]=self.add_nick_name()
        assert len(self.df_basis[self.nickn_str].unique()) == len(self.df_basis)     

        self.solution_composition_calculated=False

        self.salts={} 
        self.solution_species={}
    
    def __repr__(self):
        head=[self.conc_str,self.mass_str,self.molw_str,self.stoi_str,]
        line='salt'
        for h in head:
                line += '\t' + h
        line += '\n'
        for key, salt_di in self.salts.items():
            line += key
            for h in head:
                line += '\t'+ str(salt_di[h])
            line += '\n'
        
        return line
    
    def __add__(self,solution2):
        '''
        add two solutions, returns a new object that contains the combined solution
        i.e. c_1 V_1, and c_2 V_2 final: (c_1*V_1+c_2*V_2)/V_+V_2
        '''
        assert solution2.solution_composition_calculated
        assert self.solution_composition_calculated

        new_specie = Species(self.volume+solution2.volume)
        new_specie.salts = self.salts.copy()
#        new_specie.solution_species = self.solution_species.copy()
        
        for salt,salt_d in solution2.salts.items():
            if salt in new_specie.salts:
                new_specie.salts[salt][self.mass_str] += salt_d[self.mass_str]
            else:
                new_specie.salts[salt]=salt_d
            new_specie.salts[salt][self.conc_str]=new_specie.salts[salt][self.mass_str]/new_specie.volume
        new_specie.add_all_basis_from_reactions()
        new_specie.calculate_solution_composition()

        return new_specie
    
    def match_species(self,specie_name):
        '''
        Try to find closest match of species_name in predefined basis
        '''
        poss_spec=[]
        for specie in self.df_basis[self.basis_str]:
            if specie_name in specie:
                poss_spec.append(specie)
        return poss_spec

    
    def solution(self,only_mol=False):
        '''
        returns a data frame of the ion composition of the solution
        '''
        if not self.solution_composition_calculated:
            self.calculate_solution_composition()
        sol_dict={'Ion':[], self.conc_str:[], self.conc_mol_str:[]}
        if only_mol:
            sol_dict={'Ion':[],self.conc_mol_str:[]}
        for key,sol_d in self.solution_species.items():
            sol_dict['Ion'].append(key) 
            for h in sol_dict:
                if h == 'Ion':
                    pass
                else:
                    sol_dict[h].append(sol_d[h])
#        line += '\n Total Dissolved Solid = ' + str(1e3*self.tds/self.volume) + ' ppm\n'
        return pd.DataFrame(sol_dict)


    
    def get_unique_specie_name(self,name):
        '''
        Basis species are defined in variable self.data_basis_spec
        It is possible to omit charge in the definition 
        '''
        idx=self.get_specie(name)
        return str(self.df_basis[self.basis_str].iloc[idx])
    
    def get_specie(self,name : str):
        '''
        Check if specie is a valid basis specie
        returns index 
        '''
        val1=self.df_basis[self.basis_str]==name
        val2=self.df_basis[self.nickn_str]==name
        assert val1.sum()+val2.sum()>0, name + " not a valid basis species"
        if val1.sum() > val2.sum():
            return val1[val1].index[0]
        else:
            return val2[val2].index[0]

    def add_nick_name(self):
        '''
        remove charge from basis specie
        returns pandas.Series
        '''
        df2=self.df_basis[self.basis_str].str.split('+',n=1,expand=True)
        df3=df2[0].str.split('-',n=1,expand=True)
        return df3[0]

   
    def parse_reaction(self,reac: str):
        '''
        splits reaction into basis species and coefficients, e.g.
        Mg+2 + 2Cl- returns ['Mg+2','Cl-'], [1,2]
        '''
        vec=reac.split(' ')
        names=[] # strings
        coeff=[] # numeric
        only_number_extracted=False
        for vi in vec:
            if not only_number_extracted:
                d,n=self.extract_first_number(vi)
                coeff.append(d)
                if len(n) >0:
                    names.append(n)
                else:
                    only_number_extracted=True
            else:
                names.append(vi)
                only_number_extracted=False
        nn=[]
        cc=[]
        sign=1
        for c,n in zip(names,coeff):
            if c=='':
                pass
            elif c=='+':
                sign=1
            elif c == '-':
                sign=-1
            else:
                nn.append(sign*n)
                cc.append(c)
                sign=1    
        return cc,nn
    
    def extract_first_number(self,chars : str):
        '''
        Extracts first float/integer of a string, not only 1..9
        '''
        ret_str=''
        sign=1
        i=0
        for i,c in enumerate(chars):
            if c.isnumeric():
                ret_str += c
            else:
                if c=='+':
                    sign=1
                elif c == '-':
                    sign =-1
                elif c == '.':
                    if len(ret_str)==0:
                        ret_str += '1.'
                    else:
                        ret_str += c
                else:
                    break
        if len(ret_str) == 0:
            return 1.*sign, chars      
        return sign*float(ret_str),chars[i:]

    def add_basis_specie_from_reaction(self,reaction):
        '''
        checks if all basis species are valid, and adds to basis dic
        '''
        basis,_=self.parse_reaction(reaction)
        valid_reaction = True
        for b in basis:
            # gives error if one specie is not valid
            i=self.valid_specie(b)           
        for b in basis:
            i=self.valid_specie(b)
            bn=str(self.df_basis[self.basis_str].iloc[i])
            mw=float(self.df_basis[self.molw_str].iloc[i])
            self.solution_species[bn]={self.conc_str:0.,self.molw_str:mw}
        return valid_reaction
    
    
    def valid_specie(self,specie_name):
        ''' Checks is specie is in database '''
        try:
            i = self.get_specie(specie_name)
            return i
        except ValueError:
            print('Not a valid basis specie: ' + str(b)+'\n Check file : ' + self.data_basis_spec)

    def get_molweight_specie(self,specie):
        try:
            i = self.get_specie(specie)
        except ValueError:
            print('Not a valid basis specie: ', specie)
            print('Check file :', self.data_basis_spec)
            return -1
        mw=float(self.df_basis[self.molw_str].iloc[i])
        return mw
        

    
    def add_all_basis_from_reactions(self):
        '''
        add all basis specises consistent with the salt dict
        '''
        for _,salt_d in self.salts.items():
            self.add_basis_specie_from_reaction(salt_d[self.reac_str])

    def calc_molw(self,salt_name):
        '''
        Calculates the molweight of the salt, based on reactions
        defined in self.add_salt
        '''
        if salt_name in self.salts:
            reaction=self.salts[salt_name][self.reac_str]
            basis,coeff=self.parse_reaction(reaction)
            self.salts[salt_name][self.stoi_str]=coeff
            mw=0.
            for b,c in zip(basis,coeff):
                bn=self.get_unique_specie_name(b)
                mw += c*self.solution_species[bn][self.molw_str]
            self.salts[salt_name][self.molw_str]=mw
        else:
            print(salt_name + ' not added yet!')

    def add_salt(self,name,*,reaction, mass):
        '''
        Give the name of the salt and a speciation reaction, e.g.
        add_salt('MgCl2x6H2O',reaction='Mg+2 + 2Cl- + 6H2O',mass=9.05)
        If basis species not valid, salt will not be added
        '''
        if self.ppm:
            mass_corr=mass*1e-3*self.volume
        else:
            mass_corr=mass
        self.salts[name]={self.reac_str:reaction, self.conc_str:mass_corr/self.volume,self.mass_str:mass_corr}
        if self.add_basis_specie_from_reaction(reaction):
            print('Succesfully added ', name)
            self.calc_molw(name)
        else:
            print(name, ' not added')
            del self.salts[name]
    def add_ion(self,ion,unit='ppm'):
        pos_spec=self.valid_specie(ion)


    def calculate_solution_composition(self):
        '''
        Iterates over all salts added in add_salt and
        calculates concentration of ions
        '''
        for key,salt_d in self.salts.items():            
            reaction=salt_d[self.reac_str]
            mw_salt=salt_d[self.molw_str]
            conc=salt_d[self.conc_str]
            basis,coeff=self.parse_reaction(reaction)
            for b,c in zip(basis,coeff):
                bn=self.get_unique_specie_name(b)
                mw_b=self.solution_species[bn][self.molw_str]
                self.solution_species[bn][self.conc_str] += conc*c*mw_b/mw_salt
        for key,sol_d in self.solution_species.items():
            if sol_d[self.molw_str]>0:
                sol_d[self.conc_mol_str]=sol_d[self.conc_str]/sol_d[self.molw_str]
            else:
                sol_d[self.conc_mol_str]=0.
        self.tds=0.
        for basis in self.solution_species:
            self.tds += self.solution_species[basis][self.conc_str]
        self.tds *= 1e3/self.volume
        self.solution_composition_calculated=True

if __name__=='__main__':
    fact=.5
    factb=.5
    A=Species(fact)
    B=Species(factb)
    A.add_salt('NaCl', reaction='Na + Cl', mass=23.38*fact)
    A.add_salt('Na2SO4',reaction='2Na+ + SO4-2',mass=3.41*fact)
    A.add_salt('NaHCO3',reaction='Na+ + HCO3-',mass=0.17*fact)
    A.add_salt('KCl', reaction='K+ + Cl-',mass=0.75*fact)
    A.add_salt('MgCl2x6H2O',reaction='Mg+2 + 2Cl- + 6H2O',mass=9.05*fact)
    A.add_salt('CaCl2x2H2O',reaction='Ca+2 + 2Cl- + 2H2O',mass=1.91*fact)

    B.add_salt('NaCl', reaction='Na + Cl', mass=23.38*factb)
    B.add_salt('Na2SO4',reaction='2Na+ + SO4-2',mass=3.41*factb)
    B.add_salt('NaHCO3',reaction='Na+ + HCO3-',mass=0.17*factb)
    B.add_salt('KCl', reaction='K+ + Cl-',mass=0.75*factb)
    B.add_salt('MgCl2x6H2O',reaction='Mg+2 + 2Cl- + 6H2O',mass=9.05*factb)
    B.add_salt('CaCl2x2H2O',reaction='Ca+2 + 2Cl- + 2H2O',mass=1.91*factb)


    A.calculate_solution_composition()
    B.calculate_solution_composition()
    print(A)
    print(A.solution())
    C=A+B
    print('----------C-------------')
    print(C)
    print(C.solution())


    


    
        







# %%
