#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 17:21:45 2018

@author: wakatobi
"""
### Template from Ipthon Notebook doc_18_07_05

import sys
sys.path.append('/afs/itp.tugraz.at/user/wakatobi/Modul_neo/f90nml-0.21/')
import f90nml
import os
import yaml
from IPython import display


if __name__ == '__main__':
    print('I am main')


class Neo2_common_objects():
    """ objects appearing in both subclasses

    Attributes
    ----------
    path2code: str
        path to the repository
    path2exe:
        path to the executable
    rundir: str
        Path of the working directory
    req_files_names : dict
        dict of required filenames to to run neo2
        Only basic Set, extension of each subclass will be made
    req_files_paths : dict
        dict of paths of the required files to run neo2
        Only basic Set, extension of each subclass will be made
    template_path: str
        Path for sample files

    Methods
    -------
    set_neo2file(path):
        set a neo2.in file as sample
    set_template_path():
        Use all required Files from this directory.
    run_local():
        Run all Jobs local
    run_condor():
        Submit all Jobs to Condor
    """

    def __init__(self,wdir):


        self._neo2path=os.environ.get('NEO2PATH')
        # If NEO2PATH is not available the Attribute is NONE

        self.wdir=wdir
        if os.path.isdir(self._neo2path):
            self._path2code=self._neo2path
        else:
            self._path2code=''

        self._path2exe=''

        self.req_files_names={
        'neo2in': 'neo2.in',
        'neoin': 'neo.in'
        }
        #self._Load_default()
        self.req_files_paths=dict()

    def compile(self,overwrite=False):
        curdir=os.getcwd()
        os.chdir(self._path2code)
        exe=self._path2code+'NEO-2-QL/Build_auto/'
        if os.path.exists(exe):
            if overwrite==False:
                raise FileExistsError('Please remove folder Build_auto')
            else:
                os.chdir(exe)
                print(os.popen("cmake -DCMAKE_BUILD_TYPE=RELEASE ..").read())
                print(os.popen("make clean").read())
                print(os.popen("make -j4").read())
        else:
            os.mkdir(exe)
            os.chdir(exe)
            print(os.popen("cmake -DCMAKE_BUILD_TYPE=RELEASE ..").read())
            print(os.popen("make -j4").read())
        os.chdir(curdir)
        self.path2exe=exe+'neo_2.x'

    def _read_neo2in(self):
        try:
            self.neo2nml=Neo2File(self.req_files_paths['neo2in'])
        except:
            print('Couldn\'t read neo2.in')
            return




    def _createfiles(self,overwrite=False,link=True):
        """Create and/or link required files and folder into destination"""


        self._checkreqfiles()
        if not os.path.exists(self.path2exe):
            raise IOError('Path to Executable is missing')

        if self.singlerunpath==None:
            raise IOError('singlerunpath is missing')
        os.makedirs(self.singlerunpath,exist_ok=True)

        for i,j in self.req_files_paths.items():
            print(j)
            destfile=os.path.join(self.singlerunpath,os.path.basename(j))

            if os.path.isfile(destfile):
                if overwrite:
                    print(destfile, ' is a file and is replaced')
                    os.remove(destfile)
                else:
                    print(destfile, ' is not updated') ### Maybe check if it is the same file!!
                    continue
            if i == 'neo2in':
                if self.neo2nml.ischanged:
                    self.neo2nml.write(destfile)
                    continue
            os.symlink(j,destfile)



        i2=self.path2exe
        print('path2exe = ', i2)
        destfile=os.path.join(self.singlerunpath,os.path.basename(i2))
        print('destfile = ', destfile)
        if os.path.isfile(destfile):
            if overwrite:
                print(destfile, ' is a file and is replaced')
                os.remove(destfile)
            else:
                print(destfile, ' is not updated') ### Maybe check if it is the same file!!
                return

        os.symlink(i2,destfile)
        self._Runiscreated=True





    def _checkreqfiles(self):


        if set(self.req_files_names) != set(self.req_files_paths):
            raise ValueError('Not all paths of required Files defined')

        for i,j in self.req_files_paths.items():
            if not os.path.exists(j):
                raise RuntimeError(j,' is not a File')


        else:
            print('All required Files are existing')
            return





class Neo2Scan(Neo2_common_objects):
    """ Class for scans over one parameter


    Attributes
    ----------
    neo2in: Neo2File???
        sample inputfile for the neo2.in file
    neoin: str
        path to neo.in file. In future maybe own instance of a new class.
    unpert: str
        path to unpertubation file in Boozer coordinates. In future maybe own
        instance of a new class.
    pert: str
        path to pertubation file in Boozer cooridnates. In future maybe own
        instance of a new class.
    listofsingleruns: list
        A list of the directories of each Singlerun
    req_files_names : dict
        dict of required filenames to to run neo2
    req_files_paths : dict
        dict of paths of the required files to run neo2
    rundir: str
        Path of the working directory
    scan_param: str
        name of the parameter which will be scanned
    scan_values: list
        list of values which will be iterated
    path2code: str
        path to the repository
    path2exe:
        path to the executable
    structure: str
        string representation of folder structure


    Methods
    -------
    generate_SingleRuns()
        generate list of SingleRunInstances
    compile_code()
        compile code an link the executable
    generate_structure():
        generate structure for the singleruns
    set_neo2file(path):
        set a neo2.in file as sample
    run_local():
        Run all Jobs local
    run_condor():
        Submit all Jobs to Condor
    set_template_path():
        Use all required Files from this directory.

    """
###### Methods to be implementend: #########
#    def __init__(self,runpath):
#
#        self._neo2path=os.environ.get('NEO2PATH')
#        self.runpath=runpath
#        self.structure=self.runpath+'/lag/'
#        self._path2code=''
#        self._path2exe=''
#        self.scan_para=''
#        self.scan_values=list()
#        self.listofsingleruns=list()
#
#        self.req_files_names={
#        'neo2in': 'neo2.in',
#        'neoin': 'neo.in'
#        }
#        #self._Load_default()
#        self.req_files_paths=dict()


#    def _read_neo2in(self):
#        try:
#            self.neo2nml=Neo2File(self.req_files_paths['neo2in'])
#        except:
#            print('Couldn\'t read neo2.in')
#            return

#    def _Load_default(self):
#        pass
#
#    def set_pert_file(self):
#        """Method to set Pertubation File"""
#        pass
#
#    def set_axi_sym_file(self):
#        """Method to set axisymetric File"""
#        pass
#
#
#    def set_multispecies_file(self):
#        """Method to set Multispecies File"""
#        pass

    def __init__(self,wdir):
        super().__init__(wdir)
        self.structure=''
        self.scanparameter=''
        self.float_format='1.2e'



    def run_local(self):
        """Start Job on local machine"""
        pass

    def run_condor(self):
        """Start Jobs with Condor"""
        pass


    def _checktype(self,parameter):
        """Check data type of parameter in neo2in file"""

        if parameter in self.neo2nml:
            if isinstance(self.neo2nml[parameter],str):
                return str
            elif isinstance(self.neo2nml[parameter],list):
                return list
            elif isinstance(self.neo2nml[parameter],int):
                return int
            elif isinstance(self.neo2nml[parameter],float):
                return float
            else:
                raise AttributeError('Did not detect datatype of scanparameter')

    def _set_folder_names(self,overwrite=False):

        if self.scanparameter not in self.structure.split('/'):
            raise AttributeError('Scanparameter is not in Structure')

        for i in self.scanvalues:
            self.neo2nml[self.scanparameter]=i
            self.singlerunpath=os.path.join(self.wdir,self._set_folder_names_4_singlerun())
            self._createfiles(overwrite=overwrite)


    def _set_folder_names_4_singlerun(self):
        structure_list=self.structure.split('/')
        for i, j in enumerate(structure_list):
            if j in self.neo2nml:
                if self._checktype(j)==type(float):
                    structure_list[i]='{0}={1:{2}}'.format(j,self.neo2nml[j],self.float_format)
                elif self._checktype(j)==type(int):
                    structure_list[i]='{0}={1}'.format(j,self.neo2nml[j])
                else:
                    print(j, ' is not a int or float')

        return '/'.join(structure_list)





class SingleRun(Neo2_common_objects):
    """Class for inputs with only one neo2.in configuration
    
    So far Multispecies only,
    in Future there should be an MetaClass above these.
    Also some methods should then upmoved.


    Attributes
    ----------

    singlerunpath : str
        Path of the working directory
    req_files_names : dict
        dict of required files to to run neo2
    neo_nml : Neo2File
        is an instance of the Neo2File class with all variables of neo2.in

    Methods
    -------
    RunInitRuns()
        starting the internal neo2 init run for multispecies
        See flag lsw_multispecies = .true. in neo2.in



    Old Attributes
    --------------
    listofruns : list
        A list of runs for precomputation
    """





    ### Methods should be used but not class in total.
    
    def __init__(self,singlerunpath):


        super().__init__(singlerunpath)
        self.singlerunpath=singlerunpath
        #self.codesource=None
        #self.runsource=None # Executable
    
        #self.req_files_names={
        #'neo2in': 'neo2.in',
        #'neoin': 'neo.in'
       # }
        #self._Load_default()
        #self.req_files_paths=dict()
        #self._SetSources()
        #self._fill_req_files_paths() # including Fill required Files from neo Files, which reads neo2.in
        #self._neo2inSet=False
       # self._neo2parSaved=False
        #self._Runiscreated=False # required Files exists
        #self._Runsettingsaved=False # required Files and Settings saved
        
        ##########Info :self._sources=dict()
       
        #self.SourcePaths=dict()
        #self.SetSourcepaths()

    
    
    #### INIT DONE FOR THE FIRST!!!!####
    
    
    def _Load_default(self):
        """Load default Values stored in input_default.yaml
        
        is in the same directory as the input_class
        """


        try:
            self.__sourcedir=__file__
        except NameError:
            self.__sourcedir=''
        print(self.__sourcedir)
        self.__realpath=os.path.realpath(self.__sourcedir)
        print(self.__realpath)
        if os.path.isfile(self.__realpath):
            path = os.path.dirname(os.path.realpath(self.__sourcedir))
        elif os.path.isdir(self.__realpath):
            path=self.__realpath
    
        yamlpath = os.path.join(path, 'input_default.yaml')
        #print(path)
        try:
          with open(yamlpath) as f:
              data=yaml.load(f)
        except:
            print('Couldn\'t read default values')
            return
        print(data)
        self._runsource=data.get('runsource')
        #self.req_files_names=data.get('reqfilesnames')
        ###DOTO##
        #Cannot load self.req_files_names from YAML file
        self._sources=data.get('sources')

        
    


#    def Run_Precomp(self,overwrite=False):
#        """NOT READY YET!!!!"""
#
#        #load_parameter from an old run!
#        self._checkreqfiles()
#        self._createfiles(overwrite)
#        self.RunInitRuns()
#        #self.gen_condorinput()
#        #self.run_condor()
        
        
        
    def RunInitRuns(self):
        """Use Bash Routine from andfmar"""
        
        
        self._createfiles()
        if self.neo_nml['isw_multispecies_init'] != 1:
            raise AttributeError('isw_multispecies_init is not 1')
        else:
            curdir=os.getcwd()
            print('curdir = ',curdir)
            print('singlerunpath = ',self.singlerunpath)
            os.chdir(self.singlerunpath)
            print('Starting initrun')
            print(os.popen("./neo_2.x").read())
            for a,b,c in os.walk('./'):
                #self.listofruns=b
                break
            os.chdir(curdir)
        
        



        
    def SetNeo2file(self,neo2file): ## Method to set the path of the neo2 file and make all check(Setter of neo2.in)
        """Path to (new) neo2.in File"""
    ###Method to choose different neo2 Files###
        if os.path.isfile(neo2file):
            self.req_files_paths['neo2in']=os.path.abspath(neo2file)
            
        if os.path.isdir(neo2file):
            self.req_files_paths['neoin']=os.path.join(os.path.abspath(neo2file),'neo.in')
            self.req_files_paths['neo2in']=os.path.join(os.path.abspath(neo2file),'neo2.in')
            
            

###### Methods to be implementend: #########

#    def set_pert_file(self):
#        """Method to set Pertubation File"""
#        pass
#
#    def set_axi_sym_file(self):
#        """Method to set axisymetric File"""
#        pass
#
#
#    def set_multispecies_file(self):
#        """Method to set Multispecies File"""
#        pass
#
#
#    def run_local(self):
#        """Start Job on local machine"""
#        pass
#
        
####### Implementend Functions: #########       
        
        
    def _fill_req_files_names(self): # Get required File paths from neo(2).in Files
        """Get additional required File names from neo.in and neo2.in Files"""
        
        
        self._read_neo2in()
        try:
            self.req_files_names['multispec']=self.neo_nml['fname_multispec_in']
        except KeyError:
            print('The neo2.in File has no multispec parameter')
            return 
        
        
        self.req_files_names['in_file_pert']=self.neo_nml['in_file_pert']
        
       
        
        try:
            neo=self.req_files_paths['neoin']
            ## There must be a better catch if neoin is not a File
        except:
            print('neo.in File is missing')
            return
        
        with open(neo,'r') as f:
            for line in f:
                arg = line.split()
                if 'in_file' in arg:
                    self.req_files_names['in_file_axi']=arg[0]
                    #print(arg[0])
                    break 
    # Ensure that it is running in the correct mode
    
    #def _SetSources(self): # name has to be better chosen:
        
        ##### Uniqifi path if only a oberpath is given
       # self._sources=dict([('codesource','/proj/plasma/Neo2/MODULAR/'),
                           #     ('pert_path','/temp/andfmar/Boozer_files_perturbation_field/ASDEX_U/32169/'),                            
                #                ('singlerunsource','/temp/wakatobi/Startordner_Neo2/')])
        
        
    def _fill_req_files_paths(self,overwrite=False): # name has to be better choosen
        """Method for getting full paths from required Files"""


        try:
            files=os.listdir(self._sources['singlerunsource'])
        except:
            print('sources not correct set')
            return
        
        for fdisc,fnames in self.req_files_names.items():
            
            if fnames in files:
                if fnames in self.req_files_paths and overwrite==False:
                    #print(fnames, ' is already set to path: ' self.req_files_paths[fnames])
                    continue
                self.req_files_paths[fdisc]=os.path.join(self._sources['singlerunsource'],fnames)
                
        self._fill_req_files_names()
        #read_Neo2File
        
        
        ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
        if set(self.req_files_names) == set(self.req_files_paths):
            return
        
        try:
            files=os.listdir(self._sources['pert_path'])
        except:
            print('sources2 not correct set')
            return
        
        for fdisc,fnames in self.req_files_names.items():
            
            if fnames in files:
                if fnames in self.req_files_paths:
                    #print(fnames, ' is already set to path: ' self.req_files_paths[fnames])
                    continue
                self.req_files_paths[fdisc]=os.path.join(self._sources['pert_path'],fnames)
        
        self._fill_req_files_paths() # recursive Execution!!! maybe a problem!!!

    
#    def _read_neo2in(self):
#        try:
#            self.neo_nml=Neo2File(self.req_files_paths['neo2in'])
#        except:
#            print('Couldn\'t read neo2.in')
#            return






class Neo2File(object):
    """Class for reading and eding Neo2.in File


    Attributes
    ----------


    Methods
    -------
    chval()
        Change more values at once with keyvalue pairs
    write()
          write changes to file

    """


    def __init__(self,neo2path):
        self.ischanged=False
        self._neo2nml=f90nml.read(neo2path)
        self._neo2dict=dict()
        self._nmltodict()
        self._check_duplicate_parameter()
    
    
    def _nmltodict(self):
        for i in self:
            self._neo2dict[i]=self[i]
    
    def chval(self,**kwargs):
        
        print(self['plotting'])
        for par in kwargs:
            
            match=False
            for i in self._neo2nml: # Iterating through namelists
                if par in self._neo2nml[i]:  # if parameter in one namelist 
                    if self._neo2nml[i][par]==kwargs[par]:
                        print(par, ' in namelist', i, ' was already : ', kwargs[par])
                    else:
                        print('before change: ',par,'=',self._neo2nml[i][par])
                        self._neo2nml[i][par]=kwargs[par]
                        print(par, ' in namelist', i, ' changed to : ', kwargs[par])
                    match=True
                    break
                if match:
                    break
            else:
                print(par, ' is not in the neo2.in File')
    
    
    
    def write(self):
        return self._neo2nml.write()
    
    def _check_duplicate_parameter(self):
        te=dict(self._neo2nml)
        te2=dict(self._neo2nml)
        
        for nml in te:
            del te2[nml]
            for nml2 in te2:
                 for par in te[nml]:
                     if par in te2[nml2]:
                         raise ImportError('This version accepts only one keyword for the whole neo2.in File')
      #print(nml2,'2nd')
        
        
                            
                            
    
    
    def __repr__(self):
        display.display_pretty(self._neo2dict)
    
    
    #def __repr__(self):
     #   lines = ['{']
      #  for key, value in self._neo2dict.items():
       #     lines.append('{}:{}'.format(key, value))
        #lines.append(['}'])
        
       # return '\n'.join(lines)
    
    
    #def __getattr__(self, attr):
      #  if attr=='write':
       #     return getattr(self._neo2nml,attr)

    #def __dir__(self):
     #   print('__dir__ is called')
    #   return ['write']
        #return self._neo2nml.__dir__()
       

    
    #def __setattr__(self, attr, val):
       # print('attr = ',attr)
      #  if attr == '_neo2nml':
       #     object.__setattr__(self, attr, val)

       # return setattr(self._neo2nml, attr, val)
    
    def __getitem__(self, key):
        
        ## Unbedingt Kontrollen einführen  
        if key in self._neo2nml:
            print('Please adress your desired parameter directly')
            return None
        else:
            for i in self._neo2nml:
                if key in self._neo2nml[i]:
                    #print('parameter: ',key,'found in namelist:',i )
                    return self._neo2nml[i][key]
        val = self._neo2nml.__getitem__(key)
        if isinstance(val,f90nml.namelist.Namelist):
            valout='namelist of key'
            
        else:
            valout=val
        print('getter: called: ',key,' with value:', valout)
        return val

    def __setitem__(self, par, val):
               
        if par in self._neo2nml:
            print('Please adress your desired parameter directly')
            return None
        else: 
            for i in self._neo2nml:
                if par in self._neo2nml[i]:
                    if not self._checktype(self._neo2nml[i][par],val):
                        print('wrong type of parameter', par,type(val),'-- old type was:', type(self._neo2nml[i][par]))
                        return
                
                    if self._neo2nml[i][par]==val:
                        print(par, ' in namelist', i, ' was already : ', val)
                        return
                    else:
                        print('before change: ',par,'=',self._neo2nml[i][par])
                        self._neo2nml[i][par]=val
                        print(par, ' in namelist', i, ' changed to : ', val)
                        self.ischanged=True
                        return

        # Auch unbedingt Kotrollen einführen
        #print('setter:',key,'is set to :', val)
        #self._neo2nml.__setitem__(key, val)
        #write File when set.
        
    def _checktype(self,val1,val2):
        if type(val1)==type(val2):
            return True
        else: 
            return False
        

    def __iter__(self):
       # def chain(*iterables):
            # chain('ABC', 'DEF') --> A B C D E F
        for it in self._neo2nml:
            for element in self._neo2nml[it]:
                yield element
    
        


    

    
    #def update(self, *args, **kwargs):
     #   print('updateoe', args, kwargs)
      #  for k, v in f90nml.namelist.Namelist(*args, **kwargs).iteritems():
       #     self[k] = v
    
   # def __contains__(self, key):
    #    #print('contains called for ',key)
    #    return self._neo2nml.__contains__(key)

















