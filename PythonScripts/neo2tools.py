#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 17:21:45 2018

@author: wakatobi
"""
### Template from Ipthon Notebook doc_18_07_05

import os
import sys

import h5py
from IPython import display
from ipywidgets import widgets
import matplotlib.pyplot as plt
import numpy as np
import shutil
import subprocess
import yaml

import f90nml
import neo2post

if __name__ == '__main__':
    print('I am main')

example_files = {'neo-2-par': 'neo2.in.par-full', 'neo-2-ql': 'neo2.in.ql-full'}
code_paths =  {'neo-2-par': 'NEO-2-PAR', 'neo-2-ql': 'NEO-2-QL'}

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

    def __init__(self, wdir, templatepath=None, code_variant:str='neo-2-ql'):


        self._neo2path=os.environ.get('NEO2PATH')
        # If NEO2PATH is not available the Attribute is NONE

        self.wdir=wdir ## runpath ist wahrscheinlich besser
        if templatepath==None:
            self.templatepath=''
        else:
            self.templatepath=templatepath
        if os.path.isdir(self._neo2path):
            self._path2code=self._neo2path
        else:
            self._path2code=''

        self.path2exe=''

        self.variant = code_variant


        self.req_files_names={
        'neo2in': 'neo2.in',
        'neoin': 'neo.in'
        }
        #self._Load_default()
        self.req_files_paths=dict()

## Prefill if Run exists

        if os.path.isdir(wdir):
            if templatepath:
                self._fill_req_files_paths()
                self._compare_nml2file(self.wdir)
            else:
                try:
                    self._fill_req_files_paths(path=wdir)
                except:
                    print('path has no neo2.in')


    def _check_runpath(self):
        os.path.exists(self.wdir)
        pass

    def compile(self, overwrite: bool = False, number_of_processors: int = 1):
        """Wrapper for compiling the code.

        \attention: for now only with hardcoded build directory name.

        input:
        ------
        overwrite: boolean, if true then it is tolerated if the build
          directory already exists, if false, this leads to an
          exception. [False]
        number_of_processors: integer, number of processor make is
          allowed to use. [1]
        """
        curdir=os.getcwd()
        os.chdir(self._path2code)
        exe = os.path.join(self._path2code, code_paths[self.variant], 'Build_auto/')

        if not os.path.exists(exe):
            os.mkdir(exe)
            os.chdir(exe)
        else:
            if overwrite==False:
                # Avoid side effects: get beack to original directory
                # before raising the exception.
                os.chdir(curdir)
                raise FileExistsError('Please remove folder Build_auto')

        os.chdir(exe)
        print(os.popen("cmake -DCMAKE_BUILD_TYPE=RELEASE ..").read())
        print(os.popen("make clean").read())
        print(os.popen("make -j{0}".format(number_of_processors)).read())

        os.chdir(curdir)
        self.path2exe=exe+'neo_2.x'

    def run_local(self,save_out=''):
        """Start Job on local machine"""
        curdir=os.getcwd()
        os.chdir(self.wdir)

        if not save_out:
            with subprocess.Popen("./neo_2.x", stdout=subprocess.PIPE,
                bufsize=1,stderr=subprocess.STDOUT,universal_newlines=True) as p3:
                for line in p3.stdout:
                    print(line, end='')
        else:
            with open(save_out, "w") as outfile:
                subprocess.run('./neo_2.x', stdout=outfile,stderr=outfile)

        os.chdir(curdir)

    def _read_neo2in(self):
        try:
            self.neo2nml=Neo2File(self.req_files_paths['neo2in'], self.variant)
        except:
            print('Couldn\'t read neo2.in')
            return




    def _createfiles(self,singlerunpath='',overwrite=False,link=True):
        """Create and/or link required files and folder into destination"""

        if not singlerunpath:
            singlerunpath=self.wdir


        self._checkreqfiles()
        if not os.path.exists(self.path2exe):
            raise IOError('Path to Executable is missing')


        os.makedirs(singlerunpath,exist_ok=True)

        for i,j in self.req_files_paths.items():

            destfile=os.path.join(singlerunpath,os.path.basename(j))

            if os.path.isfile(destfile):
                if overwrite:
                    print(destfile, ' is a file and is replaced')
                    os.remove(destfile)
                else:
                    print(destfile, ' is not updated') ### TODO Maybe check if it is the same file!!
                    continue
            if i == 'neo2in':
                try:
                    if self.neo2nml.ischanged:
                        self.neo2nml.write(destfile)
                        print('neo2.in was written')
                        continue
                except AttributeError:
                    print('Neo2in was not read')

            if link:
                os.symlink(os.path.realpath(j),destfile)
            else:
                shutil.copy2(j,destfile)



        i2=self.path2exe
        destfile=os.path.join(singlerunpath,os.path.basename(i2))
        if os.path.isfile(destfile):
            if overwrite:
                print(destfile, ' is a file and is replaced')
                os.remove(destfile)
            else:
                print(destfile, ' is not updated') ### Maybe check if it is the same file!!
                return
        if link:
            os.symlink(os.path.realpath(i2),destfile)
        else:
            shutil.copy2(i2,destfile)
        self._Runiscreated=True

    def _compare_nml2file(self,path):
        self._read_neo2in()
        neonml=Neo2File(os.path.join(path,'neo2.in'), self.variant)
        if self.neo2nml==neonml:
            print('loaded neo2.in and neo2.in onn path are the same')
            return
        else:
            print('neo2in are not the same')

    def _checkreqfiles(self):


        if set(self.req_files_names) != set(self.req_files_paths):
            raise ValueError('Not all paths of required Files defined')

        if not os.path.isfile(self.path2exe):
            raise RuntimeError('Executable was not defined')

        for i,j in self.req_files_paths.items():
            if not os.path.exists(j):
                raise RuntimeError(j,' is not a File')

        print('All required Files are existing')


    def _fill_req_files_names(self):
        """Get additional magnetic File name from neo.in"""

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



    def _fill_req_files_paths(self,overwrite=True,path='',rec=True):
        """Method for getting full paths from required Files"""

        #TODO reset option if path is changing
        if not path:
            path=self.templatepath

        try:
            files=os.listdir(path)
        except:
            print('path is not set correctly')
            return

        for file,filename in self.req_files_names.items():

            if filename in files:
                if file in self.req_files_paths and overwrite==False:
                    print(file, ' is already set to path: ', self.req_files_paths[file])
                    continue
        ##TODO find orginal path to files not only links
                self.req_files_paths[file]=os.path.join(path,filename)


        if 'neo_2.x' in files:
            if not self.path2exe:
                self.path2exe=os.path.join(path,'neo_2.x')
            elif overwrite:
                self.path2exe=os.path.join(path,'neo_2.x')

        self._fill_req_files_names()
        #read_Neo2File


        ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
        if set(self.req_files_names) == set(self.req_files_paths):
            return
        elif rec:
            self._fill_req_files_paths(path=path,rec=False,overwrite=overwrite)
        else:
            print('could not fill all required paths')
            print(self.req_files_names)
            print(self.req_files_paths)




class Neo2QL(Neo2_common_objects):
    def __init__(self,wdir,templatepath=None):
        super().__init__(wdir,templatepath, 'neo-2-ql')


    def _fill_req_files_names(self):
        """Get additional required File names from neo.in and neo2.in Files"""

## Reset if file is Changing
        self._read_neo2in()
        try:
            self.req_files_names['multispec']=self.neo2nml['fname_multispec_in']
        except KeyError:
            print('The neo2.in File has no multispec parameter')
            return

        self.req_files_names['in_file_pert']=self.neo2nml['in_file_pert'].rstrip()



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


class Neo2Par(Neo2_common_objects):
    def __init__(self,wdir,templatepath=None):
        super().__init__(wdir,templatepath, 'neo-2-par')


    def run_recon(self):
        pass



class ReconPlot():
    """ Class for plotting data from Neo2 Reconstruction runs

    Attributes
    ----------
    magplot: MagneticsPlot
        instance to plot magnetic fieldline and Points of Interest

    Methods
    -------
    magneticplot(path):
        set a neo2.in file as sample
    g_plot()
        plot of the spitzer function
    interactive_plot():
        Interactiv Plot for showing g along a magnetic fieldline
    """

    def __init__(self,plotdir,rundir='',templatepath=''):

        self._templatepath=templatepath
        self.req_files_names={
                'g_lambdain':'g_vs_lambda.in',
                'spitzerinterface':'spitzerinterface.in',
                'dentf_lorentz':'dentf_lorentz.x'}
        self.req_files_paths=dict(
                dentf_lorentz='/proj/plasma/Neo2/Interface/Spitzer_Interface/Build/'
                'dentf_lorentz.x',
                spitzerinterface='/proj/plasma/Neo2/Interface/Examples/spitzerinterface.in',
                g_lambdain='/proj/plasma/Neo2/Interface/Examples/g_vs_lambda.in')

        if os.path.isdir(rundir):
            self._rundir=rundir
        else:
            raise IOError('rundir must a valid path')

        self._plotdir=os.path.join(rundir,plotdir)
        # if plotdir is absolute, join only outputs plotdir

        self._fill_req_files_paths()
        if os.path.exists(self._plotdir):
            files=os.listdir(self._plotdir)
            for file,filename in self.req_files_names.items():
                if filename not in files:
                    self._createfiles(link=False)
        else:
            self._createfiles(link=False)

    def _fill_req_files_names(self):
        pass # Check of filename inside spitzer.in has to be done.

    def _fill_req_files_paths(self,path='',replace=True):
        """Method for getting full paths from required Files"""

        #TODO reset option if path is changing
        if not path:
            if not self._templatepath:
                path=self._plotdir
            else:
                path=self._templatepath

        try:
            files=os.listdir(path)
        except:
            print('path is not set correctly:')
            print(path)
            return

        for file,filename in self.req_files_names.items():

            if filename in files:
                if file in self.req_files_paths and replace==False:
                    print(file, ' is already set to path: ', self.req_files_paths[file])
                    continue
        ##TODO find orginal path to files not only links
                self.req_files_paths[file]=os.path.join(path,filename)

        ### DOTO Implement method to iterate over all sources, otherwise problems are occuring
        if set(self.req_files_names) == set(self.req_files_paths):
            return
        else:
            print('could not fill all required paths')


    def _createfiles(self,newplotdir='',overwrite=False,link=False):
        """Create and/or link required files and folder into destination"""

        if not newplotdir:
            newplotdir=self._plotdir

        #self._checkreqfiles() #ToDO Make own _checkreqfiles()

        os.makedirs(newplotdir,exist_ok=True)

        for i,j in self.req_files_paths.items():

            destfile=os.path.join(newplotdir,os.path.basename(j))

            if os.path.isfile(destfile):
                if overwrite:
                    print(destfile, ' is a file and is replaced')
                    os.remove(destfile)
                else:
                    print(destfile, ' is not updated') ### TODO Maybe check if it is the same file!!
                    continue
            if link:
                os.symlink(os.path.realpath(j),destfile)
            else:
                shutil.copy2(j,destfile)


    def g_plot(self,*args,tags='',subplots=False,**kwargs):
        """Line plot of density function.

        Allows to plot quantities from file 'g_vs_lambda.h5'.

        Example usage:
        g_plot('g', second_dimension_index='100')

        First entry is the quantity that should be plotted. As this is
        two dimensional, a x-point (second dimension, x-grid) has to be
        selected.
        In this form you would get a line for each point of interest in
        'g_vs_lamda.h5'.

        input:
        ------
        *args: positional arguments, expected to be strings containing
          the name of the quantity to plot, e.g. 'g'. These are passed
          on to the plotting function called.
        tags: string or list of strings, select from which top level
          groups (points of interest) should be plotted, e.g. with
          'tags="p1"' you could select the first point of interest in the
          file. ['']
        subplots: boolean, if true, then each line (for the tags), gets
          is own graph. If false, then all lines are plotted into a
          single figure. [False]
        **kwargs: additional keyword arguments. These are passed on to
          the plotting function called.
        """

        _file=os.path.join(self._plotdir,'g_vs_lambda.h5')

        if not os.path.exists(_file):
            raise IOError(('Please write interested points first'))

        hdf5=h5py.File(_file,'r')
        fig=plt.figure()
        if subplots:
            for i in range(len(args)):
                fig.add_subplot(len(args),1,i+1)
        else:
            fig.add_subplot(111)

        for i in hdf5:
            if i=='version':
                continue
            if tags:
                if isinstance(tags,(list,tuple)):
                    if i not in tags:
                        continue
                else:
                    if i!=tags:
                        continue
            plot=neo2post.Neo2Plot(hdf5[i],def_x='lambda')
            plot.plot(*args,label=i,ax=fig.axes,**kwargs)
        for i in fig.axes:
            i.legend()
        plt.show()
        hdf5.close()


    def magnetic_plot(self,poi=None,write=False):


        if not hasattr(self,'magplot'):
            self.magplot=neo2post.MagneticsPlot(rundir=self._rundir,plotdir=self._plotdir)
        if poi:
            self.magplot.plot_poi(poi)
            if write:
                self.magplot.write_poi(overwrite=True)
        else:
            self.magplot.plot_magnetics()
        plt.show()

    def _plot_write(self,value):
        self.g_plot(value)
        self._run_dentf(show=False)

    def interactive_plot(self):
        if not hasattr(self,'magplot'):
            self.magplot=neo2post.MagneticsPlot(rundir=self._rundir,plotdir=self._plotdir)
        phi_range=widgets.FloatRangeSlider(
            value=[self.magplot.phi_min,self.magplot.phi_max],
            min=self.magplot.phi_min,
            max=self.magplot.phi_max,
            step=0.1,
            continuous_update=False)
        def set_range(phi_ends):
            self.magplot.phi_begin=phi_ends[0]
            self.magplot.phi_end=phi_ends[1]
        h=widgets.interactive(set_range,phi_ends=phi_range)
        f=widgets.FloatSlider(
            min=self.magplot.phi_min,
            max=self.magplot.phi_max,
            step=0.05,
            continuous_update=False)
        a=widgets.interactive(self.magnetic_plot,poi=f,write=widgets.fixed(True))
        d=widgets.interactive(self._plot_write,value=['g','gpa','gtr'])
        c=widgets.FloatText()
        widgets.jslink((c,'value'),(f,'value'))
        f.observe(d.update,'value')
        phi_range.observe(a.update,'value')
        ad=widgets.HBox((a,d))
        display.display(c,ad,h)

    def _run_dentf(self,show=True):
        """Run dentf_lorentz"""
        curdir=os.getcwd()
        os.chdir(self._plotdir)

        if show:
            with subprocess.Popen("./dentf_lorentz.x", stdout=subprocess.PIPE,
                bufsize=1,stderr=subprocess.STDOUT,universal_newlines=True) as p3:
                for line in p3.stdout:
                    print(line, end='')
        else:
             subprocess.Popen("./dentf_lorentz.x")

#        else:
#            #with open(save_out, "w") as outfile:
#             #   subprocess.run('./dentf_lorentz.x', stdout=outfile,stderr=outfile)

        os.chdir(curdir)

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
#            self.neo2nml=Neo2File(self.req_files_paths['neo2in'], self.variant)
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

    def __init__(self,wdir,templatepath=None, variant:str='neo-2-ql'):
        super().__init__(wdir,templatepath, variant)
        self.structure=''
        self.scanparameter=''#boozer_s/
        self.float_format='{:.5f}'#'1.2e'
        self.scanvalues=[]
        self.folder_name='{0}={1:{2}}'
        self.singleruns_names=dict()
        self.singleruns_instances=[]
        self.singlerun_templatepath=None




    def run_condor(self):
        """Start Jobs with Condor"""
        pass


    def _checktype(self,parameter):
        """Old method Check data type of parameter in neo2in file"""

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

    def _set_singleruns_names(self):
        '''Generate singlerun path names for boozer_s

        Should be in Class Radialscan
        '''
        if self.scanparameter=='boozer_s':
                for i in self.scanvalues:
                    dir_name=float2foldername(i,'{:.5f}')
                    self.singleruns_names[dir_name]={self.scanparameter:i}

        else:
            raise NotImplementedError()


    def _get_singleruns_names(self):

        init=False
        if self.singlerun_templatepath:
            print('singlerun_templatepath will be overriden')
        for root, dirs,files in os.walk(self.wdir):
            for i in dirs:

                f,form=foldername2float(i)
                self.scanvalues.append(f)
                #forms.add(form)

                if not init:
                    singlerun_templatepath=os.path.join(root,i,'neo2.in')
                    self.singlerun_templatepath=singlerun_templatepath
                    init=True
                    continue
                path2=os.path.join(root,i,'neo2.in')
                dict_new=compare2nml(singlerun_templatepath,path2, self.variant)
                self.singleruns_names[i]=dict_new

            self.scanparameter='boozer_s'
            self.structure='boozer_s/'
            break

    def _generate_singleruns_instances(self):
        '''Generate a singlerun instance for every singlerun name'''

        if  self.singlerun_templatepath:
            templatepath=self.singlerun_templatepath
        else:
            print('No singlerun_templatepath is set, i use multirun template path')
            templatepath=self.templatepath

        for i in self.singleruns_names:
            self.singleruns_instances.append(SingleRun(os.path.join(self.wdir,i),templatepath))
            self.singleruns_instances[-1]._fill_req_files_paths() ## Should be included in Singlerun init
            self.singleruns_instances[-1].neo2nml.chval(**self.singleruns_names[i])

        pass

    def _generate_singlerun_files(self):
        ''' Generate only the necessarry files for the singleruns'''

        pass


    def _set_folder_names(self,overwrite=False):
        '''Old method for generating SingleRuns'''

        if self.scanparameter not in self.structure.split('/'):
            raise AttributeError('Scanparameter is not in Structure')

        for i in self.scanvalues:
            self.neo2nml[self.scanparameter]=i
            singlerunpath=os.path.join(self.wdir,self._set_folder_names_4_singlerun())
            print(singlerunpath)
            self._createfiles(singlerunpath,overwrite=overwrite)


    def _set_folder_names_4_singlerun(self):
        '''Old method for generating SingleRuns'''
        structure_list=self.structure.split('/')
        for i, j in enumerate(structure_list):
            if j in self.neo2nml:
                if self._checktype(j)==float:
                    structure_list[i]=self.folder_name.format(j,self.neo2nml[j],self.float_format) #TODO individual file names
                elif self._checktype(j)==int:
                    structure_list[i]=self.folder_name.format(j,self.neo2nml[j],'')
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

    def __init__(self,singlerunpath,templatepath=None):


        super().__init__(singlerunpath,templatepath, 'neo-2-ql')
        self.singlerunpath=singlerunpath
        self.codesource=None
        self.runsource=None # Executable

        self.req_files_names={
          'neo2in': 'neo2.in',
          'neoin': 'neo.in'
        }
        self._Load_default()
        self.req_files_paths=dict()
        #self._SetSources()
        self._fill_req_files_paths() # including Fill required Files from neo Files, which reads neo2.in
        self._neo2inSet=False
        self._read_neo2in()
        self._neo2parSaved=False
        self._Runiscreated=False # required Files exists
        self._Runsettingsaved=False # required Files and Settings saved

        ##########Info :self._sources=dict()

        self.SourcePaths=dict()
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


    def Run_Precomp(self,overwrite=False):
        """NOT READY YET!!!!"""

        #load_parameter from an old run!
        self._checkreqfiles()
        self._createfiles(overwrite)
        self.RunInitRuns()
        #self.gen_condorinput()
        #self.run_condor()


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
            print(os.popen("mpirun -np {0} ./neo_2.x".format(self.neo_nml[num_spec])).read())
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

    def set_pert_file(self):
        """Method to set Pertubation File"""
        pass

    def set_axi_sym_file(self):
        """Method to set axisymetric File"""
        pass


    def set_multispecies_file(self):
        """Method to set Multispecies File"""
        pass


    def run_local(self):
        """Start Job on local machine"""
        pass

####### Implementend Functions: #########


    def _fill_req_files_paths(self,overwrite=False):
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


    def _read_neo2in(self):
        try:
            self.neo_nml=Neo2File(self.req_files_paths['neo2in'], self.variant)
            self._neo2inSet = True
        except:
            print('Couldn\'t read neo2.in')
            return



class Neo2File(object):
    """Class for reading and editing Neo2.in File


    Attributes
    ----------


    Methods
    -------
    chval()
        Change more values at once with keyvalue pairs
    write()
          write changes to file

    """


    def __init__(self,neo2path, variant: str):
        self.ischanged=False
        self._variant = variant
        self._Load(neo2path)
        self._neo2dict=dict()
        self._nmltodict()
        self._check_duplicate_parameter()


    def _Load(self, filename: str):
      """Load a namelist file, given by name(+path)..

      This function loads a namelist file for storing in this object.

      Currently this is done in a two step process, to have all
      parameters available, and not only those found in the file.

      input:
      ------
      filename: string, filename and path of the file to load.
      """
      from os.path import join
      from os import environ

      from f90nml import read, Parser

      nml = f90nml.read(filename)
      # Parser is needed as f90nml.read can not patch.
      parser = f90nml.Parser()
      neo2path = environ.get('NEO2PATH')
      # this will read in the example file, but overwrite the values with those in nml.
      self._neo2nml = parser.read(join(neo2path, 'DOC', example_files[self._variant]), nml_patch_in=nml)


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
                        self._neo2dict[par]=kwargs[par]
                        self.ischanged=True
                    match=True
                    break
                if match:
                    break
            else:
                print(par, ' is not in the neo2.in File')


    def write(self,path=''):
#Should control if write was successfull and then self.ischanged = False
        if path=='':
            return self._neo2nml.write()
        else:
            return self._neo2nml.write(path)

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
        #return self._neo2dict


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
                        self._neo2dict[par]=val
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



def float2foldername(inp, form_sp):
    return "_".join(["es", (form_sp.format(inp)).replace(".", "p")])


def foldername2float(inp):
    f = float(inp.replace("p", ".").split("_")[-1])
    precission = len(inp.split("p")[-1])
    form_sp=''.join(['{:.',str(precission),'f}'])
    return f, form_sp

def compare2nml(path1,path2, variant: str):
    '''Compare function for two fortran namelist files.

    This uses the namelist element of Neo2File class, thus the names of
    the namelists are not considered.
    Each element in file1 is compared with the element in file2, and
    differences are added to a dictionary, with value from the second
    file.
    Elements that are only in the first file are ignored. Elements that
    are only in the second file are added to the dictionary.

    Example usage:
    neo2tools.compare2nml('surface1/neo2.in',
                          'surface2/neo2.in', variant='neo-2-par')

    Assuming that these files correspond to a radial scan, then usually
    only 'boozer_s' and 'conl_over_mfp' should differ (for gernots
    variant).

    Input:
    ------
    file1, file2: strings, paths+name of the two neo2.in files to compare.
    variant: string, either of the keys of dictionary example_files.
      Determines if it is an input file for gernots or andreas version
      of the code.
    '''

    file1=Neo2File(path1, variant)
    file2=Neo2File(path2, variant)

    cc=dict()
    for i in file1:
        aa=file1._neo2dict[i]
        bb=file2._neo2dict.pop(i,None)
        if bb is None:
            continue
        if aa==bb:
            continue
        else:
            cc[i]=bb
    for i in file2._neo2dict:
        cc[i]=file2._neo2dict[i]
    return cc
