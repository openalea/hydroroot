#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 10:45:01 2019

@author: F. Bauget

Parameters class:
    :init: set the default values

    :Methods:
        - read_file: read the Yaml configuration file containing the model parameters into a Class.
        - init_calculation: perform some initialization
        - parameters_to_list: convert some parameters to list if they are not
"""

# F. Bauget 2021-05-04: added solute parameters

import numpy
import yaml
import pandas
import glob


class Parameters():
    """Init Parameters class. Set setting default values and structure

    :param archi: all inputs related to the architecture generation or reconstruction
    :param read_architecture: Bool
    :param input_dir: string
    :param input_file: list of string
    :param seed: integer or list of integer
    :param length_file: list of strings
    :param length_data: float data containing in the files above
    :param primary_length: float or list of float
    :param branching_delay: float or list of float
    :param branching_variability: float
    :param order_max: integer
    :param segment_length: float
    :param nude_length: float or list of float
    :param ref_radius: float
    :param order_decrease_factor: float
    :param hydro: inputs related to the hydrodynamics
    :param k0: float
    :param axial_conductance_data: list of 2 list of float
    :param solute: inputs related to solutes transport eather permating or not
    :param J_s: float
    :param P_s: float
    :param Cse: float
    :param Ce: float
    :param sigma: float
    :param exp: inputs related to the experimental conditions and measurements
    :param jv: float
    :param psi_e: float
    :param psi_base: float
    :param output: inputs related to the simulations and its output
    :param radfold: float or list of float
    :param axfold: float or list of float
    :param intercepts: list of float
    :param run_nb: int

    """
    def __init__(self):
        # default values
        self.archi = {
            'read_architecture': False,
            'input_dir': '',
            'input_file': [],
            'seed': None,
            'length_file': ['data/length*order2*.csv', 'data/length*order2*.csv'],
            'length_data': [],
            'primary_length': 0.13,
            'branching_delay': 2e-3,
            'branching_variability': 0.25,
            'order_max': 4,
            'segment_length': 1.e-4,
            'nude_length': 0.021,
            'ref_radius': 7.0e-5,
            'order_decrease_factor': 0.7}

        self.hydro = {
            'k0': 92.0,
            'axial_conductance_data': [[0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18],
                                    [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4, 292.6e-4, 262.5e-4, 511.1e-4]]}
        self.solute = {
            'J_s': 0.,
            'P_s': 0.,
            'Cse': 0.,
            'Ce': 0.,
            'Sigma': 1.}
        self.exp = {
            'Jv': 0.1,
            'psi_e': 0.4,
            'psi_base': 0.101325}

        self.output = {
            'radfold': 1.0,
            'axfold': 1.0,
            'intercepts': [0.01, 0.02, 0.03, 0.045, 0.06, 0.08],
            'run_nb': 1}

    def read_file(self, filename = None):
        """Read the input yaml file, set the class variables and perform some initialization
        see :func:`~init_parameter.Parameters.init_calculation`

        :param filename: string (Default value = None)
        :param return: itself
        :param Some: instance variables are converted to list if not they are not list in the yaml file
        :param branching_delay: nude_length

        """

        # read the file as a dictionary
        f = open(filename)

        parameter = yaml.load(f.read(), Loader = yaml.FullLoader)

        # pass the parameters to the class
        for pid in parameter['archi']:
            self.archi[pid] = parameter['archi'][pid]
        for pid in parameter['hydro']:
            self.hydro[pid] = parameter['hydro'][pid]
        if 'solute' in list(parameter.keys()): # because not in the actual version
            for pid in parameter['solute']:
                self.solute[pid] = parameter['solute'][pid]
        for pid in parameter['experimental']:
            self.exp[pid] = parameter['experimental'][pid]
        for pid in parameter['output']:
            self.output[pid] = parameter['output'][pid]

        # transform the parameters as list if needed
        # at this stage allow to launch a set of simulations for different set
        # see def parameters_to_list
        self.archi['primary_length'] = self.parameters_to_list(self.archi['primary_length'])
        self.archi['seed'] = self.parameters_to_list(self.archi['seed'])
        self.archi['branching_delay'] = self.parameters_to_list(self.archi['branching_delay'])
        self.archi['nude_length'] = self.parameters_to_list(self.archi['nude_length'])

        self.output['intercepts'] = self.parameters_to_list(self.output['intercepts'])
        self.output['radfold'] = self.parameters_to_list(self.output['radfold'])
        self.output['axfold'] = self.parameters_to_list(self.output['axfold'])

        # few initializations
        self.init_calculation()

    def init_calculation(self):
        """Set self.archi['length_data'] by reading the two files self.archi['length_file']
        Set the seed to None if seed is not an integer nor a list of integer
        
        :return: itself


        """

        # set the data used to calculate the length laws
        # column names in the length law files
        col_names = ('LR_length_mm', 'relative_distance_to_tip')

        for f in self.archi['length_file']:
            d_path = glob.glob(f)[0]
            pd = pandas.read_csv(d_path, sep = ';', header = 1, names = col_names)
            pd.sort_values(by = 'relative_distance_to_tip', inplace = True)

            self.archi['length_data'].append(pd)


        if type(self.archi['seed']) == str:
            # F. Bauget 2020-04-06 : case where seeds are given in a file
            d_path = glob.glob(self.archi['seed'])[0]
            self.archi['seed'] = []
            lineList = [line.rstrip('\n') for line in open(d_path)]
            for x in lineList:
                try:
                    value = int(x)
                    self.archi['seed'].append(value)
                except ValueError:
                    break
        elif type(self.archi['seed']) != int and type(self.archi['seed']) != list:
            #set seed to None if not integer
            self.archi['seed'] = [None]

    def parameters_to_list(self, parameter):
        """transform parameter to a list

        :param parameter: the parameter to transform to a list
        :param return: 
        :param parameter: list
        :param In: the yaml file it is possible to use the following syntaxe
        :param this: syntax and the corresponding list is calculated
        :param to: a list of one single element
        :param For: example
        :param range: 
        :param 0: 02 is converted to
        :param see: func

        """
        if type(parameter) != list:
            if type(parameter) == str:
                if parameter.find('range') == 0:
                    parameter = parameter.replace('range(', '')
                    parameter = parameter.replace(')', '')
                    var = parameter.split(',')
                    parameter = numpy.arange(float(var[0]), float(var[1]), float(var[2])).tolist()
            elif type(parameter) == float or type(parameter) == int :
                parameter = [parameter]

        return parameter
