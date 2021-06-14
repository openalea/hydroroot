#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 10:45:01 2019

@author: F. Bauget

Work in progress

Parameters class:
    init: set the default values
    Methods: read_file read the Yaml configuration file containing the model parameters into a Class.
            Calls other methods to perform some initialization
"""

import numpy
import yaml
import pandas
import glob


class Parameters():
    """
    Init Parameters class. Set setting default values and structure
    Parameters contains 4 dictionaries grouping inputs in 4 gategories:
        -   archi: all inputs related to the architecture generation or reconstruction
                *  read_architecture: Bool, True = read the architecture from a data (i.e. reconstructed MTG)
                *  input_dir:  string, input directory of the architecture file
                *  input_file: list of string,  architecture files
                *  seed: integer or list of integer, seed for the Markov generator, if none seed is generated
                *  length_file: list of strings, files name (included relative path) containing data used to calculate 
                            length laws for the laterals generation
                *  length_data: float data containing in the files above
                *  primary_length: float or list of float, length of the primary root
                *  branching_delay: float or list of float, average distance between to branching
                *  branching_variability: float, add variability in the delay
                *  order_max: integer, maximum laterals order 
                *  segment_length: float, length of the vertices
                *  nude_length: float or list of float, length from the tip without any branching
                *  ref_radius: float, radius of the primary root
                *  order_decrease_factor: float, decrease factor apply to the radius to account for its decrease 
                            with lateral order
        - hydro: inputs related to the hydrodynamics
                *  k0: float, the radial conductivity
                *  axial_conductance_data: list of 2 list of float, the axial conductance vs distance from the tip
        - exp: inputs related to the experimental conditions and measurements
                *  jv: float, measured flux at the base
                *  psi_e: float, water potential surrounding the root
                *  psi_base: float, water potential at the base
        - output: inputs related to the simulations and its output
                *  radfold: float or list of float, in factor to k0
                *  axfold: float or list of float, in factor to axial_conductance_data
                *  intercepts: list of float, distance from base to calculate the number of intercepts
                *  run_nb: int, number of simulation with the same set of input
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
        """
        Read the input yaml file and perform some initialization
            see init_calculation
        :param filename: the input yaml file with the parameters
        :return: itself
        """

        # read the file as a dictionary
        f = open(filename)

        parameter = yaml.load(f.read(), Loader = yaml.FullLoader)

        # pass the parameters to the class
        for pid in parameter['archi']:
            self.archi[pid] = parameter['archi'][pid]
        for pid in parameter['hydro']:
            self.hydro[pid] = parameter['hydro'][pid]
        if 'solute' in parameter.keys(): # because not in the actual version
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
        """
        Set self.archi['length_data'] by reading the two files self.archi['length_file']
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
        """
        transform parameter to a list
            *   in the yaml file it is possible to use the following syntaxe "range(start, end, step)"
                then a test checks this syntax and the corresponding list is calculated
            *   if it is a float or integer, then it is transform to a list of one single element
        :param parameter: the parameter to transform to a list
        :return:
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