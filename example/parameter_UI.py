import copy
import sys
import io
import yaml
import glob
import codecs

sys.path.extend(['../src'])

import matplotlib.pyplot as plt
import ipywidgets as widgets
import pandas as pd
import numpy as np
from IPython.display import display, clear_output
# from oawidgets.plantgl import PlantGL
# from openalea.plantgl.algo.view import view
from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.display import mtg_scene
from miscelenous_fct import water_solute_model

### GLOBAL
Flag_archi = False
param = Parameters()
default_param = Parameters()
df_archi = pd.DataFrame()
df_law = [pd.DataFrame(), pd.DataFrame()]
df_data = [None, None]

### output
# plot
out_plot = widgets.Output()

out_stdout = widgets.Output(layout = {
    'width': '100%',
    'height': '160px',
    'border': '1px solid black'})

### for correspondance between widgets description and parameters keys
d = {'L PR:': 'primary_length', 'LR delay:': 'branching_delay', 'variability:': 'branching_variability',
     'L nude:': 'nude_length', 'seg length:': 'segment_length', 'R ref:': 'ref_radius',
     'beta:': 'order_decrease_factor',
     'Js:': 'J_s', 'Ps:': 'P_s', 'Cse:': 'Cse', 'Ce:': 'Ce', 'Sigma:': 'Sigma',
     'Jv:': 'Jv', 'P ext:': 'psi_e', 'P base:': 'psi_base',
     'axfold:': 'axfold', 'radfold:': 'radfold'}

dict_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K_exp', 'axial flux': 'J_out',
             'solute concentration': 'C', 'non-permeant concentration': 'Cpeg'}

def dict_to_Parameter(input):
    output = Parameters()
    # pass the parameters to the class
    for pid in input['archi']:
        output.archi[pid] = input['archi'][pid]
    for pid in input['hydro']:
        output.hydro[pid] = input['hydro'][pid]
    if 'solute' in list(input.keys()):  # because not in the actual version
        for pid in input['solute']:
            output.solute[pid] = input['solute'][pid]
    for pid in input['experimental']:
        output.exp[pid] = input['experimental'][pid]
    for pid in input['output']:
        output.output[pid] = input['output'][pid]

    # transform the parameters as list if needed
    # at this stage allow to launch a set of simulations for different set
    # see def parameters_to_list
    output.archi['primary_length'] = output.parameters_to_list(output.archi['primary_length'])
    output.archi['seed'] = output.parameters_to_list(output.archi['seed'])
    output.archi['branching_delay'] = output.parameters_to_list(output.archi['branching_delay'])
    output.archi['nude_length'] = output.parameters_to_list(output.archi['nude_length'])

    output.output['intercepts'] = output.parameters_to_list(output.output['intercepts'])
    output.output['radfold'] = output.parameters_to_list(output.output['radfold'])
    output.output['axfold'] = output.parameters_to_list(output.output['axfold'])

    return output


class my_UI(object):

    def my_widget_FileUpload(self, description, on_label = False):
        c = widgets.Label(value = 'no file selected')
        w = widgets.FileUpload(
            description = description,
            accept = '',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
            multiple = False  # True to accept multiple files upload else False
        )

        def handle_change(change):
            name = w.value[0]['name']
            c.value = name

        def handle_label_change(change):
            global default_param
            filename = c.value
            f = self.wd_load_yaml.children[1].value[0]['content']
            f = codecs.decode(f, encoding = "utf-8")
            d = yaml.load(f, Loader = yaml.FullLoader)
            default_param = dict_to_Parameter(d)
            self.set_widgets_from_Parameters(default_param)

        if on_label:
            c.observe(handle_label_change, names = 'value')

        w.observe(handle_change, names = 'value')

        return widgets.Box([c, w])

    def __init__(self):
        self.wd_load_yaml = self.my_widget_FileUpload('load parameters', on_label = True)

        ################
        ## TABS ########
        ################

        # **************************
        ### archi
        # **************************
        ### Do we use an existing architecture or do we generate it
        if Flag_archi:
            value = 'use existing architecture'
        else:
            value = 'generate architecture'
        self.wd_archi_radio = widgets.RadioButtons(options = ['use existing architecture', 'generate architecture'],
                                                   value = value, description = 'Which architecture:', disabled = False)

        def handle_radio_change(change):
            global Flag_archi
            Flag_archi = (change.new == 'use existing architecture')

        self.wd_archi_radio.observe(handle_radio_change, names = 'value')

        ### File with architecture
        self.wd_archi_file_name = self.my_widget_FileUpload('archi file')

        self.wd_length_law = []
        self.wd_length_law.append(self.my_widget_FileUpload('length law 0'))
        self.wd_length_law.append(self.my_widget_FileUpload('length law 1'))

        ### generated architecture parameters
        if default_param.archi['seed']:
            seed = default_param.archi['seed'][0]
        else:
            seed = ''
        self.wd_archi_int = []
        self.wd_archi_int.append(widgets.Text(description = 'Seed:', value = str(seed),
                                              tooltip = 'if empty will be generated randomly',
                                              disabled = False, layout = widgets.Layout(width = '200px')))
        self.wd_archi_int.append(widgets.IntText(description = 'order_max:', value = default_param.archi['order_max'],
                                                 tooltip = 'maximum of lateral order',
                                                 disabled = False, layout = widgets.Layout(width = '200px')))
        archi_float = [('L PR', default_param.archi['primary_length'], 'primary root length', False),
                       ('LR delay', default_param.archi['branching_delay'], 'averaged distance between laterals (m)',
                        False),
                       ('variability', default_param.archi['branching_variability'],
                        'variability on LR length and delay [0;1]', False),
                       ('L nude', default_param.archi['nude_length'], 'distance from tip without laterals (m)',
                        False),
                       (
                           'seg length', default_param.archi['segment_length'], 'representative element length (m)',
                           False),
                       ('R ref', default_param.archi['ref_radius'], 'radius of primary root or radius of reference',
                        False),
                       ('beta', default_param.archi['order_decrease_factor'],
                        'factor of radius decrease between root order', False)]
        self.wd_archi_float = [widgets.FloatText(description = n + ':',
                                                 tooltip = t, value = v,
                                                 disabled = f, step = v / 10.0,
                                                 layout = widgets.Layout(width = '200px')) for n, v, t, f in
                               archi_float]
        wgridarchi = widgets.GridBox(self.wd_archi_int + self.wd_archi_float,
                                     layout = widgets.Layout(grid_template_columns = "repeat(3, 200px)"))

        wd = widgets.HBox([self.wd_archi_file_name, self.wd_archi_radio])
        vboxfile = widgets.Box([wd, widgets.VBox(self.wd_length_law)])
        wboxarchi = widgets.VBox([vboxfile, wgridarchi])

        # **************************
        ### HYDRO
        # **************************
        self.wd_caption_K = widgets.Label(value = 'axial conductance (K) versus distance to tip (x)')
        # Kx_default ='0 4.84e-05  \n0.0557 4.32e-04 \n0.108 5.63e-02 \n0.1489 6.54e-02 \n0.1803 7.55e-02 \n0.2204 1.07e-01 \n0.235 1.20e-01 \n0.2895 1.43e-01'
        dK = pd.DataFrame(default_param.hydro['axial_conductance_data'])
        dK = dK.transpose()
        Kx_default = dK.to_string(index = False, header = False)
        self.wd_K_float = widgets.Textarea(description = 'K:', tooltip = '1st col: x (m), 2nd col: K (10-9 m4/Mpa/s)',
                                           value = Kx_default, layout = widgets.Layout(width = '400px'))
        self.wd_K = widgets.VBox([self.wd_caption_K, self.wd_K_float])
        self.wd_k0 = widgets.FloatText(description = 'k0:',
                                       tooltip = 'radial conductivity (10-9 m/MPa/s)',
                                       value = default_param.hydro['k0'],
                                       step = 0.1,
                                       layout = widgets.Layout(width = '200px'))
        self.wd_hydro = widgets.HBox([self.wd_K, self.wd_k0])

        # **************************
        ### SOLUTE
        # **************************
        solute_float = [
            ('Js:', default_param.solute['J_s'], 'active pumping rate (mol/m2/s)'),
            ('Ps:', default_param.solute['P_s'], 'permeability coefficient (m/s)'),
            ('Cse:', default_param.solute['Cse'], 'concentration of permeating solutes (mol/m3 or mM)'),
            ('Ce:', default_param.solute['Ce'], 'concentration of non-permeating solutes (mol/m3 or mM)'),
            ('Sigma:', default_param.solute['Sigma'], 'reflexion coefficient')
        ]
        self.wd_solute_float = [widgets.FloatText(description = n,
                                                  tooltip = t, value = v,
                                                  step = v / 10.0,
                                                  layout = widgets.Layout(width = '200px')) for n, v, t in solute_float]

        self.wd_solute = widgets.GridBox(self.wd_solute_float,
                                         layout = widgets.Layout(grid_template_columns = "repeat(2, 300px)"))

        # **************************
        ### EXPERIMENTAL
        # **************************
        exp_float = [('Jv:', default_param.exp['Jv'], 'experimental sap flux (uL/s)'),
                     ('P ext:', default_param.exp['psi_e'], 'hydrostatic pressure in the pressure chamber (MPa)'),
                     ('P base:', default_param.exp['psi_base'], 'pressure at the root base (MPa)')
                     ]
        self.wd_exp_float = [widgets.FloatText(description = n,
                                               tooltip = t, value = v,
                                               step = v / 10.0,
                                               layout = widgets.Layout(width = '200px')) for n, v, t in exp_float]
        self.wd_data_to_fit = []
        self.wd_data_to_fit.append(self.my_widget_FileUpload('cut and flow data'))
        self.wd_data_to_fit.append(self.my_widget_FileUpload('JvP data'))
        self.wd_exp = widgets.Box([widgets.VBox(self.wd_exp_float), widgets.VBox(self.wd_data_to_fit)])

        # **************************
        ### OUTPUT
        # **************************
        output_float = [('axfold:', default_param.output['axfold'], 'factor to explore a K range'),
                        ('radfold:', default_param.output['radfold'], 'factor to explore a k0 range')
                        ]
        self.wd_output_float = [widgets.FloatText(description = n,
                                                  tooltip = t, value = v,
                                                  step = v / 10.0,
                                                  layout = widgets.Layout(width = '200px')) for n, v, t in output_float]
        self.wd_output = widgets.Box(self.wd_output_float)

        # **************************
        ### FNAL TABS
        # **************************
        tab_contents = ['archi', 'hydro', 'solute', 'experimental', 'output']
        tab_children = [wboxarchi, self.wd_hydro, self.wd_solute, self.wd_exp, self.wd_output]
        self.tab = widgets.Tab()
        self.tab.children = tab_children
        [self.tab.set_title(i, title) for i, title in enumerate(tab_contents)]

        # **************************
        ### sim result variables
        # **************************
        self.dresults = None
        self.g = None

        # **************************
        #### RUNNING SIMULATION ###
        # **************************

        self.b_sim = widgets.Button(
            description = 'run',
            disabled = False,
            button_style = '',  # 'success', 'info', 'warning', 'danger' or ''
            tooltip = 'run a direct simulation or an adjustment',
            icon = 'check'
        )

        def on_b_sim_clicked(b):
            self.sim_run()

        self.b_sim.on_click(on_b_sim_clicked)

        style = {'description_width': 'initial'}
        self.wd_opt_check = []
        list_description = ['K', 'k', 'Js', 'Ps', 'sigma']
        self.wd_opt_label = widgets.Label(value = 'parameters to adjust')
        for i in range(5):
            self.wd_opt_check.append(widgets.Checkbox(
                value = False,
                description = list_description[i],
                disabled = False,
                indent = False, layout = widgets.Layout(width = '75px', height = '50px')
            ))
        wd_opt_line_check = widgets.HBox(self.wd_opt_check)
        wd_opt_line_1 = widgets.VBox([self.wd_opt_label, wd_opt_line_check])

        self.wd_opt_verbose = widgets.Checkbox(
            value = False,
            description = 'verbose',
            disabled = False,
            indent = False, layout = widgets.Layout(width = '75px', height = '50px')
        )
        self.wd_data = widgets.Dropdown(
            options = ['all', 'JvP', 'cnf'],
            value = 'all',
            description = 'data to fit:',
            disabled = False,
        )
        self.wd_output_file = widgets.Text(
            value = '',
            placeholder = 'relative or full path',
            description = 'Output file:',
            disabled = False, style = style
        )
        self.wd_graph_clear = widgets.Checkbox(value = True, description = 'clear the output graph or not',
                                               disabled = False, indent = False,
                                               layout = widgets.Layout(width = '75px', height = '50px'))
        wd_opt_line_2 = widgets.HBox([self.wd_opt_verbose, self.wd_data, self.wd_output_file])
        wd_opt_line_3 = widgets.HBox([self.b_sim, self.wd_graph_clear])
        self.wd_opt = widgets.VBox([wd_opt_line_1, wd_opt_line_2, wd_opt_line_3],
                                   layout = widgets.Layout(border = 'solid 1px'))

        self.outBox = widgets.VBox([out_stdout, out_plot])

    def display(self):
        display(self.wd_load_yaml)
        display(self.tab)
        display(self.wd_opt)
        display(self.outBox)

    @out_stdout.capture()
    def sim_run(self):
        global param, df_archi, df_data, fig, axs
        _df_archi = None
        param = self.set_parameters(param)

        param_to_opt = []
        for w in self.wd_opt_check:
            if w.value: param_to_opt.append(w.description)
        if len(param_to_opt) == 0.0: param_to_opt = None
        output = None
        if self.wd_output_file.value: output = self.wd_output_file.value
        if Flag_archi: _df_archi = df_archi

        with out_stdout:
            clear_output(wait = True)
        
        self.dresults, self.g = water_solute_model(parameter = param, df_archi = _df_archi, df_law = None,
                                                   df_cnf = df_data[0], df_JvP = df_data[1],
                                                   Data_to_Optim = param_to_opt,
                                                   Flag_verbose = self.wd_opt_verbose.value,
                                                   data_to_use = self.wd_data.value, output = output,
                                                   optim_method = 'COBYLA', Flag_debug = False,
                                                   Flag_radius = True, Flag_Constraint = True, dK_constraint = -3.0e-2,
                                                   Flag_w_Lpr = False,
                                                   Flag_w_cnf = False)
        # some plots
        ############
        _dresults = self.dresults.replace(r'^s*$', float('NaN'), regex = True)
        with out_plot:
            clear_output(wait = True)
            if self.wd_graph_clear.value and not out_plot.outputs:
                fig, axs = plt.subplots(nrows = 1, ncols = 2, constrained_layout = True)
                axs[0].set_xlabel('DP')
                axs[0].set_ylabel('Jv(P)')
                axs[1].set_xlabel('max_length')
                axs[1].set_ylabel('J cnf')
            elif self.wd_graph_clear.value:
                for i in range(2): axs[i].clear()

            if 'Jexp(P)' in list(_dresults.columns):
                d = _dresults[['dp', 'Jexp(P)', 'Jv(P)']].dropna()
                d.sort_values(['dp'], inplace = True)
                d.plot.scatter('dp', 'Jexp(P)', ax = axs[0])
                d.plot.line('dp', 'Jv(P)', ax = axs[0])
                j = np.array(d.loc[:, ['Jv(P)', 'Jexp(P)']])
                axs[0].set_ylim(j.min(), j.max())
                axs[0].set_xlabel('DP')
                axs[0].set_ylabel('Jv(P)')

            if 'Jexp cnf (uL/s)' in list(_dresults.columns):
                d = _dresults[['max_length', 'Jexp cnf (uL/s)', 'Jv cnf (uL/s)']].dropna()
                d.plot.scatter('max_length', 'Jexp cnf (uL/s)', ax = axs[1])
                d.plot.line('max_length', 'Jv cnf (uL/s)', ax = axs[1])
                j = np.array(d.loc[:, ['Jv cnf (uL/s)', 'Jexp cnf (uL/s)']])
                axs[1].set_ylim(j.min(), j.max())
                axs[1].set_xlabel('max_length')
                axs[1].set_ylabel('J cnf')

            for i in range(1): display(axs[i].figure)

    def set_parameters(self, parameter = Parameters()):
        global df_archi, df_data, df_law, default_param

        #### archi
        # **************
        parameter.archi['input_dir'] = default_param.archi['input_dir']
        parameter.archi['input_file'] = default_param.archi['input_file']
        parameter.archi['read_architecture'] = Flag_archi
        if self.wd_archi_file_name.children[1].value:
            name = self.wd_archi_file_name.children[1].value[0].name
            parameter.archi['input_file'] = [name]  # must be a list of string
            df_archi = pd.read_csv(io.BytesIO(self.wd_archi_file_name.children[1].value[0].content), delimiter = '\t',
                                   dtype = {'order': str})
            df_archi['db'] = df_archi['distance_from_base_(mm)'] * 1.e-3
            df_archi['lr'] = df_archi['lateral_root_length_(mm)'] * 1.e-3
            if 'averaged_diameter_(mm)' in df_archi:
                df_archi['radius'] = df_archi['averaged_diameter_(mm)'] * 0.5e-3
        if len(self.wd_archi_int[0].value) > 0:
            seed = int(self.wd_archi_int[0].value)
        else:
            seed = None
        parameter.archi['seed'] = seed
        parameter.archi['order_max'] = self.wd_archi_int[1].value

        for i in range(2):
            if len(list(self.wd_length_law[i].children[1].value)) > 0:
                name = list(self.wd_length_law[i].children[1].value)[0]
                df_law[i] = pd.read_csv(io.BytesIO(self.wd_length_law[i].children[1].value[0]['content']),
                                        delimiter = ';')
        if len(list(self.wd_length_law[0].children[1].value)) > 0:
            if df_law[1].empty:
                df_law[1] = copy.deepcopy(df_law[0])

        for w in self.wd_archi_float:
            parameter.archi[d[w.description]] = w.value

        # overwrite unwanted default value in Parameters
        if df_law[0].empty:
            parameter.archi['length_data'] = None
        else:
            parameter.archi['length_data'] = df_law
        parameter.archi['length_file'] = None

        ### hydro
        # **************

        # K(x) is entered as 1st col x and 2nd col K
        # below transform the string in 2 lists of floats according to the separator
        var = (self.wd_K_float.value).split('\n')
        for sep in ['\t', ' ', ';', ',']:  # try several separator
            if (var[0].find(sep) > 0) and (len(var[0].split(sep)) >= 2):  # >=2 if trailing space
                break
        x = []
        K = []
        for v in var:
            x.append(float(v.split(sep)[0]))
            K.append(float(v.split(sep)[1]))
        parameter.hydro['axial_conductance_data'] = [x, K]
        parameter.hydro['k0'] = self.wd_k0.value

        ### solute
        # **************

        for w in self.wd_solute_float:
            parameter.solute[d[w.description]] = w.value

        ### experimental
        # **************

        for w in self.wd_exp_float:
            parameter.exp[d[w.description]] = w.value

        for i in range(2):
            if len(list(self.wd_data_to_fit[i].children[1].value)) > 0:
                name = list(self.wd_data_to_fit[i].children[1].value)[0]
                df_data[i] = pd.read_csv(io.BytesIO(self.wd_data_to_fit[i].children[1].value[0]['content']), sep = ',',
                                         keep_default_na = True)

        ### output
        # *********
        for w in self.wd_output_float:
            parameter.output[d[w.description]] = w.value

        # some parameters are list by default so convert number to list if necessary
        for key in ['primary_length', 'seed', 'branching_delay', 'nude_length']:
            if type(parameter.archi[key]) != list:
                parameter.archi[key] = [parameter.archi[key]]
        for key in ['radfold', 'axfold']:
            if type(parameter.output[key]) != list:
                parameter.output[key] = [parameter.output[key]]

        return parameter

    def set_widgets_from_Parameters(self, param = None):
        global df_law, df_archi

        Flag_archi = param.archi['read_architecture']
        if Flag_archi:
            value = 'use existing architecture'
        else:
            value = 'generate architecture'
        self.wd_archi_radio.value = value

        archi_f = glob.glob(param.archi['input_dir'] + param.archi['input_file'][0])
        archi_f = archi_f[0]
        if archi_f:
            df_archi = read_archi_data(archi_f) if param.archi['read_architecture'] else None
        self.wd_archi_file_name.children[0].value = param.archi['input_dir'] + param.archi['input_file'][0]

        i = 0
        for f in param.archi['length_file']:
            col_names = ('LR_length_mm', 'relative_distance_to_tip')
            d_path = glob.glob(f)[0]
            df = pd.read_csv(d_path, sep = ';', header = 1, names = col_names)
            df.sort_values(by = 'relative_distance_to_tip', inplace = True)
            df_law[i] = copy.deepcopy(df)
            self.wd_length_law[i].children[0].value = f
            i += 1

        seed = ''
        if param.archi['seed']:
            seed = str(param.archi['seed'][0])

        self.wd_archi_int[0].value = seed
        self.wd_archi_int[1].value = param.archi['order_max']

        for w in self.wd_archi_float:
            if type(param.archi[d[w.description]]) == list:
                p = param.archi[d[w.description]][0]
            else:
                p = param.archi[d[w.description]]
            w.value = p

        dK = pd.DataFrame(param.hydro['axial_conductance_data'])
        dK = dK.transpose()
        Kx_default = dK.to_string(index = False, header = False)
        self.wd_K_float.value = Kx_default
        self.wd_k0.value = param.hydro['k0']

        for w in self.wd_solute_float:
            p = param.solute[d[w.description]]
            w.value = p

        for w in self.wd_exp_float:
            p = param.exp[d[w.description]]
            w.value = p

        for w in self.wd_output_float:
            if type(param.output[d[w.description]]) == list:
                p = param.output[d[w.description]][0]
            else:
                p = param.output[d[w.description]]
            w.value = p

# class plot_architecture_UI(object):

#         def __init__(self, g):
#             self.g = g

#         def display(self):
#             def my_view(cut: int = 0, prop: str = 'j', imgsize: tuple = (800, 800), perspective: bool = True,
#                         zoom: float = 1,
#                         azimuth: float = 0, elevation: float = 0, line_width = 1.0):
#                 # list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}
#                 global dict_prop
#                 key = list(self.g.keys())[cut]
#                 g = self.g[key].copy()

#                 keys = list(g.property('radius').keys())
#                 radius = np.array(list(g.property('radius').values()))
#                 new_radius = radius * line_width
#                 g.properties()['radius'] = dict(list(zip(keys, new_radius)))

#                 s = mtg_scene(g, prop_cmap = dict_prop[prop], has_radius = True)
#                 return view(scene = s, imgsize = imgsize, perspective = perspective, zoom = zoom, azimuth = azimuth,
#                             elevation = elevation)

#             _list = [i for i in range(len(self.g))]
#             _list_prop = list(dict_prop.keys())

#             widgets.interact(my_view, cut = _list, prop = _list_prop, imgsize = widgets.fixed((800, 800)),
#                      perspective = False, zoom = (0.01, 1), azimuth = (-180, 180), elevation = (-90, 90),
#                      line_width = (1, 5))

# class plot_3D_ui(object):

#     def __init__(self, g):
#             self.g = g

#     def display(self):
#         def my_view(cut: int = 0, prop: str = 'j', line_width = 1.0):
#             # list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}
#             global dict_prop

#             key = list(self.g.keys())[cut]
#             g = self.g[key].copy()

#             keys = list(g.property('radius').keys())
#             radius = np.array(list(g.property('radius').values()))
#             new_radius = radius * line_width
#             g.properties()['radius'] = dict(list(zip(keys, new_radius)))

#             s = mtg_scene(g, prop_cmap = dict_prop[prop], has_radius = True)

#             return PlantGL(s)

#         _list = [i for i in range(len(self.g))]
#         _list_prop = list(dict_prop.keys())

#         widgets.interact(my_view, cut = _list, prop = _list_prop, line_width = (1, 5))
