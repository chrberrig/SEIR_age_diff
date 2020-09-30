#!/bin/python3
import csv
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os

# === Definitions and functions ===

# -- meta --

status = "n"

save_figs = True # False
plot_heatmap = False
plot_model = False
plot_format = "pdf"
fig_dir = "../fig"
cm_dir =  "../contact_matrices" 
meta_dir =  "../meta" 

# -- wheights --
ws_d = 1
ww_d = 0.5
wh_d = 1
wo_d = 0.2

wp = 1
wc = 0.2


def heatmap(title, data, show_heatmap=plot_heatmap):
    plt.figure()
    plt.title(title)
    plt.imshow(data, cmap='hot', interpolation='nearest')
    if save_figs:
        save_name = 'heatmap_{0}.{1}'.format(title, plot_format)
        save_string = '{0}/{1}'.format(fig_dir, save_name)
        plt.savefig(save_string)
        print("heatmap '{0}' has been created and saved.".format(save_name))
    if show_heatmap:
        plt.show()

def get_headers (filename):
    f = open(filename, 'r')
    reader = csv.reader(f)
    headers = next(reader, None)
    f.close()
    return headers

def get_data(directory, plot_bool=False, nan_bool=True):
    data = {}
    for filename in os.listdir(directory):
        if filename.endswith(".csv"): 
            fn = filename.split(".")[0]
            filename_path = os.path.join(directory, filename)
            temp = np.genfromtxt(filename_path, delimiter=",")
            if nan_bool:
                temp = temp[~np.isnan(temp).all(axis=1)]
                temp = temp[:,~np.isnan(temp).all(axis=0)]
            # sets nan to zero instead of deleting rows/columns
            # data[fn] = np.nan_to_num(np.genfromtxt(filename_path, delimiter=','))
            data[fn] = temp
            if plot_bool:
                heatmap(fn, data[fn])
            # heatmap(fn, data[fn])
            continue
    return data

# def pop_binning(contact_matrix, pop_vector):
#     bins = []
#     numbins = len(contact_matrix[0])
#     for i in range(numbins - 1):
#         bins.append(sum(pop_vector[10*i:10*i + 10]))
#     bins.append(sum(pop_vector[10*(numbins-1):]))
#     return bins
# 

# ========== This is only executed if called as main script ==========

# if __name__ == "__main__":
#     
#     # === Data import ===
#     
#     # -- Contact matrices --
#     headers = get_headers("{0}/LE_ph.csv".format(cm_dir))
#     title_list = headers[1:]
#     
#     contact_data = get_data("{0}".format(cm_dir))
#     # massage:
#     contact_matrix = {}
#     probabilities = {}
#     keylist = []
#     for key,value in contact_data.items(): 
#         keylist.append(key)
#         contact_matrix[key] = contact_data[key][0:-1]
#         probabilities[key] = contact_data[key][-1]
#         # heatmap(key, contact_matrix[key])
# 
#     prob = probabilities["LE_ph"]
#     n_bins = len(prob)
#     N = 1 # 3*1e6
#     pop_bins = N*prob
# 
#     # -- weights and table from data --
#     w_array = gen_w_array_dict(status, n_bins)
#     print(w_array)
#     
#     savestring_w_table = "{0}/w_table_{1}.txt".format(meta_dir, status)
#     gen_w_table(title_list, w_array, savestring_w_table)
#     
#     
#     # === modelling ===
#     
#     compartments = ["S", "E", "I", "R"]
#     
#     # -- parameters --
#     # N = sum(pop_bins)
#     i0 = 1e-5 #initial infection
#     r0 = 1.2 # 2.68 #basic reproduction number
#     alpha = 1/3. #reciprocal incubation time.
#     gamma = 1/4. #reciprocal recovery time 
#     beta = gamma*r0# 1/3. #reciprocal time between contacts
#     phi = 1/(7*50)
#     # eta = 0.26 # infection efficiency
#     
#     
#     # -- initializasion --
#     # initial condition
#     I0 = np.ones(n_bins)*i0*1/pop_bins
#     S0 = - I0 + 1 # +pop_bins
#     E0 = np.zeros(n_bins)
#     R0 = np.zeros(n_bins)
#     P0 = np.concatenate((S0, E0, I0, R0),axis=None)
#     # print("P0: " + str(P0))
#     
#     #initializasion of variable vectors
#     for c in compartments:
#         init_compart = "{} = np.zeros(n_bins)".format(c)
#         init_compart_deriv = "d{}dt = np.zeros(n_bins)".format(c)
#         exec(init_compart)
#         exec(init_compart_deriv)
#     IM = np.zeros(n_bins)
#     dPdt = np.zeros(n_bins*len(compartments))
#     
#     # -- modellling --
#     solve(model, )    
# 
#     
#     
#     # === Plotting ===
#     
#     # plot results: individual dynamics of the different age-groups
#     for i in range(n_bins):
#         for ci in range(len(compartments)):
#             exec("{0} = P[:, i + {1}*n_bins]".format(compartments[ci], ci))
#         state_dict = {"Susceptible":S, "Exposed":E, "Infected":I, "Removed":R }
#         model_plot(time, state_dict, "{0}_{1}".format(title_list[i], status), 
#                     "{0} for ages: {1} years".format(model_name, title_list[i]))
#         # generating figure text:
#         save_file = "{0}/{1}_{2}_figtext_{3}.txt".format(fig_dir, model_name, title_list[i], status)
#         gen_figtext(S, E, I, R, save_file)
#     
#     # plot results: dynamics of the total population
#     for ci in range(len(compartments)):
#         exec("{0} = np.sum(np.multiply(pop_bins,P[:, {1}*n_bins:(1 + {1})*n_bins]), axis=1)".format(compartments[ci], ci))
#     state_dict = {"Susceptible":S, "Exposed":E, "Infected":I, "Removed":R}
#     model_plot(time, state_dict, "total_population_mix_{0}".format(status), 
#                 "{} for total population with mixing".format(model_name))
#     # generating figure text:
#     save_file = "{0}/{1}_total_population_mix_figtext_{2}.txt".format(fig_dir, model_name, status)
#     gen_figtext(S, E, I, R, save_file)
#     
