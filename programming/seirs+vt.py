#!/bin/python3
import csv
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import os
from contactmatrix import *

# === Definitions and functions ===

# -- meta --

model_name = "SEIRS"

status = "n"

save_figs = True # False
plot_heatmap = False
plot_model = False
plot_format = "pdf"
fig_dir = "../fig"
cm_dir =  "../contact_matrices/LE" 
bd_dir =  "../befolkningsdata" 
meta_dir =  "../meta" 

# -- wheights --
ws_d = 1
ww_d = 0.5
wh_d = 1
wo_d = 0.2

wp = 1
wc = 0.2



def model_plot(time, data_dict, save_name, titl="", show_model=plot_model):
    plt.figure()
    for key, value in data_dict.items():
        plt.plot(time, value, label=key)
    if titl == "":
        titl = save_name
    plt.title(titl)
    plt.legend(loc="best")
    plt.xlabel('time')
    plt.ylabel('P(t)')
    if save_figs:
        file_name = '{0}_{1}.{2}'.format(model_name, save_name, plot_format)
        save_string = '{0}/{1}'.format(fig_dir, file_name)
        plt.savefig(save_string)
        print("plot '{0}' has been created and saved.".format(file_name))
    if show_model:
        plt.show()

# def heatmap(title, data, show_heatmap=plot_heatmap):
#     plt.figure()
#     plt.title(title)
#     plt.imshow(data, cmap='hot', interpolation='nearest')
#     if save_figs:
#         save_name = 'heatmap_{0}.{1}'.format(title, plot_format)
#         save_string = '{0}/{1}'.format(fig_dir, save_name)
#         plt.savefig(save_string)
#         print("heatmap '{0}' has been created and saved.".format(save_name))
#     if show_heatmap:
#         plt.show()
# 
# def get_headers (filename):
#     f = open(filename, 'r')
#     reader = csv.reader(f)
#     headers = next(reader, None)
#     f.close()
#     return headers
# 
# def get_data(directory, plot_bool=False, nan_bool=True):
#     data = {}
#     for filename in os.listdir(directory):
#         if filename.endswith(".csv"): 
#             fn = filename.split(".")[0]
#             filename_path = os.path.join(directory, filename)
#             temp = np.genfromtxt(filename_path, delimiter=",")
#             if nan_bool:
#                 temp = temp[~np.isnan(temp).all(axis=1)]
#                 temp = temp[:,~np.isnan(temp).all(axis=0)]
#             # sets nan to zero instead of deleting rows/columns
#             # data[fn] = np.nan_to_num(np.genfromtxt(filename_path, delimiter=','))
#             data[fn] = temp
#             if plot_bool:
#                 heatmap(fn, data[fn])
#             # heatmap(fn, data[fn])
#             continue
#     return data
# 
# def pop_binning(contact_matrix, pop_vector):
#     bins = []
#     numbins = len(contact_matrix[0])
#     for i in range(numbins - 1):
#         bins.append(sum(pop_vector[10*i:10*i + 10]))
#     bins.append(sum(pop_vector[10*(numbins-1):]))
#     return bins

# -- wheights --

def gen_w_dict(status):
    if status == "n":
        ws, ww, wh, wo = ws_d, 0, wh_d, wo_d
    if status == "q":
        ws, ww, wh, wo = 0, 0, wh_d, wo_d
    w_dict = {}
    for k in keylist:
        # print("round(w{0}*w{1}, 3)".format(k[-2], k[-1]))
        w_dict[k] = eval("round(w{0}*w{1}, 3)".format(k[-2], k[-1]))
    return w_dict

def gen_w_array_dict(status, n_bins):
    # this section constructs dicts, that contain arrays, 
    # - defining what scenario is put into which agegroup:
    w_array_dict = {}
    for k in keylist:
        w_array_dict[k] = []
        # status list is used to design the different wheight scenarios, defined above, 
        # that each age group is whaigted with.
        status_list = [status, "q"] + [status]*(n_bins - 2)
        for s in status_list:
            w_array_dict[k].append(gen_w_dict(s)[k])
        w_array_dict[k] = np.array(w_array_dict[k])
    return w_array_dict

def gen_w_table(title_list, w_array, savestring):
    outfile = open(savestring, "w")
    structure = "\\begin{{tabular}}{{l | *{{{0}}}{{c}}}}".format(len(title_list))
    print(structure, file=outfile)
    print("\\hline", file=outfile)
    header = "Weigts & "+ " & ".join(title_list) + " \\\\"
    print(header, file=outfile)
    print("\\hline", file=outfile)
    count = 0
    for k,v in w_array.items():
        count += 1
        row = "$w_{{{0}}}$".format(k[-2:]) + " & " + " & ".join([str(i) for i in v]) + " \\\\"
        if count == len(w_array.items()):
            row = row[:-3]
        print(row, file=outfile)
    print("\\end{tabular}", file=outfile)
    outfile.close()

# -- fig caption --
def crit_values(s, e, i, r):
    s_crit = (min(s), time[np.argmin(s)])
    e_crit = (max(e), time[np.argmax(e)])
    i_crit = (max(i), time[np.argmax(i)])
    r_crit = (max(r), time[np.argmax(r)])
    crit_list = [s_crit, e_crit, i_crit, r_crit] 
    return crit_list

def gen_figtext(S, E, I, R, savefile):
    crit_list = crit_values(S, E, I, R)
    fig_text = model_name + " model of COVID-19 with contact mixing. "\
                "Note that {0}\% becomes infected, "\
                "with a maximum infectionratio at {1}\%, "\
                "at time {2} days.".format(round(100*crit_list[3][0],2), 
                                          round(100*crit_list[2][0],2), 
                                          round(crit_list[2][1]))
    outfile = open(savefile, "w")
    print(fig_text, file=outfile)
    outfile.close()

# Force of infection:
def FI(I,params):
    w_array, contact_matrix = params
    wcTI = np.array([np.multiply(w_array[k], np.matmul(contact_matrix[k].T, I)) for k in keylist])
    IM = np.sum(wcTI, axis=0)
    return IM

# function that returns dP/dt
def model(P,t):
    P = np.reshape(P, (len(compartments), n_bins))
    S, E, I, R = P
    IM = FI(I,(w_array, contact_matrix))
    # wcTI = np.array([np.multiply(w_array[k], np.matmul(contact_matrix[k].T, I)) for k in keylist])
    # IM = np.sum(wcTI, axis=0)
    for i in range(n_bins):
        # pop = 1 # pop_bins[i]
        dSdt[i] = -beta*IM[i]*S[i] + phi*R[i] 
        dEdt[i] = beta*IM[i]*S[i] - alpha*E[i]
        dIdt[i] = alpha*E[i] - gamma*I[i] 
        dRdt[i] = gamma*I[i] - phi*R[i]
    dPdt = np.concatenate((dSdt, dEdt, dIdt, dRdt), axis=None)
    return dPdt

def solve():
    # time points
    time = np.linspace(0,int(1e3), int(1e3))
    # solve ODE by integration:
    P = odeint(model,P0,time)

# ========== This is only executed if called as main script ==========

if __name__ == "__main__":
    
    # === Data import ===
    
    # -- Contact matrices --
    headers = get_headers("{0}/LE_ph.csv".format(cm_dir))
    title_list = headers[1:]
    
    contact_data = get_data("{0}".format(cm_dir))
    # massage:
    contact_matrix = {}
    probabilities = {}
    keylist = []
    for key,value in contact_data.items(): 
        keylist.append(key)
        contact_matrix[key] = contact_data[key][0:-1]
        probabilities[key] = contact_data[key][-1]

    prob = probabilities["LE_ph"]
    n_bins = len(prob)
    N = 1 # 3*1e6
    pop_bins = N*prob

    # -- weights and table from data --
    w_array = gen_w_array_dict(status, n_bins)
    # print(w_array)
    
    savestring_w_table = "{0}/w_table_{1}.txt".format(meta_dir, status)
    gen_w_table(title_list, w_array, savestring_w_table)
    
    
    # === modelling ===
    
    compartments = ["S", "E", "I", "R"]
    
    # -- parameters --
    # N = sum(pop_bins)
    i0 = 1e-5 #initial infection
    r0 = 1.2 # 2.68 #basic reproduction number
    alpha = 1/3. #reciprocal incubation time.
    gamma = 1/4. #reciprocal recovery time 
    beta = gamma*r0# 1/3. #reciprocal time between contacts
    phi = 1/(7*50)
    # eta = 0.26 # infection efficiency
    
    
    # -- initializasion --
    # initial condition
    I0 = np.ones(n_bins)*i0*1/pop_bins
    S0 = - I0 + 1 # +pop_bins
    E0 = np.zeros(n_bins)
    R0 = np.zeros(n_bins)
    P0 = np.concatenate((S0, E0, I0, R0),axis=None)
    # print("P0: " + str(P0))
    
    #initializasion of variable vectors
    for c in compartments:
        init_compart = "{} = np.zeros(n_bins)".format(c)
        init_compart_deriv = "d{}dt = np.zeros(n_bins)".format(c)
        exec(init_compart)
        exec(init_compart_deriv)
    IM = np.zeros(n_bins)
    dPdt = np.zeros(n_bins*len(compartments))
    
    # -- modellling --
    solve(model, )    

    
    
    # === Plotting ===
    
    # plot results: individual dynamics of the different age-groups
    for i in range(n_bins):
        for ci in range(len(compartments)):
            exec("{0} = P[:, i + {1}*n_bins]".format(compartments[ci], ci))
        state_dict = {"Susceptible":S, "Exposed":E, "Infected":I, "Removed":R }
        model_plot(time, state_dict, "{0}_{1}".format(title_list[i], status), 
                    "{0} for ages: {1} years".format(model_name, title_list[i]))
        # generating figure text:
        save_file = "{0}/{1}_{2}_figtext_{3}.txt".format(fig_dir, model_name, title_list[i], status)
        gen_figtext(S, E, I, R, save_file)
    
    # plot results: dynamics of the total population
    for ci in range(len(compartments)):
        exec("{0} = np.sum(np.multiply(pop_bins,P[:, {1}*n_bins:(1 + {1})*n_bins]), axis=1)".format(compartments[ci], ci))
    state_dict = {"Susceptible":S, "Exposed":E, "Infected":I, "Removed":R}
    model_plot(time, state_dict, "total_population_mix_{0}".format(status), 
                "{} for total population with mixing".format(model_name))
    # generating figure text:
    save_file = "{0}/{1}_total_population_mix_figtext_{2}.txt".format(fig_dir, model_name, status)
    gen_figtext(S, E, I, R, save_file)
    
