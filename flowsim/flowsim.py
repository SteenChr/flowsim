# -*- coding: utf-8 -*-
"""
Created on Tue May 18 08:40:41 2021

@author: au156185
"""

# import os 
import yaml as yml
import pandas as pd
import numpy as np
from math import ceil
from datetime import timedelta, datetime
import re
import copy
#import cj_response_funcs as cj
#import flowsim.cj_response_funcs as cj
import flowsim.cj_response_funcs.cj_response_funcs as cj
# from importlib.machinery import SourceFileLoader

import matplotlib.pyplot as plt
import matplotlib.dates as dates
# import matplotlib.ticker as ticker

############################################################################### 
"""
Taken from https://stackoverflow.com/questions/4628122/how-to-construct-a-timedelta-object-from-a-simple-string
"""

regex = re.compile(r'^((?P<days>[\.\d]+?)d)? *'r'((?P<hours>\d+?)h)? *'\
                   r'((?P<minutes>\d+?)m)? *'r'((?P<seconds>\d+?)s)?', 
                   re.IGNORECASE)
def parse_time(time_str):
    parts = regex.match(time_str)
    if not parts:
        return(timedelta(0))
    parts = parts.groupdict()
    time_params = {}
    for name, param in parts.items():
        if param:
            try:
                time_params[name] = int(param)
            except:
                continue
    tdel = timedelta(**time_params)
    return(tdel)
############################################################################### 
    
def resample(df, rule, sampler, datecol = 'date', valcol = 'val', scalar = 1.):
    
    df.set_index(datecol,inplace=True)
    if sampler == 'sum':
        df = df.resample(rule).sum()
    elif sampler == 'mean':
        df = df.resample(rule).mean()
    elif sampler == 'interpolate':
        df = df.resample(rule).interpolate(method='linear')
    elif sampler == 'nearest':
        df = df.resample(rule).nearest()
    elif sampler == 'uniform':
        df = df.resample(rule).ffill()
        df[valcol] = df[valcol]*scalar
    df.reset_index(inplace=True)
    
    return(df)
############################################################################### 

def read_obs(fo, plotname, obs, **kwargs):
    
    plotobs = True
    try:
        plotobs = obs["plotobs"]
        if not isinstance(plotobs, bool):
            plotobs = False
    except:
        plotobs = False
    try:
        file = obs["file"]
    except:
        print_msg("\n\"file\" not defined for obs in plot = %s"%plotname, fo)
        # print("\n\"file\" not defined for obs in plot = %s"%plotname)
        # fo.write("\n\"file\" not defined for obs in plot = %s"%plotname)
        plotobs = False
    try:
        header = obs["header"]
    except:
        header = None
    if not isinstance(header,int):
        if isinstance(header,str):
            if header.lower() == 'none':
                header = None
            else:
                msg="\n\"header\" can only be int>=0, list of int>=0, or "+\
                    "'None' for obs in plot = %s"
                print_msg(msg%plotname, fo)
                # print(msg%plotname)
                # fo.write(msg%plotname)
                plotobs = False
    elif header < 0:
        msg="\n\"header\" can only be int>=0, list of int>=0, or "+\
            "'None' for obs in plot = %s"
        print_msg(msg%plotname, fo)
        # print(msg%plotname)
        # fo.write(msg%plotname)
        plotobs = False
    try:
        colnames = obs["colnames"]
    except:
        colnames = None
    if header == None and colnames == None:
        fmt="\nNeither \"header\" nor \"colnames\" defined for obs in plot = %s"
        print_msg(fmt%plotname, fo)
        # print(fmt%plotname)
        # fo.write(fmt%plotname)
        plotobs = False
    try:
        datecol = obs["date"]
        if colnames != None:
            if not datecol in colnames:
                fmt = ("\n\"date\" given as %s which is not found"+
                       " in \"colnames\" for obs in plot = %s")
                print_msg(fmt%(datecol,plotname), fo)
                # print(fmt%(datecol,plotname))
                # fo.write(fmt%(datecol,plotname))
                plotobs = False
    except:
        print_msg("\n\"date\" not defined for obs in plot = %s"%plotname, fo)
        # fo.write("\n\"date\" not defined for obs in plot = %s"
                      # %plotname)
        plotobs = False
    try:
        dtformat = obs["dtformat"]
    except:
        print_msg("\n\"dtformat\" not defined for obs in plot = %s"%plotname, 
                  fo)
        # fo.write("\n\"dtformat\" not defined for obs in plot = %s"
                      # %plotname)
        plotobs = False
    try:
        qnam = obs["val"]
        if isinstance(qnam,str):
            qnam = [qnam]
        if isinstance(qnam, list):
            for nam in qnam:
                if not isinstance(nam, str):
                    msg = ("\nAn element specified for \"val\" is not a"+
                           "string for obs in plot = %s")
                    print_msg(msg%(plotname), fo)
                    # fo.write(msg%(plotname))
                    plotobs = False
                elif colnames != None:
                    if not nam in colnames:
                        print_msg(("\n\"val\" (%s) given for obs in plot = %s" + 
                                   " is not found in \"colnames\"")
                                  %(nam,plotname), fo)
                        # print(" is not found in \"colnames\"")
                        # fo.write("\n\"val\" (%s) given for obs in plot = %s"
                        #               %(nam,plotname))
                        # fo.write(" is not found in \"colnames\"")
                        plotobs = False
        else:
            msg = ("\nFor \"val\", specify either a string or a" +
                   "list of strings for obs in plot = %s")
            print_msg(msg%(plotname), fo)
            # fo.write(msg%(plotname))
            plotobs = False
    except:
        print_msg("\n\"name\" not defined for obs in plot = %s"%plotname, fo)
        # fo.write("\n\"name\" not defined for obs in plot = %s"%plotname)
        plotobs = False
    try:
        convfact = obs["convfact"]
        try:
            convfact = float(convfact)
        except:
            print_msg("\n\"convfact\" is not a float for obs in plot = %s"
                      %plotname, fo)
            # fo.write("\n\"convfact\" is not a float for obs in plot = %s"
            #       %plotname)
            plotobs = False
    except:
        print_msg("\n\"convfact\" not defined for obs in plot = %s"%plotname,
                  fo)
        # fo.write("\n\"covfact\" not defined for obs in plot = %s"
        #               %plotname)
        plotobs = False
    try:
        dividewith = obs["dividewith"]
        if not isinstance(dividewith,dict):
            print_msg("\n\"dividewith\" not a dictionary for obs in plot = %s"
                      %plotname, fo)
            # fo.write("\n\"dividewith\" not a dictionary for obs in plot = %s"
            #               %plotname)
            plotobs = False
        else:
            for nam in qnam:
                try:
                    dividewith[nam]
                except:
                    print_msg("\n\"dividewith\" misses entry for '%s' for obs in plot = %s"
                              %(nam,plotname), fo)
                    # fo.write("\n\"dividewith\" misses entry for '%s' for obs in plot = %s"
                    #       %(nam,plotname))
                    plotobs=False
                    continue
    except:
        print_msg("\n\"dividewith\" not defined for obs in plot = %s"
                  %plotname, fo)
        # fo.write("\n\"dividewith\" not defined for obs in plot = %s"
        #               %plotname)
        plotobs = False
    try:
        sep = obs["sep"]
    except:
        sep = ","
    try:
        decimal = obs["decimal"]
    except:
        decimal = "."
    try:
        skiprows = obs["skiprows"]
        if not isinstance(skiprows,int):
            print_msg("\n\"skiprows\" is not an integer for obs in plot = %s"
                      %plotname, fo)
            # fo.write("\n\"skiprows\" is not an integer for obs in plot = %s"
            #       %plotname)
            plotobs = False
        elif skiprows < 0:
            print_msg("\n\"skiprows\" is not >= 0 for obs in plot = %s"%plotname,
                      fo)
            # fo.write("\n\"skiprows\" is not >= 0 for obs in plot = %s"%plotname)
            plotobs = False
    except:
        skiprows = 0
        
    resample = False
    for key in kwargs:
        if key == 'resamp_rule':
            # rule = kwargs[key]
            resample = True
            
    if plotobs:
        try:
#       Read observations (in l/s), resample to daily values
            col = copy.deepcopy(qnam)
            col.append(datecol)
            if colnames == None:
                Q = pd.read_csv(file, header=header, sep=sep, decimal = decimal,
                                skiprows=skiprows, usecols=col)
            else:
                Q = pd.read_csv(file, header=header, sep=sep, decimal = decimal,
                                names=colnames, skiprows=skiprows, usecols=col)
            if datecol != "date":
                Q.rename(columns={datecol: "date"}, inplace=True)
            Q.date = pd.to_datetime(Q.date, format=dtformat)
            for nam in qnam:
                Q[nam] = Q[nam]*convfact/dividewith[nam]
            # Resample to daily values
            if resample:
                difdays = Q["date"].iloc[1] - Q["date"].iloc[0]
                rule = timedelta(days=1)
                if difdays < rule: # downsample
                    Q.set_index('date', inplace=True)
                    Q = Q.resample(rule).mean()
                    Q.reset_index(inplace=True)
                elif difdays > rule: # upsample
                    Q.set_index('date', inplace=True)
                    Q = Q.resample(rule).interpolate(method='linear')
                    Q.reset_index(inplace=True)
        except:
            print_msg("\nCould not open or read input file (%s) for obs in plot = %s correctly."
                      %(file,plotname), fo)
            # fo.write("\nCould not open or read input file (%s) for obs in plot = %s correctly."
            #               %(file,plotname))
            Q = None
            plotobs = False

    return(plotobs, qnam, Q)

    ############################################################################### 

def plot_response(self, plotname, Sim, yax = 'lin', **kwargs):

    plot = self.plot[plotname]
    # Settings for plot of series
    try:
        plotseries = plot["plotseries"]
    except:
        plotseries = list()
    try:
        yax = plot["yaxis"]
        if yax != "log":
            yax = "lin"
    except:
        yax = "lin"
    set_ylim = False
    try:
        ylim = plot["ylim"]
        if isinstance(ylim, list) and len(ylim)==2:
            if isinstance(ylim[0],(float,int)) and isinstance(ylim[0],(float,int)):
                if ylim[0] < ylim[1]:
                    ylim = (ylim[0],ylim[1])
                    set_ylim = True
    except:
        pass
    try:
        ytitle = plot["ytitle"]
    except:
        ytitle = ""
    # Settings for observation to be read and plotted
    try:
        obs = plot["obs"]
        plotobs, qnam, Q = read_obs(self.fo, plotname, obs, **kwargs)
    except:
        plotobs = False

#   Prepare plot consisting of column of subplots
    plt.close()
    # fig = plt.figure()
    # plt.axis("off")
    nsim= len(Sim)
    row = nsim
    col = 1
    num = 0
    fig, axs = plt.subplots(row, col)
    # if len(axs) == 1:

#   Make subplot for each simulation period
    for sim in Sim:
        if nsim > 1:
            ax = axs[num]
            num += 1
        else:
            ax = axs

        period = [sim["date"][0], sim["date"][len(sim)-1]]

        if plotobs:
            Q.set_index("date",inplace=True)
            indx1 = Q.index.get_loc(period[0], method='nearest')
            indx2 = Q.index.get_loc(period[1], method='nearest')
            Q.reset_index(inplace=True)
            df = Q.iloc[indx1:indx2+1]
            resample = False
            for key in kwargs:
                if key == 'resamp_rule':
                    rule = kwargs[key]
                    resample = True
            if resample:
                df.set_index("date",inplace=True)
                df = df.resample(rule).mean()
                df.reset_index(inplace=True)
            for nam in qnam:
                ax.plot(df.date, df[nam], label = nam+"-obs.")

        for s in plotseries:
            try:
                ax.plot(sim.date, sim[s], label = s, linestyle="--")
#                ax.plot(sim.date, sim[s]/self.tscale, label = s, linestyle="--")
            except:
                print_msg("\nIn plot %s, could not recognize %s as a simulated series."
                      %(plotname,s), self.fo)
                # self.fo.write("\nIn plot %s, could not recognize %s as a simulated series."
                #               %(plotname,s))

        locator = dates.AutoDateLocator()#(minticks=3, maxticks=7)
        formatter = dates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
##  If matplot too old to contain ConciseDateFormatter
#        ax.xaxis.set_major_locator(dates.YearLocator())
#        for tick in ax.xaxis.get_major_ticks():
#            tick.tick1line.set_markersize(5)
#            tick.tick2line.set_markersize(0)
#        ax.xaxis.set_minor_locator(dates.MonthLocator(bymonthday=1))

        ax.set_xlim(period)
        ax.set_ylabel(ytitle)
        if yax == 'log':
            ax.set_yscale('log')
        if set_ylim:
            ax.set_ylim(ylim)
        ltitle = str(period[0].year)
        if period[1].year != period[0].year:
            ltitle = ltitle+" - "+str(period[1].year)
        ax.legend(title = ltitle)
    if nsim > 1:
        axs[0].set_title("Plot: "+plotname)
    else:
        axs.set_title("Plot: "+plotname)

    if len(plotseries) > 0:
        plotsim=True
    else:
        plotsim = False
        
    if plotobs and plotsim:
        fname = "Sim_Obs_plot_"+plotname+".png"
    elif plotsim:
        fname = "Sim_plot_"+plotname+".png"
    elif plotobs:
        fname = "Obs_plot_"+plotname+".png"
    if plotsim or plotobs:
        plt.savefig(fname,dpi=600, facecolor='w', edgecolor='w',\
        orientation='landscape', format=None,\
        transparent=False, bbox_inches="tight", pad_inches=0.1,\
        ) 
    
    plt.show()
        
    return()


# ############################################################################### 
        
class model:   
        
    def response_function_parameters(self, f):
        par = {
               "lin_res": tuple(['TC']),
               "lin_res_sens_tc": tuple(['TC']),
               
               "solid_step": ("T", "S"),
               "solid_step_sens_d": ("T", "S"),
               "solid_step_sens_x": ("T", "S"),
               
               "solid_lin":  ("T", "S"),
                
               "solid_rad":  ("T", "S", "C"),
   #            "solid_rad_int":  sinf_solid_radiation_unit_response_function_int,
               "solid_rad_sens_d": ("T", "S", "C"),
               "solid_rad_sens_h": ("T", "S", "C"),
               "solid_rad_sens_x": ("T", "S", "C"),
    
               "solid_prod": ("T", "S"),
               "solid_prod_sens_d": ("T", "S"),
               "solid_prod_sens_k": ("T", "S"),
               "solid_prod_sens_x": ("T", "S"),
                
               "solid_rad_prod": ("T","S","C"),
               "solid_rad_prod_sens_d": ("T","S","C"),
               "solid_rad_prod_sens_k": ("T","S","C"),
               "solid_rad_prod_sens_h": ("T","S","C"),
               "solid_rad_prod_sens_x": ("T","S","C"),
            
               "slab_step": ("T","S","L"),
               "slab_step_sens_d": ("T","S","L"),
               "slab_step_sens_x": ("T","S","L"),
            
               "slab_prod":  ("T", "S", "L"),
               "slab_prod_sens_d":  ("T", "S", "L"),
               "slab_prod_sens_k":  ("T", "S", "L"),
               "slab_prod_sens_x":  ("T", "S", "L"),
                
               "slab_rad": ("T", "S", "C", "L"),
               "slab_rad_sens_d": ("T", "S", "C", "L"),
               "slab_rad_sens_h": ("T", "S", "C", "L"),
               "slab_rad_sens_l": ("T", "S", "C", "L"),
               "slab_rad_sens_x": ("T", "S", "C", "L"),
                
               "slab_rad_prod": ("T", "S", "C", "L"),
               "slab_rad_prod_sens_d": ("T", "S", "C", "L"),
               "slab_rad_prod_sens_h": ("T", "S", "C", "L"),
               "slab_rad_prod_sens_k": ("T", "S", "C", "L"),
               "slab_rad_prod_sens_x": ("T", "S", "C", "L")
               }
        return(par[f])

    def response_function_kwargs(self, f):
        kwargs =  {
               "lin_res": {'aggregate': 'sum', 'steady_state_init' : False},
               "lin_res_sens_tc": {'aggregate': 'sum', 'steady_state_init' : False},
               
               "solid_step": {'aggregate': 'step', 'steady_state_init' : True},
               "solid_step_sens_d": {'aggregate': 'step', 'steady_state_init' : True},
               "solid_step_sens_x": {'aggregate': 'step', 'steady_state_init' : True},
               
               "solid_lin":  ("T", "S"),
                
               "solid_rad":  {'aggregate': 'step', 'steady_state_init' : True},
   #            "solid_rad_int":  {'aggregate': 'step', 'steady_state_init' : False},
               "solid_rad_sens_d": {'aggregate': 'step', 'steady_state_init' : True},
               "solid_rad_sens_h": {'aggregate': 'step', 'steady_state_init' : True},
               "solid_rad_sens_x": {'aggregate': 'step', 'steady_state_init' : True},
    
               "solid_prod": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_prod_sens_d": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_prod_sens_k": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_prod_sens_x": {'aggregate': 'step', 'steady_state_init' : False},
                
               "solid_rad_prod": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_rad_prod_sens_d": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_rad_prod_sens_k": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_rad_prod_sens_h": {'aggregate': 'step', 'steady_state_init' : False},
               "solid_rad_prod_sens_x": {'aggregate': 'step', 'steady_state_init' : False},
            
               "slab_step": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_step_sens_d": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_step_sens_x": {'aggregate': 'step', 'steady_state_init' : True},
            
               "slab_prod":  {'aggregate': 'step', 'steady_state_init' : True},
               "slab_prod_sens_d":  {'aggregate': 'step', 'steady_state_init' : True},
               "slab_prod_sens_k":  {'aggregate': 'step', 'steady_state_init' : True},
               "slab_prod_sens_x":  {'aggregate': 'step', 'steady_state_init' : True},
                
               "slab_rad": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_sens_d": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_sens_h": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_sens_l": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_sens_x": {'aggregate': 'step', 'steady_state_init' : True},
                
               "slab_rad_prod": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_prod_sens_d": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_prod_sens_h": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_prod_sens_k": {'aggregate': 'step', 'steady_state_init' : True},
               "slab_rad_prod_sens_x": {'aggregate': 'step', 'steady_state_init' : True}
               } 
        return(kwargs[f])
    
    def cj_func(self, f):
        cj_func = {
                  "sinf_head_perf": "solid_step",
                  "sinf_head_leak": "solid_rad",
                  "sinf_rech_perf": "solid_prod",
                  "sinf_rech_leak": "solid_rad_prod",
                  "fin_head_perf":  "slab_step",
                  "fin_head_leak":  "slab_rad",
                  "fin_rech_perf":  "slab_prod",
                  "fin_rech_leak":  "slab_rad_prod",
                  }
        return(cj_func[f])
        
    def aquifer_specifications_all_right(self):
        try:
            self.aquifers
        except:
            print_msg("\n\"aquifers\" block missing in %s!"
                      %self.yamlfile, self.fo)
            # self.fo.write("\n\"aquifers\" block missing in flowsim.yaml!")
            return(False)
        no_error = True
        for aqf in self.aquifers:
            aq = self.aquifers[aqf]
            func = aq["func"]
            for f in func:
                try:
                    par = self.response_function_parameters(self.cj_func(f))
                    if not isinstance(par, tuple):
                        par = tuple(par)
                    for p in par:
                        try:
                            aq[p]
                        except:
                            print_msg("\nAquifer %s misses formation about %s"%(aqf,p),
                                      self.fo)
                            # self.fo.write("\nAquifer %s misses formation about %s"%(aqf,p))
                            no_error = False
                    bcname = func[f]        
                    try:
                        self.bcfiles[bcname]#bc[ibc]]
                    except:
                        msg = "\nFor aquifer \"%s\", the bc \"%s\" is not specified "+\
                              "under \"boundaryconditions\" in %s!"
                        print_msg(msg%(aqf,bcname,self.yamlfile), self.fo)#bc[ibc]))
                        # self.fo.write(msg%(aqf,bcname))#bc[ibc]))
                        no_error = False
                except:
                    print_msg(
                        "\nUnknown function \"%s\" specified for aquifer \"%s\""
                        %(f,aqf), self.fo)
                    # self.fo.write("\nUnknown function \"%s\" specified for aquifer \"%s\""%(f,aqf))
                    no_error = False
        return(no_error)

    def bc_read(self):
        resampler = ['sum', 'mean', 'interpolate', 'nearest', 'uniform']
#        to_m = {'m': 1., 'cm': .01, 'mm': .001}
        self.bcdict = dict()
#        self.bc_to_meter = dict()
        no_error = True
        for key in self.bcfiles.keys():
            bc = self.bcfiles[key]
            try:
                file = bc["file"]
            except:
                print_msg("\n\"file\" not defined for bc = %s"%key, self.fo)
                # self.fo.write("\n\"file\" not defined for bc = %s"%key)
                no_error = False
            try:
                header = bc["header"]
            except:
                header = None
            if not isinstance(header,int):
                if isinstance(header,str):
                    if header.lower() == 'none':
                        header = None
                    else:
                        msg="\n\"header\" can only be int>=0, list of int>=0, or 'None'"\
                        + "for bc = %s"
                        print_msg(msg%key, self.fo)
                        # self.fo.write(msg%key)
                        no_error = False
                else:
                    msg="\n\"header\" can only be intt>=0, list of intt>=0, or 'None'"\
                    + "for bc = %s"
                    print_msg(msg%key, self.fo)
                    # self.fo.write(msg%key)
                    no_error = False
            elif header < 0:
                msg="\n\"header\" can only be int>=0, list of int>=0, or 'None'"\
                + "for bc = %s"
                print_msg(msg%key, self.fo)
                # self.fo.write(msg%key)
                no_error = False
            try:
                colnames = bc["colnames"]
            except:
                colnames = None
            if header == None and colnames == None:
                fmt="\nNeither \"header\" nor \"colnames\" defined for bc = %s"
                print_msg(fmt%key, self.fo)
                # self.fo.write(fmt%key)
                no_error = False
            try:
                datecol = bc["date"]
                if colnames != None:
                    if not datecol in colnames:
                        fmt = ("\n\"date\" given as %s which is not found"+
                               " in \"colnames\" for bc = %s")
                        print_msg(fmt%(datecol,key), self.fo)
                        # self.fo.write(fmt%(datecol,key))
                        no_error = False
            except:
                print_msg("\n\"date\" not defined for bc = %s"%key, self.fo)
                # self.fo.write("\n\"date\" not defined for bc = %s"%key)
                no_error = False
            try:
                valcol = bc["val"]
                if colnames != None:
                    if not valcol in colnames:
                        fmt = ("\n\"val\" given as %s which is not found in "+
                               "\"colnames\" for bc = %s")
                        print_msg(fmt%(valcol,key), self.fo)
                        # self.fo.write(fmt%(valcol,key))
                        no_error = False
            except:
                print_msg("\n\"val\" not defined for bc = %s"%key, self.fo)
                # self.fo.write("\n\"val\" not defined for bc = %s"%key)
                no_error = False
            try:
                dtformat = bc["dtformat"]
            except:
                print_msg("\n\"dtformat\" not defined for bc = %s"%key, self.fo)
                # self.fo.write("\n\"dtformat\" not defined for bc = %s"%key)
                no_error = False
            try:
                bctype = bc["type"]
                if not bctype in ('head', 'flux'):
                    print_msg("\n%s not valid as \"type\" bc = %s"
                              %(bctype,key), self.fo)
                    # self.fo.write("\n%s not valid as \"type\" bc = %s"
                    #       %(bctype,key))
                    no_error = False
            except:
                print_msg("\n\"type\" not specified for bc = %s"%key, self.fo)
                # self.fo.write("\n\"type\" not specified for bc = %s"%key)
                no_error = False
            try:
                convfact = bc["convfact"]
                try:
                    to_right_units = float(convfact)
                except:
                    print_msg("\n\"convfact\" is not a float for bc = %s"%key, 
                              self.fo)
                    # self.fo.write("\n\"covfact\" is not a float for bc = %s"%key)
                    no_error = False
            except:
                print_msg("\n\"convfact\" not defined for bc = %s"%key, self.fo)
                # self.fo.write("\n\"covfact\" not defined for bc = %s"%key)
                no_error = False
            try:
                sep = bc["sep"]
            except:
                sep = "'\s+|[,\s*|;\s*'"
                # '\s+|,\s*|;\s*' means: either "one or more white spaces (\s+), 
                # or (|) comma followed by zero or more white spaces (,\s*),
                # or (|) semicolon followed by zero or more white spaces (;\s*)
                sep = "'\s+|[,;]\s*|;\s*'"
            try:
                decimal = bc["decimal"]
            except:
                decimal = "."
            try:
                skiprows = bc["skiprows"]
                if not isinstance(skiprows,int):
                    print_msg("\n\"skiprows\" is not an integer for bc = %s"
                              %key, self.fo)
                    # self.fo.write("\n\"skiprows\" is not an integer for bc = %s"%key)
                    no_error = False
                elif skiprows < 0:
                    print_msg("\n\"skiprows\" is not >= 0 for bc = %s"%key, 
                              self.fo)
                    # self.fo.write("\n\"skiprows\" is not >= 0 for bc = %s"%key)
                    no_error = False
            except:
                skiprows = 0
            if no_error:
                try:
                    col = [datecol, valcol]
                    if colnames == None:
                        bcd = pd.read_csv(file, header=header, usecols=col,
                                          sep = sep, decimal = decimal, 
                                          skiprows=skiprows, engine='python')
                    else:
                        bcd = pd.read_csv(file, header=header, usecols=col, sep=sep, 
                                          decimal = decimal, names=colnames,
                                          skiprows=skiprows, engine='python')
                    bcd.rename(columns = {col[0]: 'date', col[1]: 'val'}, inplace=True)
                    # Change to meter as unit
                    bcd["val"] = bcd["val"]*to_right_units
                    try:
                        bcd.date = pd.to_datetime(bcd.date, format=dtformat)
                        
                        # If needed, resample to steplength
                        timdel = bcd['date'][1] - bcd['date'][0]
                        steplength = self.steplength
                        if timdel < steplength: # Downsample
                            try:
                                if bctype == 'flux':
#                                    downsampler = 'sum'
                                    downsampler = 'mean'
                                else:
                                    downsampler = 'mean'
                                if not downsampler in resampler:
                                    msg = "\n%s not valid \"downsampler\" for bc = %s"%(downsampler,key)
                                    print_msg(msg, self.fo)
                                    # self.fo.write("\n%s not valid \"downsampler\" for bc = %s"%(downsampler,key))
                                    no_error = False
                                else:
                                    bcd = resample(bcd, steplength, downsampler)
                            except:
                                print_msg("\n\"downsampler\" not defined for bc = %s"%key, self.fo)
                                # self.fo.write("\n\"downsampler\" not defined for bc = %s"%key)
                                no_error = False
                        elif timdel > steplength: # Upsample
                            try:
                                if bctype == 'flux':
                                    upsampler = 'uniform'
                                else:
                                    upsampler = 'interpolate'
                                if not upsampler in resampler:
                                    print_msg("\n%s not valid \"upsampler\" for bc = %s"%(upsampler,key), self.fo)
                                    # self.fo.write("\n%s not valid \"upsampler\" for bc = %s"%(upsampler,key))
                                    no_error = False
                                else:
                                    bcd = resample(bcd, steplength, upsampler, 
                                                   scalar = (steplength/timdel))
                            except:
                                print_msg("\n\"downsampler\" not defined for bc = %s"%key, self.fo)
                                # self.fo.write("\n\"downsampler\" not defined for bc = %s"%key)
                                no_error = False
                    except:
                        msg = "\nFor bc \"%s\", could not convert date using format %s!"
                        print_msg(msg%(key, self.dtformat), self.fo)
                        # self.fo.write(msg%(key, self.dtformat))
                        no_error = False                        
####                    # Convert from length (m) to flux (m/d)
#                    if bctype == 'flux':
#                        bcd["val"] = bcd["val"]/(steplength.days +
#                           steplength.seconds/86400.)
                        
                    self.bcdict[key] = bcd
                except:
                    fmt="\nFor bc \"%s\", could not read csv file %s correctly!"
                    print_msg(fmt%(key, file), self.fo)
                    # self.fo.write(fmt%(key, file))
                    no_error = False
                    
        return(no_error)


    def bc_specifications_all_right(self):
        try:
            self.bcfiles
            return(self.bc_read())
        except:
            print_msg("\n\"boundaryconditions\" block missing in %s!"%self.yamlfile, self.fo)
            # self.fo.write("\n\"boundaryconditions\" block missing in flowsim.yaml!")
            return(False)

    def simulation_periods_all_right(self):
        sim_per = list()
        no_error = True
        try:
            for p in self.simulation_periods:
                try:
                    begin_str = p['begin']
                    end_str = p['end']
                except:
                    print_msg("\n'begin' or 'end' not defined for 'simulation_periods'!"
                          , self.fo)
                    # self.fo.write("\n'begin' or 'end' not defined for 'simulation_periods'!")
                    no_error = False
                    continue
                try:
                    begin = datetime.strptime(begin_str,self.dtformat)
                    end = datetime.strptime(end_str,self.dtformat)
                    if not begin < end:
                        print_msg(("\nFor 'simulation_periods', 'begin' value (%s)"+
                              " must be smaller than 'end' value (%s)!")%
                              (begin_str, end_str), self.fo)
                        # self.fo.write("\n'begin' or 'end' not defined for 'simulation_periods'!")
                        no_error = False
                        continue
                    sim_per.append((begin, end))
                except:
                    print_msg("\nPeriod [%s, %s] could not be converted using format %s!"\
                          %(p[0], p[1], self.dtformat), self.fo)
                    # self.fo.write("\nPeriod [%s, %s] could not be converted using format %s!"\
                    #       %(p[0], p[1], self.dtformat))
                    no_error = False
            self.simulation_periods = sim_per
        except:
            print_msg(("\n\"simulation_periods\" block missing in %s!"
                    %self.yamlfile), self.fo)
            # self.fo.write("\n\"simulation_periods\" block missing in flowsim.yaml!")
            no_error = False

        return(no_error)

    def response_setting_all_right(self):
        response = dict([('hydro', ('head', 'flux')), ('heat', ('temp', 'flux'))])
        try:
             r = response[self.responsetype]
        except:
            print_msg("\nWrong \"responsetype\" ('%s') given in %s!"
                  %(self.responsetype, self.yamlfile), self.fo)
            # self.fo.write("\nWrong \"responsetype\" ('%s') given in flowsim.yaml!"
            #               %self.responsetype)
            return(False)
        if not self.response in r:
            print_msg("\nWrong \"response\" ('%s') given for \"responsetype\" ('%s') in %s!"
                  %(self.response, self.responsetype, self.yamlfile), self.fo)
            # self.fo.write("\nWrong \"response\" ('%s') given for \"responsetype\" ('%s') in flowsim.yaml!"
            #       %(self.response, self.responsetype))
            return(False)
        return(True)
        
    def x_all_right(self):
        if not isinstance(self.x, list):
            print_msg("\n\"In %s, x must be a list of positive floats!"
                      %self.yamlfile, self.fo)
            # self.fo.write("\n\"In flowsim.yaml, x must be a list of positive floats!")
            return(False)
        for x in self.x:
            if not (isinstance(x,int) or isinstance(x,float)):
                print_msg("\n\"In %s, a value in x is not float or integer!"
                          %self.yamlfile, self.fo)
                # self.fo.write("\n\"In flowsim.yaml, a value in x is not float or integer!")
                return(False)
            if x < 0.0:
                print_msg("\n\"In %s, a value in x is negative!"%self.yamlfile,
                          self.fo)
                # self.fo.write("\n\"In flowsim.yaml, a value in x negative!")
                return(False)
        return(True)
            
        
    def initialize_read(self):
        self.responsetype = "hydro" # alternative is "heat"
        self.response = "flux" # alternative is "head" (for "head") or "temp" (for "heat")
        self.x = [0.0]
        self.steplength = timedelta(1)
        self.dtformat = "%Y-%m-%d"
        self.warmup_days = 0. # 3652.5
        self.makeplot = False
        self.userdef_output = False
        self.plot = dict()

        no_error = True
        for doc in self.input:
            # yaml_obj = yml.dump(doc, indent=4, default_flow_style=False)
            # f = open("yaml.dump","w")
            # f.write(yaml_obj)
            # f.close()
            for key, value in doc.items():
                if key == "aquifers":
                    self.aquifers = value
                elif key == "boundaryconditions":
                    self.bcfiles = value
                elif key == "simulation_periods":
                    self.simulation_periods = value
                elif key == "dtformat":
                    self.dtformat = value
                elif key == "warmup_days":
                    self.warmup_days = value
                elif key == 'steplength':
                    steplength = parse_time(value)
                    if steplength > timedelta(0):
                        self.steplength = steplength
#                elif key == "responsetype":
#                    self.responsetype = value.lower()
                elif key == "response":
                    self.response = value.lower()
                elif key == 'x':
                    self.x = value
#                elif key == "makeplot":
#                    if isinstance(value,bool):
#                        self.makeplot = value
                elif key == "plot":
                    self.makeplot = True
                elif key == "userdef_output":
                    if isinstance(value,bool):
                        self.userdef_output = value
                elif key == 'plot':
                    self.plot = value
        if not self.bc_specifications_all_right():
            no_error = False
        if not self.aquifer_specifications_all_right():
            no_error = False
        if not self.simulation_periods_all_right():
            no_error = False
        if not self.response_setting_all_right():
            no_error = False
        if not self.x_all_right():
            no_error = False

        return(no_error)
        
    def load_input_file(self, yaml):
        no_error = True

        try:
            self.yamlfile = yaml
            fi = open(self.yamlfile,'r')
            try:
                self.input = yml.load_all(fi, Loader=yml.FullLoader)
            except:
                print_msg("\nCould not load input file %s"%self.yamlfile, 
                          self.fo)
                # self.fo.write("\nCould not load input file %s"%self.yamlfile)
                no_error = False
        except:
            print_msg("\nCould not open input file %s"%self.yamlfile, self.fo)
            # self.fo.write("\nCould not open input file %s"%self.yamlfile)
            no_error = False
            
        return(no_error)
        
    def simulation(self):
        
        steplength = self.steplength
        no_warmup_per = ceil(self.warmup_days/(steplength.days +
                                               steplength.seconds/86400.))
        x = self.x
        nx = len(x)
        nsys = len(self.aquifers)
        
        # Determine beginning and end dates common to all bc time series
        ibc = 0
        for bc in self.bcdict.values():
            end = len(bc)-1
            if ibc < 1:
                bc_date_begin = bc["date"].iloc[0]
                bc_date_end   = bc["date"].iloc[end]
                ibc += 1
            else:
                if bc_date_begin < bc["date"].iloc[0]:
                    bc_date_begin = bc["date"].iloc[0]
                if bc_date_end > bc["date"].iloc[end]:
                    bc_date_end = bc["date"].iloc[end]
            
        # Set extended simulation periods to begin "warmup_days" earlier than
        # beginning of period given in input; also check that bc time series
        # extend over the simulation periods - otherwise adjust.
        simulation_periods = list()
        for sp in self.simulation_periods:
            d1 = sp[0] - steplength*no_warmup_per
            while d1 < bc_date_begin:
                d1 = d1 + steplength
            d2 = sp[1]
            while d2 > bc_date_end:
                d2 = d2 - steplength
            simulation_periods.append((d1, d2))
        
        # Make list of DataFrames for extended simulation periods;
        # resample to simulation stress periods
        Sim = list()
        for period in simulation_periods:
            #Make list of dates:
            d = pd.date_range(period[0], period[1]).tolist()
            Days = np.arange(0.0, float(len(d)), 1.0)
            d = pd.date_range(period[0], period[1], freq=steplength).tolist()
            Days = list()
            for i in range(0,len(d)):
                timdelt = d[i] - d[0]
                Days.append(timdelt.days+timdelt.seconds/86400.)
            df = pd.DataFrame(list(zip(d,Days)),columns=["date","tday"])
            for key, val in self.bcdict.items():
                # Locate which part of dataframe to use and insert to df                
                val.set_index("date",inplace=True)
                indx1 = val.index.get_loc(period[0], method='nearest')
                indx2 = indx1 + len(df)
                val.reset_index(inplace=True)
                df[key] = val["val"].iloc[indx1:indx2].tolist()
            zeros = np.zeros_like(Days)
            for ix in range(0,nx):
                for aqf in self.aquifers:
                    if self.response == 'flux':
                        df["q_"+aqf+"_x="+str(x[ix])] = zeros
                    elif self.responsetype == 'heat':
                        df["T_"+aqf+"_x="+str(x[ix])] = zeros
                    else:
                        df["h_"+aqf+"_x="+str(x[ix])] = zeros
                if nsys > 1:        
                    if self.response == 'flux':            
                        df["q_tot_x="+str(x[ix])] = zeros
                    elif self.responsetype == 'heat':
                        df["T_tot_x="+str(x[ix])] = zeros
                    else:
                        df["h_tot_x="+str(x[ix])] = zeros
            Sim.append(df)
            
        # Make simulation of response at x
        for isim in range(0,len(Sim)):
            sim = Sim[isim]
            if nsys > 1 and self.response == 'flux':
                resp_tot = np.zeros((len(sim), len(x)), dtype=float)
            for aqf in self.aquifers:
                aq = self.aquifers[aqf]
                func = aq["func"]
                resp_aqf = np.zeros((len(sim), len(x)), dtype=float)
                for f in func:
                    # Make list of bc values
                    bcname = func[f]
                    try:
                        bcfac = aq['bcfac'][bcname]
                    except:
                        bcfac = 1.0
                    bc = np.array((bcfac*sim[bcname]).tolist())
                    # cj_response_funcs name
                    cjf = self.cj_func(f)
                    # Make list of function arguments (parameter values)
                    par = list()
                    pnames = cj.response_function_parameters(cjf)
                    for p in pnames:
                        if p == 'D': # Diffusivity
                            D = aq["T"]/aq["S"]
                            par.append(D)
                        elif p == 'K':
                            K = aq["T"]
                            par.append(K)
                        elif p == 'h':
                            h = aq["C"]/aq["T"]
                            par.append(h)
                        elif p == 'l':
                            l = aq["L"]
                            par.append(l)
                        elif p == 'TC':
                            TC = aq["TC"]
                            par.append(TC)
                    par = tuple(par)
                    t = np.array((sim.tday).tolist())
                    kwargs = self.response_function_kwargs(cjf)
                    if kwargs['steady_state_init']:
                        bcinit = np.average(bc)
                        kwargs['bcinit'] = bcinit
                    if self.response == 'flux':
                        if cjf == 'lin_res':
                            resp_aqf += cj.cts_response(cjf, bc, x, t, *par,
                                                        **kwargs)/TC
                        else:
                            fluxfac = 1./aq["L"] # factor to convert flux from m2/d to m/d
        
                            # Simulate flux at bc (at x = 0)
                            dhdx = cj.cts_response(cjf+"_sens_x", bc, x, t, *par,
                                                   **kwargs)
                            T = aq["T"]
                            resp_aqf += fluxfac * T * dhdx#[:,0]
                    else:
                        resp_aqf += cj.cts_response(cjf, bc, x, t, *par, **kwargs)
                        
                if nsys > 1 and self.response == 'flux':
                    resp_tot += resp_aqf
                if self.response == 'flux':
                    for ix in range(0,len(x)):
                        Sim[isim]["q_"+aqf+"_x="+str(x[ix])] = resp_aqf[:,ix]
                elif self.responsetype == 'heat':
                    for ix in range(0,len(x)):
                        Sim[isim]["T_"+aqf+"_x="+str(x[ix])] = resp_aqf[:,ix]
                else:
                    for ix in range(0,len(x)):
                        Sim[isim]["h_"+aqf+"_x="+str(x[ix])] = resp_aqf[:,ix]
            if nsys > 1:
                if self.response == 'flux':            
                    for ix in range(0,len(x)):
                        Sim[isim]["q_tot_x="+str(x[ix])] = resp_tot[:,ix]

        # Drop rows from extended periods to match input simulation periods
        isim = 0
        for period in self.simulation_periods:
            df = Sim[isim]    
            df.set_index("date", inplace=True)
            indx = df.index.get_loc(period[0])
            df.reset_index(inplace=True)
            droplist = list(range(0,indx))
            df.drop(droplist, axis=0, inplace=True)
            df.reset_index(drop=True, inplace=True)
            isim += 1

        return(Sim)

###############################################################################    
def write_pest_read_file(Sim):
    f_pest = open("flowsim-PEST.csv","w")
    simsum = 0.0
    summer_sum = 0.0
    winter_sum = 0.0
    ndata=0
    for sim in Sim:
        sim.to_csv(f_pest, index=False, columns=['q_tot_x=0.0'], 
                   float_format=" %12.6f", line_terminator="\n")
        difdays = (sim["date"].iloc[1] - sim["date"].iloc[0]).days
        # Sum of summer and winter flow; assumes time unit is DAY
        if difdays < 1:
            sim.set_index('date', inplace=True)
            sim1 = sim.resample('1D').mean()
            sim.reset_index(inplace=True)
            sim1.reset_index(inplace=True)
        elif difdays > 1: # The following has not been tested
#            print("\nUpsampling scheme has not been tested!")
            sim.set_index('date', inplace=True)
            sim1 = sim.resample('1D').bfill()
            sim.reset_index(inplace=True)
            sim1.reset_index(inplace=True)
        else:
            sim1 = sim
        first_year = sim1.date[0].year
        last_year = sim1.date[len(sim1)-1].year
        sim1.set_index('date', inplace=True)
        # Summer flow
        for year in range(first_year, last_year+1):
            start_date = str(year)+'-06-01'
            end_date = str(year)+'-08-31'
            sim2 = sim1.loc[start_date : end_date]
            summer_sum += sim2['q_tot_x=0.0'].sum()
        # Winter flow
        for year in range(first_year, last_year+1):
            start_date = str(year)+'-01-01'
            end_date = str(year)+'-03-01'
            sim2 = sim1.loc[start_date : end_date]
            winter_sum += sim2['q_tot_x=0.0'].sum()
            start_date = str(year)+'-12-01'
            end_date = str(year)+'-12-01'
            sim2 = sim1.loc[start_date : end_date]
            winter_sum += sim2['q_tot_x=0.0'].sum()
        sim1.reset_index(inplace=True)
        simsum += sim1['q_tot_x=0.0'].sum()
        ndata += len(sim1)
    f_pest.write(" %12.6f\n"%simsum)
    f_pest.write(" %12.6f\n"%summer_sum)
    f_pest.write(" %12.6f\n"%winter_sum)
    f_pest.write(" %12.6f"%(simsum/ndata))
    f_pest.close()
    
def write_userdef_file(Sim, fo):
    
#    import os
    from importlib.machinery import SourceFileLoader
#    cwd = os.getcwd()
    name2="write_userdef_file.py"
#    name2=os.path.join(cwd,name)
#    print(cwd)
#    print(name)
#    print(name2)
    try:
        wpf = SourceFileLoader("wpf",name2).load_module()
        wpf.write_userdef_file(Sim)
    except:
        print_msg("\nCould not find or run script named %s"%name2, fo)
        # print("\nCould not find or run script named %s"%name2)
        # fo.write("\nCould not find or run script named %s"%name2)

def print_msg(msg, fl):
    """ prints msg on screen and in log file """
    print(msg)
    fl.write("\n"+msg)        
    
###############################################################################    

def run_model(yaml='flowsim.yaml', log='flowsim.log'):

#    dir_path = os.path.dirname(os.path.realpath(__file__))
#    cj_module = os.path.join(dir_path, "cj_response_funcs.py")
#    cj = SourceFileLoader("cj", cj_module).load_module()

    from flowsim import version
    __version__ = version.__version__

    m = model()

    try:
        m.fo = open(log,'w')
        print_msg("\nRunning flowsim version %s\n"%__version__, m.fo)
    except:
        print("Error: Cannot open logfile: %s"%log)
        return
    
    ### Load input, initialize, and check that all required input is available
    if m.load_input_file(yaml):      # Load yaml file
        if m.initialize_read():  # Check input and set values
            Sim = m.simulation() # Make simulation results
            f_res = open("flowsim-res.csv","w")
    #        f_pest = open("flowsim-PEST.csv","w")
    #        simsum = 0.0
            for sim in Sim:
                colnams = list()
                for nam in sim.columns:
                    colnams.append(nam)
                colnams[0] = colnams[0].rjust(10)
                for i in range(1,len(colnams)):
                    colnams[i] = colnams[i].rjust(15)
                sim.to_csv(f_res, index=False, header=colnams, 
                           float_format=" %14.6f", line_terminator="\n")
    #            sim.to_csv(f_pest, index=False, columns=['q_Shallow_x=0.0'], 
    #                       float_format=" %12.6f", line_terminator="\n")
    #            simsum += sim['q_tot_x=0.0'].sum()
    #        f_pest.write("%12.2f"%simsum)
            f_res.close()
    #        f_pest.close()
            if m.userdef_output:
                write_userdef_file(Sim, m.fo)
                # import write_pest_file as wpf
                # wpf.write_pest_read_file(Sim)
            
            if m.makeplot:
                for key in m.plot:
                    plot_response(m, key, Sim, yax = 'lin', 
                                  resamp_rule = m.steplength)
                
                plt.close()
    m.fo.close()
