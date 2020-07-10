# Pulsar-Catalogue-Summer-2020
# Contains the updated Fermi Catalogue for Summer 2020, including the addition of the Fermi Detected but Not Flagged Pulsars highlighted in red.


import numpy as np
from scipy.special import gammaln
import matplotlib.pyplot as plt , plt.rcParams["figure.figsize"]=[40,40]
from scipy import misc
from psrqpy import QueryATNF
query=QueryATNF()
import os 
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import QTable, Table, Column, MaskedColumn, join
t = Table.read('https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/gll_psc_v23.fit') #Gets FERMI source catalog
#Clean the Fermi Table and add a corresponding PSRJ column
t = Table(t, masked=True, copy=False) 
t['PSRJ'] = ''
t['PSRJ'] = MaskedColumn(t['PSRJ'],dtype='S12')
t['PSRJ'].mask = True
for ind,source in enumerate(t):
    class1 = source['CLASS1']
    class2 = source['CLASS2']
    t['CLASS1'][ind] = class1.strip()
    t["CLASS2"][ind] = class2.strip()
    if t['CLASS1'][ind] == 'PSR' or t['CLASS1'][ind] == 'psr':
        name = source['ASSOC1']
        parts = name.split(' ')
        t['PSRJ'][ind] = parts[1]
        t['PSRJ'].mask[ind] = False
    if t["CLASS2"][ind] == 'psr': #Note that there are no CLASS2 = PSR sources
        name = source['ASSOC2']
        parts = name.split(' ')
        t['PSRJ'][ind] = parts[1]
        t['PSRJ'].mask[ind] = False
print(set(t['CLASS1']))
print(set(t['CLASS2']))
print(t[t['CLASS2'] == 'glc'])
print(sum(t['PSRJ'].mask == True))
print(sum(t['PSRJ'].mask == False))
combined_table = join(query.table, t.filled(), keys='PSRJ',join_type='outer') #joins the tables
combined_table = Table(combined_table, masked=True, copy=False)
colors=["grey","lime","magenta","cyan"]
pulsar_periodderivatives=query.table['P1']
Not_Fermi=[np.logical_and([query.table['SURVEY']!='FermiAssoc'],np.logical_and([query.table['SURVEY']!='FermiBlind'],[query.table['SURVEY']!='FermiAssoc,ghrss']))]
pulsar_periods=query.table['P0']
Not_Fermi=np.array(Not_Fermi)
Not_Fermi=Not_Fermi.flatten()
def rotate(lst, x):    
    return lst[-x:] + lst[:-x]
mpl.rcParams.update({'font.size': 36})
fig = plt.figure(figsize=(40, 40))  
scatter_axes = plt.subplot2grid((3, 3), (1, 0), rowspan=2, colspan=2, fig=fig)
x_hist_axes = plt.subplot2grid((3, 3), (0, 0), colspan=2,
                               sharex=scatter_axes,fig=fig )
y_hist_axes = plt.subplot2grid((3, 3), (1, 2), rowspan=2,
                               sharey=scatter_axes, fig=fig)



scatter_axes.plot(query.table['P0'],query.table['P1'] , '.',color='grey',label='Non-FermiPulsars',markersize=20)
scatter_axes.plot(query.table['P0'][query.table['SURVEY']=='FermiBlind'],
            query.table['P1'][query.table['SURVEY']=='FermiBlind'],'.',color=colors[3],label='FermiBlind',markersize=20)
scatter_axes.set_xscale('log')
scatter_axes.set_yscale('log')
scatter_axes.set_ylim(1e-22, 1e-9)
scatter_axes.set_xlim(0.001,50)
scatter_axes.set_xlabel("period (s)")
scatter_axes.set_ylabel("period derivative (s/s)")



scatter_axes.plot(query.table['P0'][query.table['SURVEY']=='FermiAssoc'],
            query.table['P1'][query.table['SURVEY']=='FermiAssoc'],'.',color=colors[2],label='FermiAssoc',markersize=20)
#x_hist_axes.hist(positive['P0'],bins=np.logspace(np.log10(1e-3),np.log10(3.0),30))
#y_hist_axes.hist(positive['P1'],bins=np.logspace(np.log10(1e-22),np.log10(1e-9),50),orientation='horizontal')
#x_hist_axes.set_xscale('log')
#y_hist_axes.set_xscale('log')
#y_hist_axes.set_xlim(1e-22, 1e-9)
#x_hist_axes.set_xlim(0.001,50)
#x_hist_axes.set_xlabel("period (s)")
#y_hist_axes.set_xlabel("period derivative (s/s)")

scatter_axes.plot(query.table['P0'][query.table['SURVEY']=='FermiAssoc,ghrss'],
            query.table['P1'][query.table['SURVEY']=='FermiAssoc,ghrss'],'.',color=colors[1],label='FermiAssoc,ghrss',markersize=20)
scatter_axes.legend()


x_hist_axes.hist([pulsar_periods[query.table['SURVEY']=='FermiAssoc,ghrss'],
         pulsar_periods[query.table['SURVEY']=='FermiAssoc'],
          pulsar_periods[query.table['SURVEY']=='FermiBlind'],pulsar_periods[Not_Fermi]],
          bins=np.logspace(np.log10(1e-3),np.log10(50),30),color=rotate(colors[0:4],-1),stacked=True)
y_hist_axes.hist([pulsar_periodderivatives[query.table['SURVEY']=='FermiAssoc,ghrss'],
         pulsar_periodderivatives[query.table['SURVEY']=='FermiAssoc'],
          pulsar_periodderivatives[query.table['SURVEY']=='FermiBlind'],pulsar_periodderivatives[Not_Fermi]],
          bins=np.logspace(np.log10(1e-22),np.log10(1e-9),50),orientation='horizontal',color=rotate(colors[0:4],-1),stacked=True)
x_hist_axes.set_yscale("log")
y_hist_axes.set_xscale("log")
x_hist_axes.set_ylim(0.5,500)
y_hist_axes.set_xlim(0.5,300)
from astropy import units as u
from astropy.coordinates import SkyCoord
combined_table['FERMIcoord'] = SkyCoord(combined_table['RAJ2000'], combined_table['DEJ2000'], unit=(u.deg, u.deg))
combined_table['ATNFcoord'] = SkyCoord(combined_table['RAJ'], combined_table['DECJ'], unit=(u.hourangle, u.deg))
a = np.array([combined_table['PSRJ'] == 'N'])
a = a.flatten()
print(sum(a))
combined_table.mask['PSRJ'] = Column(a)
print(sum(combined_table.mask['PSRJ']))
combined_table.mask['FERMIcoord'] = combined_table.mask['RAJ2000']
combined_table.mask['ATNFcoord'] = combined_table.mask['RAJ']
print(sum(combined_table.mask['FERMIcoord']))
print(sum(combined_table.mask['ATNFcoord']))

colors=["grey","lime","magenta","cyan","red"]
pulsar_periodderivatives=combined_table['P1']
Not_Fermi_ATNF=[np.logical_and([combined_table['SURVEY']!='FermiAssoc'],np.logical_and([combined_table['SURVEY']!='FermiBlind'],[combined_table['SURVEY']!='FermiAssoc,ghrss']))]
Not_Fermi_ATNF=np.array(Not_Fermi_ATNF)
Not_Fermi_ATNF=Not_Fermi_ATNF.flatten()
Fermi_Detected=[np.logical_and(combined_table.mask['FERMIcoord']==False,combined_table.mask['ATNFcoord']==False)]
Fermi_Detected=np.array(Fermi_Detected)
Fermi_Detected=Fermi_Detected.flatten()
Fermi_Detected_Not_Flagged=[np.logical_and(Not_Fermi_ATNF,Fermi_Detected)]
Fermi_Detected_Not_Flagged=np.array(Fermi_Detected_Not_Flagged)
Fermi_Detected_Not_Flagged=Fermi_Detected_Not_Flagged.flatten()
Not_Fermi=[np.logical_and(combined_table.mask['FERMIcoord']==True,combined_table.mask['ATNFcoord']==False)]
pulsar_periods=combined_table['P0']
Not_Fermi=np.array(Not_Fermi)
Not_Fermi=Not_Fermi.flatten()
def rotate(lst, x):    
    return lst[-x:] + lst[:-x]
mpl.rcParams.update({'font.size': 36})
fig = plt.figure(figsize=(40, 40))  
scatter_axes = plt.subplot2grid((3, 3), (1, 0), rowspan=2, colspan=2, fig=fig)
x_hist_axes = plt.subplot2grid((3, 3), (0, 0), colspan=2,
                               sharex=scatter_axes,fig=fig )
y_hist_axes = plt.subplot2grid((3, 3), (1, 2), rowspan=2,
                               sharey=scatter_axes, fig=fig)

scatter_axes.plot(combined_table['P0'],combined_table['P1'] , '.',color='grey',label='Non-FermiPulsars',markersize=20)
scatter_axes.plot(combined_table['P0'][combined_table.mask['FERMIcoord']==False],combined_table['P1'][combined_table.mask['FERMIcoord']==False],'.',color='red',label='Detected Pulsars',markersize=20)
scatter_axes.plot(combined_table['P0'][combined_table['SURVEY']=='FermiBlind'],
            combined_table['P1'][combined_table['SURVEY']=='FermiBlind'],'.',color=colors[3],label='FermiBlind',markersize=20)
scatter_axes.set_xscale('log')
scatter_axes.set_yscale('log')
scatter_axes.set_ylim(1e-22, 1e-9)
scatter_axes.set_xlim(0.001,50)
scatter_axes.set_xlabel("period (s)")
scatter_axes.set_ylabel("period derivative (s/s)")



scatter_axes.plot(combined_table['P0'][combined_table['SURVEY']=='FermiAssoc'],
            combined_table['P1'][combined_table['SURVEY']=='FermiAssoc'],'.',color=colors[2],label='FermiAssoc',markersize=20)
#x_hist_axes.hist(positive['P0'],bins=np.logspace(np.log10(1e-3),np.log10(3.0),30))
#y_hist_axes.hist(positive['P1'],bins=np.logspace(np.log10(1e-22),np.log10(1e-9),50),orientation='horizontal')
#x_hist_axes.set_xscale('log')
#y_hist_axes.set_xscale('log')
#y_hist_axes.set_xlim(1e-22, 1e-9)
#x_hist_axes.set_xlim(0.001,50)
#x_hist_axes.set_xlabel("period (s)")
#y_hist_axes.set_xlabel("period derivative (s/s)")

scatter_axes.plot(combined_table['P0'][combined_table['SURVEY']=='FermiAssoc,ghrss'],
            combined_table['P1'][combined_table['SURVEY']=='FermiAssoc,ghrss'],'.',color=colors[1],label='FermiAssoc,ghrss',markersize=20)
scatter_axes.legend()


x_hist_axes.hist([pulsar_periods[combined_table['SURVEY']=='FermiAssoc,ghrss'],
         pulsar_periods[combined_table['SURVEY']=='FermiAssoc'],
          pulsar_periods[combined_table['SURVEY']=='FermiBlind'],pulsar_periods[Fermi_Detected_Not_Flagged],pulsar_periods[Not_Fermi]],
          bins=np.logspace(np.log10(1e-3),np.log10(50),30),color=rotate(colors[0:5],-1),stacked=True)
y_hist_axes.hist([pulsar_periodderivatives[combined_table['SURVEY']=='FermiAssoc,ghrss'],
         pulsar_periodderivatives[combined_table['SURVEY']=='FermiAssoc'],
          pulsar_periodderivatives[combined_table['SURVEY']=='FermiBlind'],pulsar_periodderivatives[Fermi_Detected_Not_Flagged],pulsar_periodderivatives[Not_Fermi]],
          bins=np.logspace(np.log10(1e-22),np.log10(1e-9),50),orientation='horizontal',color=rotate(colors[0:5],-1),stacked=True)
x_hist_axes.set_yscale("log")
y_hist_axes.set_xscale("log")
x_hist_axes.set_ylim(0.5,500)
y_hist_axes.set_xlim(0.5,300)
combined_table[Fermi_Detected_Not_Flagged]











