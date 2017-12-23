#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDpred, Special script:

Here we produce an interactive web based plot to analyze various set of
parameters to calculate SASA.

#Script requires bokeh library and dependencies
#Install them with "pip install -r H-SASA_params_req.txt"

"""
from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, show, output_notebook
from bokeh.plotting import figure
from bokeh.charts import Bar
from bokeh.io import output_file, show
from bokeh.layouts import widgetbox, row, column, layout
from bokeh.models.widgets import CheckboxGroup, Slider, RadioGroup, RadioButtonGroup, Div
from bokeh.models import CustomJS, ColumnDataSource, Range1d
from bokeh.models.annotations import Label
from bokeh.models import SingleIntervalTicker, LinearAxis

import pandas as pd
import numpy as np

dfNdist=pd.read_csv("results/nucl_H_Ndist_multparam_H-SASA.csv.gz")
dfELdist=pd.read_csv("results/nucl_H_ELdist_multparam_H-SASA.csv.gz")
dfMD1=pd.read_csv("results/nucl_H_MD1ns_multparam_H-SASA.csv.gz")
dfMD2=pd.read_csv("results/nucl_H_MD2ns_multparam_H-SASA.csv.gz")
dfMDav=pd.read_csv("../exampleMD/results/nucl_MD_multparam_H-SASA.csv.gz")
# df1najm1=pd.read_csv("results/1NAJm1_multparam_H-SASA.csv")
# dfscmin=pd.read_csv("results/SCmin_multiparam_H-SASA.csv")
# dfscall=pd.read_csv("results/SCall_multiparam_H-SASA.csv")
dfNdist['struct']='X-H nuclear'
dfELdist['struct']='X-H electron'
dfMD1['struct']='MD 1ns'
dfMD2['struct']='MD 2ns'
dfMDav['struct']='MD averaged'
# df1najm1['struct']='1NAJ model 1'
# dfscmin['struct']='NMR+LAXS min'
# dfscall['struct']='NMR+LAXS 50'
df=pd.concat([dfNdist,dfELdist,dfMD1,dfMD2,dfMDav])
df=df[(df['resid_0']<=60)&(df['resid_0']>=-60)]

output_file("results/nucl_H-SASA.html")
# output_notebook()

cbox = CheckboxGroup(
        labels=["H1\'", "H2\' H2\'\'", "H3\'","H4\'","H5\' H5\'\'"], active=[0,1,2,3,4],name='Contributing atoms')
slider = Slider(start=1.0, end=5.0, value=1.4, step=.1, title="Probe radius, A",orientation='horizontal')
rb_group = RadioGroup(
        labels=["FreeSASA Default", "CHARMM36-rmin", "AMBER10-rmin"], active=1,name='Atom radii')
rb_group2 = RadioGroup(
        labels=["X-H nuclear","X-H electron","MD 1ns","MD 2ns",'MD averaged'], active=4,name='Structure')

x=list(df[(df['struct']=='MD averaged')&(df['vdw_set']=='charmm36-rmin')&(df['probe_size']==1.4)]['resid_0'].values)

probe_size=ColumnDataSource(data=dict(probe_size=[1.4]))
vdw_set=ColumnDataSource(data=dict(vdw_set=['charmm36-rmin']))
struct=ColumnDataSource(data=dict(struct=['MD averaged']))

yt=df[(df['struct']==struct.data['struct'][0])&(df['vdw_set']==vdw_set.data['vdw_set'][0])&(df['probe_size']==probe_size.data['probe_size'][0])][['H-SASA_0','H-SASA_1','H-SASA_2','H-SASA_3','H-SASA_4']].sum(axis=1).values
yb=df[(df['struct']==struct.data['struct'][0])&(df['vdw_set']==vdw_set.data['vdw_set'][0])&(df['probe_size']==probe_size.data['probe_size'][0])][['H-SASA_5','H-SASA_6','H-SASA_7','H-SASA_8','H-SASA_9']].sum(axis=1).values

yh=df[(df['struct']==struct.data['struct'][0])&(df['vdw_set']==vdw_set.data['vdw_set'][0])&(df['probe_size']==probe_size.data['probe_size'][0])][['H-SASA_0','H-SASA_1','H-SASA_2','H-SASA_3','H-SASA_4','H-SASA_5','H-SASA_6','H-SASA_7','H-SASA_8','H-SASA_9']]

ym=df[['struct','probe_size','vdw_set','H-SASA_0','H-SASA_1','H-SASA_2','H-SASA_3','H-SASA_4','H-SASA_5','H-SASA_6','H-SASA_7','H-SASA_8','H-SASA_9']]


source = ColumnDataSource(data=dict(x=x, yt=yt,yb=yb))
origf = ColumnDataSource(data=yh)
orig = ColumnDataSource(data=ym)

p1 = figure(plot_width=1500, plot_height=400,x_range=Range1d(-60,60),y_range=Range1d(0, 90),toolbar_location=None,\
    x_axis_label='DNA sequence',y_axis_label='H-SASA, A^2',title=None,\
    x_axis_type=None)
p1.line('x', 'yt', source=source, line_color="green", line_width=3, legend="H-SASA Top strand")
p1.circle('x', 'yt', source=source,size=5, fill_color="green", line_color="green", line_width=3, )
p1.legend.location = "top_right"
p1.line('x', 'yb', source=source, line_color="blue", line_width=3, legend="H-SASA Bottom strand")
p1.circle('x', 'yb', source=source,size=5, fill_color="blue", line_color="blue", line_width=3, )

ticker = SingleIntervalTicker(interval=5, num_minor_ticks=10)
xaxis = LinearAxis(ticker=ticker,axis_label='DNA sequence')
p1.add_layout(xaxis, 'below')

cb = CustomJS(args=dict(source=source,origf=origf,cbox=cbox), code="""
        var data = source.data;
        var origfd = origf.data;
        var val = cbox.active;
        
        yt = data['yt']
        yb = data['yb']
        for (i = 0; i < yt.length; i++) {
        yt[i]=0
        yb[i]=0
        for (var k = 0, len = val.length; k < len; k++) {
        yt[i] = yt[i]+origfd['H-SASA_'.concat(val[k].toString())][i]
        yb[i] = yb[i]+origfd['H-SASA_'.concat((val[k]+5).toString())][i]
        }
        }

        source.trigger('change');
        
    """)

cb_slider = CustomJS(args=dict(source=source,slider=slider,origf=origf,orig=orig,struct=struct,probe_size=probe_size,vdw_set=vdw_set,cbox=cbox), code="""
        var val = cbox.active;
        var data = source.data;
        var origfd = origf.data;
        var origd = orig.data;
        var valp = slider.value;
        probe_size.data['probe_size'][0]=valp;
        probe_size.trigger('change');
        var probe_sized = valp;
        var vdw_set = vdw_set.data['vdw_set'][0];
        var struct = struct.data['struct'][0];
        
        
        var k=0;
        for (i = 0; i < origd['H-SASA_0'].length; i++){
        if(origd['probe_size'][i]===probe_sized){
        if(origd['vdw_set'][i]===vdw_set){
        if(origd['struct'][i]===struct){
        origfd['H-SASA_0'][k]=origd['H-SASA_0'][i]
        origfd['H-SASA_1'][k]=origd['H-SASA_1'][i]
        origfd['H-SASA_2'][k]=origd['H-SASA_2'][i]
        origfd['H-SASA_3'][k]=origd['H-SASA_3'][i]
        origfd['H-SASA_4'][k]=origd['H-SASA_4'][i]
        origfd['H-SASA_5'][k]=origd['H-SASA_5'][i]
        origfd['H-SASA_6'][k]=origd['H-SASA_6'][i]
        origfd['H-SASA_7'][k]=origd['H-SASA_7'][i]
        origfd['H-SASA_8'][k]=origd['H-SASA_8'][i]
        origfd['H-SASA_9'][k]=origd['H-SASA_9'][i]
        k=k+1;
        }
        }
        }
        }

        origf.trigger('change');
        
        yt = data['yt']
        yb = data['yb']
        for (i = 0; i < yt.length; i++) {
        yt[i]=0
        yb[i]=0
        for (var k = 0, len = val.length; k < len; k++) {
        yt[i] = yt[i]+origfd['H-SASA_'.concat(val[k].toString())][i]
        yb[i] = yb[i]+origfd['H-SASA_'.concat((val[k]+5).toString())][i]
        }
        }

        source.trigger('change');
    """)

cb_ff = CustomJS(args=dict(source=source,rb_group=rb_group,slider=slider,origf=origf,orig=orig,struct=struct,probe_size=probe_size,vdw_set=vdw_set,cbox=cbox), code="""
        var val = cbox.active;
        var data = source.data;
        var origfd = origf.data;
        var origd = orig.data;
        var valp = slider.value;
        var valff = rb_group.active;
        var struct = struct.data['struct'][0];
        var probe_sized = valp;
        if(valff==0){vdw_set.data['vdw_set'][0]='None'; vdw_setd='None';}
        if(valff==1){vdw_set.data['vdw_set'][0]='charmm36-rmin'; vdw_setd='charmm36-rmin';}
        if(valff==2){vdw_set.data['vdw_set'][0]='amber10-rmin'; vdw_setd='amber10-rmin';}
        vdw_set.trigger('change');


        var k=0;
        for (i = 0; i < origd['H-SASA_0'].length; i++){
        if(origd['probe_size'][i]===probe_sized){
        if(origd['vdw_set'][i]===vdw_setd){
        if(origd['struct'][i]===struct){
        origfd['H-SASA_0'][k]=origd['H-SASA_0'][i]
        origfd['H-SASA_1'][k]=origd['H-SASA_1'][i]
        origfd['H-SASA_2'][k]=origd['H-SASA_2'][i]
        origfd['H-SASA_3'][k]=origd['H-SASA_3'][i]
        origfd['H-SASA_4'][k]=origd['H-SASA_4'][i]
        origfd['H-SASA_5'][k]=origd['H-SASA_5'][i]
        origfd['H-SASA_6'][k]=origd['H-SASA_6'][i]
        origfd['H-SASA_7'][k]=origd['H-SASA_7'][i]
        origfd['H-SASA_8'][k]=origd['H-SASA_8'][i]
        origfd['H-SASA_9'][k]=origd['H-SASA_9'][i]
        k=k+1;
        }
        }
        }
        }

        origf.trigger('change');
        
        yt = data['yt']
        yb = data['yb']
        for (i = 0; i < yt.length; i++) {
        yt[i]=0
        yb[i]=0
        for (var k = 0, len = val.length; k < len; k++) {
        yt[i] = yt[i]+origfd['H-SASA_'.concat(val[k].toString())][i]
        yb[i] = yb[i]+origfd['H-SASA_'.concat((val[k]+5).toString())][i]
        }
        }

        source.trigger('change');
    """)

cb_struct = CustomJS(args=dict(source=source,rb_group2=rb_group2,slider=slider,origf=origf,orig=orig,struct=struct,probe_size=probe_size,vdw_set=vdw_set,cbox=cbox), code="""
        var val = cbox.active;
        var origd = orig.data;
        var data = source.data;
        var origfd = origf.data;
        var valp = slider.value;
        var valstruct = rb_group2.active;
        var vdw_setd = vdw_set.data['vdw_set'][0];
        var probe_sized = valp;
        if(valstruct==0){struct.data['struct'][0]='X-H nuclear'; structd='X-H nuclear';}
        if(valstruct==1){struct.data['struct'][0]='X-H electron'; structd='X-H electron';}
        if(valstruct==2){struct.data['struct'][0]='MD 1ns'; structd='MD 1ns';}
        if(valstruct==3){struct.data['struct'][0]='MD 2ns'; structd='MD 2ns';}
        if(valstruct==4){struct.data['struct'][0]='MD averaged'; structd='MD averaged';}
        struct.trigger('change');

        
        var k=0;
        for (i = 0; i < origd['H-SASA_0'].length; i++){
        if(origd['probe_size'][i]===probe_sized){
        if(origd['vdw_set'][i]===vdw_setd){
        if(origd['struct'][i]===structd){
        origfd['H-SASA_0'][k]=origd['H-SASA_0'][i]
        origfd['H-SASA_1'][k]=origd['H-SASA_1'][i]
        origfd['H-SASA_2'][k]=origd['H-SASA_2'][i]
        origfd['H-SASA_3'][k]=origd['H-SASA_3'][i]
        origfd['H-SASA_4'][k]=origd['H-SASA_4'][i]
        origfd['H-SASA_5'][k]=origd['H-SASA_5'][i]
        origfd['H-SASA_6'][k]=origd['H-SASA_6'][i]
        origfd['H-SASA_7'][k]=origd['H-SASA_7'][i]
        origfd['H-SASA_8'][k]=origd['H-SASA_8'][i]
        origfd['H-SASA_9'][k]=origd['H-SASA_9'][i]
        k=k+1;
        }
        }
        }
        }

        origf.trigger('change');
        
        yt = data['yt']
        yb = data['yb']
        for (i = 0; i < yt.length; i++) {
        yt[i]=0
        yb[i]=0
        for (var k = 0, len = val.length; k < len; k++) {
        yt[i] = yt[i]+origfd['H-SASA_'.concat(val[k].toString())][i]
        yb[i] = yb[i]+origfd['H-SASA_'.concat((val[k]+5).toString())][i]
        }
        }

        source.trigger('change');
    """)


cbox.js_on_change('active', cb)
slider.js_on_change('value', cb_slider)
rb_group.js_on_change('active', cb_ff)
rb_group2.js_on_change('active', cb_struct)

div = Div(text="""<b>Figure: H-SASA profiles calculated with different parameters for a
    3D structural model of <i>S. cerevisiae</i> centromeric nucleosome with 601TA DNA sequence </b>
<br><br>
This interactive plot allows to explore dependence of H-SASA profiles on different paramters.
Use contols above to choose between following paramters:<br><br>
<b>Structure used for calculations</b>: X-H nuclear - orginal X-ray derived structure, hydrogen atoms added via REDUCE program with nuclear distances for X-H bond length;
X-H electron - original X-ray derived structure, hydrogen atoms added via REDUCE program with electron cloud distances for X-H bond length;
MD 1 ns - X-ray derived structure after 1 ns of molecular dynamics simulations with CHARMM36 force field with NAMD;
MD 2 ns - X-ray derived structure after 2 ns of molecular dynamics simulations with CHARMM36 force field with NAMD;
MD average - average profile for 50 frames spaced 1 ns apart from molecular dynamics simulations with CHARMM36 force field with NAMD.
<br><br>
<b>Radii of atoms used for SASA calculations</b>: FreeSASA Default - default radii in FreeSASA program; CHARMM36-rmin - van der Waals radii (Rmin paramter) of atoms as defined by CHARMM36 force field;
AMBER10-rmin - van der Waals radii (Rmin paramter) of atoms as defined by AMBER10 force field.
<br><br>
<b>Contributions of deoxyribose hydrogen atoms</b>: SASA of deoxyribose hydrogen atoms selected will be included in calculation of H-SASA profile.
<br><br>
<b>Probe radius</b>: radius of the probe sphere used to calculate SASA, 1.4 A is a default for SASA calculations, OH-radical size is around 1.2 A, higher values are more sensitive to the geometry of the complex.

    """,
width=1000, height=100)



# l1=layout([[p1,column(widgetbox(rb_group2),widgetbox(rb_group),widgetbox(cbox),widgetbox(slider))]])
l1=layout([[p1],[widgetbox(rb_group2),widgetbox(rb_group),widgetbox(cbox),widgetbox(slider)],[widgetbox(div)]])
show(l1)