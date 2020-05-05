""" Run the functions"""
# --------  Run Model and Plot it -------- #
[Totalendlist, VIRendlist, gamma_list2] =  equation_run(temp_list, mux, betx, phix, gamx, time, runs)

temp_list2=list(itertools.chain.from_iterable(itertools.repeat(x, runs) for x in temp_list))
size_list = [i*500 for i in gamma_list2]
lytic = [i*100 for i in gamma_list2]
source = ColumnDataSource(data=dict(x=Totalendlist, y=VIRendlist, temp_list2 = temp_list2, size_list = size_list, lytic = lytic))

hover = HoverTool(tooltips=[("Ice Temp (C)","@temp_list2"),("Bacterial Population", "$x"),("Viral Population", "$y"), ("Lytic %", "@lytic")])
mapper = LinearColorMapper(palette=viridis(256), low=min(temp_list2), high=max(temp_list2))

p = figure(plot_height=500, plot_width=500, x_axis_type = 'log', y_axis_type = 'log',x_range = (1e4, 1e8), y_range = (1e5,1e10), tools = [hover])
p.xaxis.axis_label = 'Bacteria per mL'
p.yaxis.axis_label = 'Virus per mL'


# plot the ratio lines 
Ratiolines = np.linspace(1e4, 1e8,1000)
Lines = p.line(Ratiolines, 1*Ratiolines)
Lines2 = p.line(Ratiolines, 10*Ratiolines)
Lines3 = p.line(Ratiolines, 100*Ratiolines)
Lines4 =p.line(Ratiolines, 1000*Ratiolines)
Lines5= p.line(Ratiolines, 10000*Ratiolines)
toggle2 = Toggle(label="Ratio Lines", button_type="success", active=True)
toggle2.js_link('active', Lines , 'visible')
toggle2.js_link('active', Lines2 , 'visible')
toggle2.js_link('active', Lines3 , 'visible')
toggle2.js_link('active', Lines4 , 'visible')
toggle2.js_link('active', Lines5 , 'visible')

#plot the data
p.circle('x', 'y', fill_color= transform('temp_list2', mapper), size = 'size_list', fill_alpha=0.6, line_color=None, source = source)
layout = row(
    p, toggle2
)


show(layout)
