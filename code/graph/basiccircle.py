from bokeh.plotting import figure, output_file, show

p = figure(width=400, height=400)
p.circle(2, 3, radius=.5, alpha=0.5)
output_file("out.html")
show(p)
