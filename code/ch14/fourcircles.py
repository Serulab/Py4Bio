from bokeh.plotting import figure, output_file, show

p = figure(width=500, height=500)
x = [1, 1, 2, 2]
y = [1, 2, 1, 2]
p.circle(x, y, radius=.35, alpha=0.5, color='red')
output_file("out.html")
show(p)
