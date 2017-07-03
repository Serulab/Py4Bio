color_code = [('red', 1), ('green', 2), ('blue', 3), ('black', 4)]
name = 'blue'
for color_pair in color_code:
    if name == color_pair[0]:
        code = color_pair[1]
print(code)
