color_code = [('red',1), ('green',2), ('blue',3), ('black',4)]
name = 'blue'
i = 0
while name != color_code[i][0]:
    i += 1
code = color_code[i][1]
print(code)
