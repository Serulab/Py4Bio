def save_list(input_list, file_name='temp.txt'):
    """A list (input_list) is saved to a file (file_name)"""
    with open(file_name, 'w') as fh:
        print(*input_list, sep='\n', file=fh)
    return None
