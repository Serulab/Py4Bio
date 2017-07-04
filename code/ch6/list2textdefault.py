def save_list(input_list, file_name='temp.txt'):
    """A list (input_list) is saved in a file (file_name)"""
    with open(file_name, 'w') as fh:
        for item in input_list:
            fh.write('{0}\n'.format(item))
    return None
