class Zdic(dict):
    """ A dictionary-like object that return 0 when a user
        request a non-existent key.
    """
    
    def __missing__(self,x):
        return 0
