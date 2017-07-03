import pickle
sp_dict = {'one':'uno', 'two':'dos', 'three':'tres'}
with open('spdict.data', 'wb') as fh:
    pickle.dump(sp_dict, fh)
