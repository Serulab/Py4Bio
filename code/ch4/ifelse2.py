three_letter_code = {'A':'Ala', 'N':'Asn', 'D':'Asp', 'C':'Cys'}
aa = input('Enter one letter: ')
if aa in three_letter_code:
    print('The three letter code for {0} is {1}'.format(aa,
          three_letter_code[aa]))
else:
    print("Sorry, I don't have it in my dictionary")
