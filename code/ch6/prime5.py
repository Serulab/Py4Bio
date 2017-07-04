def is_prime(n):
   """ Check if n is a prime number.
   Sample usage:
   >>> is_prime(0)
   False
   >>> is_prime(1)
   True
   >>> is_prime(2)
   True
   >>> is_prime(3)
   True
   >>> is_prime(4)
   False
   >>> is_prime(5)
   True
   """

   if n <= 0:
       # This is only for numbers > 0.
       return False
   for x in range(2, n):
       if n%x == 0:
           return False
   return True

def _test():
   import doctest
   doctest.testmod()

if __name__ == '__main__':
   _test()
