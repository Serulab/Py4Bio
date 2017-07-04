def is_prime(n):
   """Returns True is n is prime, False if not"""
   for i in range(2,n-1):
       if n%i == 0:
           return False
   return True

def all_primes(n):
   primes = []
   for number in range(1,n):
       if isprime(number):
           primes.append(number)
   return p
