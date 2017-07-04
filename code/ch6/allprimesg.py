def g_all_primes(n):
    for number in range(1,n):
        if is_prime(number):
            yield number
