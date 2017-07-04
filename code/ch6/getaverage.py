def average(*numbers):
    if len(numbers)==0:
        return None
    else:
        total = sum(numbers)
        return total / len(numbers)
