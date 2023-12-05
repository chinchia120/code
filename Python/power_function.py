def power(base_num, pow_num):
    sum = 1
    for i in range(pow_num):
        sum *= base_num
        print(i)

    return sum


print(power(2, 6))