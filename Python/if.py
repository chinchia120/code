hungry = True

if hungry:
    print("eat")
elif hungry == False:
    print("don't eat")

score = 40

if score == 100:
    print("Nice")
elif score >= 80:
    print("Good")
elif score >= 60:
    print("Not Bad")
else:
    print("Bad")


def max_num(num1, num2, num3):
    if(num1 >= num2 and num1 >= num3):
        return num1
    elif(num2 >= num1 and num2 >= num3):
        return num2
    else:
        return num3

print(max_num(10, 20, 15))