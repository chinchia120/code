num1 = float(input("enter num1 ==> ")) 
ope = input("enter operator ==> ")
num2 = float(input("enter num2 ==> ")) 

if ope == "+":
    print(num1 + num2)
elif ope == "-":
    print(num1 - num2)
elif ope == "*":
    print(num1 * num2)
elif ope == "/":
    print(num1 / num2)
else:
    print("error\n") 