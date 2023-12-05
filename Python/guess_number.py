guess = None
answer = 77
cnt = 0
limit = 3

while guess != answer and cnt < limit:
    cnt += 1
    guess = int(input("enter your guess ==> "))
    if guess > answer:
        print("smaller")
    elif guess < answer:
        print("bigger")

if cnt <= limit and guess == answer:
    print("you win")
else:
    print("you lose")