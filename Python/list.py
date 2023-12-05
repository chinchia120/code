from ast import Num
from cmath import phase
from pickle import TRUE
from statistics import mode
from tabnanny import check
from telnetlib import PRAGMA_HEARTBEAT, SGA
from math import *

phase = "Hello"
print(phase.index("o"))
print(phase.replace("H", "h"))
print(round(3.66))
print(floor(4.6))
print(ceil(5.1))

num1 = input("enter num1 ==> ")
num2 = input("enter num2 ==> ")
print(float(num1) + float(num2))

score = [100, "AAA", 80, 70, True]
friend = ["BBB", "CCC", "DDD"]
print(score)
print(score[: -1])

score.extend(friend)
print(score)

score.append(30)
print(score)

score.insert(2, 30)
print(score)

score.remove(100)
print(score)

score.pop()
print(score)

score = [100, 90, 80, 70, 10, 30, 90]
score.sort()
print(score)

score.reverse()
print(score)

print(score.index(90))

print(score.count(90))

score.clear()
print(score)