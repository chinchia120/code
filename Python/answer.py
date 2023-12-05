#!/usr/bin/python
# -*- coding: utf-8 -*-
from question import Question

test = [
    "1 + 3 = ?\n(a) 2\n(b) 3\n(c) 4\n\n",
    "1公尺等於幾公分 ?\n(a) 10\n(b) 100\n(c) 1000\n\n",
    "香蕉是什麼顏色 ?\n(a) 黑色\n(b) 黃色\n(c) 白色\n\n"
]

questions = [
    Question(test[0], "c"),
    Question(test[1], "b"),
    Question(test[2], "b")
]

def run_test(questions):
    score = 0
    for question in questions:
        answer = input(question.description)
        if answer == question.answer:
            score += 1
            
    print("score = " + str(score) + "/" + str(len(questions)))
    
run_test(questions)
            