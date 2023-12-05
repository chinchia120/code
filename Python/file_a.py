file = open("C:\\Users\\user\\Documents\\code\\python\\123.txt", mode = "a", encoding = "utf-8")

file.close()

with open("C:\\Users\\user\\Documents\\code\\python\\123.txt", mode = "a", encoding = "utf-8") as file:
    file.write(" AAA")