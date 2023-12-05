class Phone:
    def __init__(self, os, number, is_WaterProof):
        self.os = os
        self.number = number
        self.isWaterProof = is_WaterProof

    def is_ios(self):
        if self.os == "ios":
            return True
        else:
            return False

    def add(self, num1, num2):
        return num1 + num2


phone1 = Phone("ios", 123, True)
phone2 = Phone("android", 456, False)

print(phone1.os)
print(phone1.is_ios())
print(phone1.add(3, 4))
