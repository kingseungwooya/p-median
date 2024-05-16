import random


num_lines = 10


with open('numbers.txt', 'w') as file:
    for _ in range(num_lines):

        x = 9
        y = random.randint(0, 9)
        z = random.randint(0, 9)
        numbers = [x, y, z]

        for num in numbers:
            file.write(str(num) + ' ')
        file.write('\n')
