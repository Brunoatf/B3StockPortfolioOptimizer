import numpy as np
import matplotlib.pyplot as plt

def read_reaching_values(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    sequences = []
    i = 0
    while i < len(lines):
        count = int(lines[i].strip())
        i += 1
        values = list(map(float, lines[i].strip().split()))
        values.sort()
        sequences.append(values)
        i += 1
    
    return sequences

file_path = 'uni_tour_ttt.txt'

sequences = read_reaching_values(file_path)

plt.figure(figsize=(10, 8))

for idx, seq in enumerate(sequences):
    n = len(seq)
    cumulative_proportion = (np.arange(1, n + 1, 1) - 1/2) / 50 
    print(cumulative_proportion)
    cumulative_time_on_test = seq
    print(cumulative_time_on_test)
    name = ""
    if idx % 3 == 0 :
        name += "Uniform "
    elif idx % 3 == 1 :
        name += "Flat "
    else :
        name += "Blend "
    if idx < 3 :
        name += "tournament"
    else :
        name += "ranking"
    plt.plot(cumulative_time_on_test, cumulative_proportion, marker='o', linestyle='-', label=name)

# Add labels and title
plt.xlabel('Total de ciclos')
plt.ylabel('Probabilidade de atingir o valor alvo')
plt.title('TTT Plot')

# Add legend
plt.legend()

# Show grid
plt.grid(True)

# Show the plot
plt.show()
