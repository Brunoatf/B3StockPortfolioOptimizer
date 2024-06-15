import numpy as np
import matplotlib.pyplot as plt

def read_performance_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    num_algorithms = len(lines) - 1  # Last line is the optimal value
    num_problems = 50

    performance_data = np.zeros((num_algorithms, num_problems))

    for i in range(num_algorithms):
        performance_data[i, :] = list(map(float, lines[i].split()))

    optimal_value = float(lines[-1].strip())
    
    return performance_data, optimal_value

file_path = 'performance_profile.txt'

performance_data, optimal_value = read_performance_data(file_path)

num_algorithms = performance_data.shape[0]
num_problems = performance_data.shape[1]

performance_ratios = np.zeros_like(performance_data)
for i in range(num_algorithms):
    performance_ratios[i, :] =  optimal_value / performance_data[i, :]

taus = np.linspace(1, np.max(performance_ratios), 100)
cdf_values = np.zeros((num_algorithms, len(taus)))

for j in range(num_algorithms):
    for k, tau in enumerate(taus):
        cdf_values[j, k] = np.mean(performance_ratios[j, :] <= tau)

plt.figure(figsize=(10, 6))
for j in range(num_algorithms):
    name = ""
    if j % 3 == 0 :
        name += "Uniform "
    elif j % 3 == 1 :
        name += "Flat "
    else :
        name += "Blend "
    if j < 3 :
        name += "tournament"
    else :
        name += "ranking"
    plt.plot(taus, cdf_values[j], label=name)

plt.xlabel('Desempenho')
plt.ylabel('Probabilidade')
plt.title('Performance Profile Plot')
plt.legend()
plt.grid(True)
plt.show()
