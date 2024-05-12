import random

class Chromosome():

    def __init__(self, chromosome_length, selected_stocks = None, stock_proportions = None):
        
        self.chromossome_length = chromosome_length

        if selected_stocks is None:
            self.selected_stocks = [1] * self.num_available_stocks + [0] * (self.portfolio_size - self.num_available_stocks)
            self.selected_stocks = random.shuffle(self.selected_stocks)

        if stock_proportions is None:
            self.stock_proportions = [random.random() for i in range(0, self.portfolio_size)]
            for i in range(self.population_size):
                if self.selected_stocks[i] == 0:
                    self.stock_proportions[i] = 0
            self.stock_proportions = [round(i / sum(self.stock_proportions)) for i in self.stock_proportions]

        self.fitness = self.compute_fitness()

    def __len__(self):
        return self.chromossome_length

    def compute_fitness(self):
        pass
        # calcular o fitness com base nos dados históricos de cada ação
        # o fitness será (1 - lambda) * retorno do portfolio - lambda * risco do portfolio
        # onde lambda é o coeficiente de tolerância ao risco
        # o retorno do portfolio é a média dos retornos de cada ação ponderada pela proporção de cada ação no portfolio
        # o risco do portfolio é a soma das covariâncias dos retornos de cada par de ações ponderada pela proporção de cada ação no portfolio
    
class GeneticAlgorithm():

    def __init__(
        self,
        population_size,
        mutation_rate,
        crossover_rate,
        elitism_count,
        tournament_size,
        portfolio_size,
        num_available_stocks,
        risk_tolerance_coefficent
    ):
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.elitism_count = elitism_count,
        self.tournament_size = tournament_size,
        self.portfolio_size = portfolio_size
        self.num_available_stocks = num_available_stocks
        self.risk_tolerance_coefficent = risk_tolerance_coefficent

    def eval_population(self, population):
        population_fitness = 0
        for i in range(0, self.population_size):
            population_fitness += population[i].fitness
        return population_fitness / self.population_size
    
    def get_best_chromosome(self, population):
        population.sort(key=lambda x: x.fitness, reverse=True)
        return population[0]

    def crossover_population(self, population):
        pass

    def mutate_population(self, population):
        pass

    def select_chromossomes(self, population):

        parents = []
        population.sort(key=lambda x: x.fitness, reverse=True)
    
        for i in range(0, self.elitism_count):
            parents.append(population[i])
            del population[i]

        for _ in range(self.tournament_size):
            population = random.shuffle(population)
            tournaments = [tournaments[i:i + self.tournament_size] for i in range(0, len(tournaments), self.tournament_size)]
            for tournament in tournaments:
                tournament.sort(key=lambda x: x.fitness, reverse=True)
                parents.append(tournament[0])
        
        return parents

    def __call__(self, num_generations):

        population = self.init_population()
        for i in range(0, num_generations):

            parents = self.select_chromossomes(population)

            population = self.crossover_population(parents)
            population = self.mutate_population(population)

            population_fitness = self.eval_population(population)
            print(f"Generation {i} mean fitness: {population_fitness}")

            best_chromosome = self.get_best_chromosome(population)
            print(f"Generation {i} best chromossome fitness: {best_chromosome.fitness}")
