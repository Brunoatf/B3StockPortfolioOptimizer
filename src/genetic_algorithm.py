import random

class Chromosome():

    def __init__(
        self,
        chromosome_length,
        portfolio_size,
        minimum_stock_proportion = 0.05,
        maximum_stock_proportion = 0.5,
        selected_stocks = None,
        stock_proportions = None
    ):
        
        self.chromossome_length = chromosome_length
        self.minimum_stock_proportion = minimum_stock_proportion
        self.maximum_stock_proportion = maximum_stock_proportion

        if selected_stocks is None:
            self.selected_stocks = [random.randint(0, 1) for i in range(self.chromossome_length)]
        else:
            self.selected_stocks = selected_stocks

        if stock_proportions is None:
            self.stock_proportions = [random.random() for i in range(self.chromossome_length)]
        else:
            self.stock_proportions = stock_proportions

        # Inicialização contando apenas soluções válidas:
        # Talvez seja interessante inicializar com soluções inválidas para explorar o espaço de busca
        self.repair(portfolio_size)

        self.fitness = self.compute_fitness()

    def __len__(self):
        return self.chromossome_length
    
    def repair(self, portfolio_size):

        # Primeiramente ajustar o tamanho do portfolio:

        selected_stocks_count = sum(self.selected_stocks)

        # Se o número de ações selecionadas for maior que o tamanho do portfolio, remover as ações de menor proporção que já foram selecionadas:
        if selected_stocks_count > portfolio_size:

            sorted_proportions = sorted(self.stock_proportions)

            for proportion in sorted_proportions[:selected_stocks_count - portfolio_size]:

                stock_index = self.stock_proportions.index(proportion)
                self.selected_stocks[stock_index] = 0
                self.stock_proportions[stock_index] = 0

        # Se o número de ações selecionadas for menor que o tamanho do portfolio, adicionar as ações de maior proporção que ainda não foram selecionadas:
        # Aqui talvez seja interessante permitir que o algoritmo monte portfólios com menos ações. Vai que o melhor portfólio não tem exatamente o tamanho do portfolio_size
        elif selected_stocks_count < portfolio_size:

            sorted_proportions = sorted(self.stock_proportions, reverse=True)

            for proportion in sorted_proportions[:portfolio_size - selected_stocks_count]:
                    
                    stock_index = self.stock_proportions.index(proportion)
                    self.selected_stocks[stock_index] = 1
                    self.stock_proportions[stock_index] = random.uniform(self.minimum_stock_proportion, self.maximum_stock_proportion)

        
        # Em seguida ajustar as proporções das ações selecionadas:

        # Ações não escolhidas recebem proporção 0:
        self.stock_proportions = [proportion * selected for proportion, selected in zip(self.stock_proportions, self.selected_stocks)]

        # Aplica-se a restrição de proporção mínima:
        unnormalized_weight_sum = sum(self.stock_proportions)
        minimum_new_weight_sum = self.minimum_stock_proportion * self.portfolio_size
        for i in range(self.chromossome_length):
            if self.selected_stocks[i] == 1:
                # Essa fórmula garante que a soma das proporções das ações selecionadas seja igual a 1:
                # Só functiona se o portfolio_size for igual ao número de ações selecionadas
                # A ideia é começar com a proporção mínima e distribuir o restante proporcionalmente. Formulação melhor está no artigo related work 1
                self.stock_proportions[i] = self.minimum_stock_proportion + (self.stock_proportions[i] / unnormalized_weight_sum) * (1 - minimum_new_weight_sum)
        
        # Restrição de proporção máxima descrita no arquivo upper_bound_reference:

        over_bound_stocks = set()
        initial_chosen_stocks = set()
        for i in range(len(self.selected_stocks)):
            if self.selected_stocks[i] == 1:
                initial_chosen_stocks.add(i)

        # Enqunato ouver alguma ação em initial_chosen_stocks que esteja acima do limite de proporção, ajustar as proporções:
        while remaining_over_bound := len([i for i in initial_chosen_stocks - over_bound_stocks if self.stock_proportions[i] > self.maximum_stock_proportion]) > 0:
            
            # calcula a porcentagem total ainda presente em inital_chosen_stocks
            # ações com porcentagem acima do permitido passam pra over_bound_stocks com propoção corrigida para a máxima, 1 por iteração
            total_percentage = 0
            moved_to_upper_bound = False
            for stock in initial_chosen_stocks - over_bound_stocks:
                total_percentage += self.stock_proportions[stock]
                if not moved_to_upper_bound and self.stock_proportions[stock] > self.maximum_stock_proportion:
                    moved_to_upper_bound = True
                    over_bound_stocks.add(stock)
                    self.stock_proportions[stock] = self.maximum_stock_proportion
                    initial_chosen_stocks.remove(stock)

            # Calcula o quanto daria pra distribuir tomando que todas as stocks tem pelo menos a proporção mínima
            # E também que todas as em over_bound_stocks tem a proporção máxima
            free_percentage = 1 - (len(initial_chosen_stocks - over_bound_stocks) * self.minimum_stock_proportion + self.maximum_stock_proportion * len(over_bound_stocks))

            for stock in initial_chosen_stocks - over_bound_stocks:
                self.stock_proportions[stock] = self.minimum_stock_proportion + (self.stock_proportions[stock] / total_percentage) * free_percentage

    def compute_portfolio_return(self):
        pass

    def compute_portfolio_risk(self):
        pass

    def compute_fitness(self):

        self.repair()

        fitness = (1 - self.risk_tolerance_coefficent) * self.compute_portfolio_return() - self.risk_tolerance_coefficent * self.compute_portfolio_risk()
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
        risk_tolerance_coefficent,
        minimum_stock_proportion,
        maximum_stock_proportion
    ):
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.elitism_count = elitism_count
        self.tournament_size = tournament_size
        self.portfolio_size = portfolio_size
        self.num_available_stocks = num_available_stocks
        self.risk_tolerance_coefficent = risk_tolerance_coefficent
        self.minimum_stock_proportion = minimum_stock_proportion
        self.maximum_stock_proportion = maximum_stock_proportion

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
            population = self.repair_population(population)

            population_fitness = self.eval_population(population)
            print(f"Generation {i} mean fitness: {population_fitness}")

            best_chromosome = self.get_best_chromosome(population)
            print(f"Generation {i} best chromossome fitness: {best_chromosome.fitness}")