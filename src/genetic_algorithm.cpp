#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <set>

using namespace std;


/**
 * @class Chromosome
 * @brief Classe que representa um cromossomo.
 */
class Chromosome {

public:

    /**
     * @brief Construtor que inicializa o cromossomo.
     * @param num_available_stocks Número de ações disponíveis.
     * @param portfolio_size Tamanho do portfólio.
     * @param minimum_stock_proportion Proporção mínima de uma ação no portfólio.
     * @param maximum_stock_proportion Proporção máxima de uma ação no portfólio.
     * @param selected_stocks Vetor de inteiros com 1 indicando que uma ação foi escolhida e 0 indicado o contrário
     * @param stock_proportions Vetor de doubles com cada elemento i indicando a porcentagem da ação i no portfólio
     */
    Chromosome(
        int chromosome_length,
        int portfolio_size,
        double minimum_stock_proportion = 0.05,
        double maximum_stock_proportion = 0.5,
        vector<int> selected_stocks = {},
        vector<double> stock_proportions = {}
    ) : 
        chromosome_length(chromosome_length),
        minimum_stock_proportion(minimum_stock_proportion),
        maximum_stock_proportion(maximum_stock_proportion),
        portfolio_size(portfolio_size) 
    {

        if (selected_stocks.empty()) {
            selected_stocks.resize(chromosome_length);
            for (int &stock : selected_stocks) {
                stock = rand() % 2;
            }
        }
        this->selected_stocks = selected_stocks;

        if (stock_proportions.empty()) {
            stock_proportions.resize(chromosome_length);
            for (double &proportion : stock_proportions) {
                proportion = static_cast<double>(rand()) / RAND_MAX;
            }
        }
        this->stock_proportions = stock_proportions;

        // Inicialização contando apenas soluções válidas:
        // Talvez seja interessante inicializar com soluções inválidas para explorar o espaço de busca
        repair(portfolio_size);

        fitness = compute_fitness();
    }

    /**
     * @brief Obtém o comprimento do cromossomo.
     * @return O comprimento do cromossomo.
     */
    int length() const {
        return chromosome_length;
    }
    
    /*
    * @brief Faz a mutação de um gene do cromossomo
    * @param gene O gene que deseja-se aplicar uma mutação
    */
    void mutate_gene(int gene) { 
        if (selected_stocks[gene]) {
            selected_stocks[gene] = 0;
        } 
        else {
            uniform_real_distribution<> dis(0.0, 1.0);
            mt19937 gen(random_device());
            selected_stocks[gene] = 1;
            stock_proportions[gene] = dis(gen);
        }
    }

    /**
     * @brief Repara o cromossomo para garantir que ele é válido.
     * @param portfolio_size Tamanho do portfólio.
     */
    void repair(int portfolio_size) {

        // Primeiramente ajustar o tamanho do portfolio:

        int selected_stocks_count = accumulate(selected_stocks.begin(), selected_stocks.end(), 0);

        // Se o número de ações selecionadas for maior que o tamanho do portfolio, remover as ações de menor proporção que já foram selecionadas:
        if (selected_stocks_count > portfolio_size) {
            auto sorted_proportions = stock_proportions;
            sort(sorted_proportions.begin(), sorted_proportions.end());

            for (int i = 0; i < selected_stocks_count - portfolio_size; ++i) {
                auto it = find(stock_proportions.begin(), stock_proportions.end(), sorted_proportions[i]);
                int stock_index = distance(stock_proportions.begin(), it);
                selected_stocks[stock_index] = 0;
                stock_proportions[stock_index] = 0.0;
            }
        }
        // Se o número de ações selecionadas for menor que o tamanho do portfolio, adicionar as ações de maior proporção que ainda não foram selecionadas:
        // Aqui talvez seja interessante permitir que o algoritmo monte portfólios com menos ações. Vai que o melhor portfólio não tem exatamente o tamanho do portfolio_size
        else if (selected_stocks_count < portfolio_size) {
            auto sorted_proportions = stock_proportions;
            sort(sorted_proportions.rbegin(), sorted_proportions.rend());

            for (int i = 0; i < portfolio_size - selected_stocks_count; ++i) {
                auto it = find(stock_proportions.begin(), stock_proportions.end(), sorted_proportions[i]);
                int stock_index = distance(stock_proportions.begin(), it);
                selected_stocks[stock_index] = 1;
            }
        }

        // Em seguida ajustar as proporções das ações selecionadas:

        // Ações não escolhidas recebem proporção 0:
        for (int i = 0; i < chromosome_length; ++i) {
            if (selected_stocks[i] == 0) {
                stock_proportions[i] = 0.0;
            }
        }

        // Aplica-se a restrição de proporção mínima:
        double unnormalized_weight_sum = accumulate(stock_proportions.begin(), stock_proportions.end(), 0.0);
        double minimum_new_weight_sum = minimum_stock_proportion * portfolio_size;
        for (int i = 0; i < chromosome_length; ++i) {
            if (selected_stocks[i] == 1) {
                // Essa fórmula garante que a soma das proporções das ações selecionadas seja igual a 1:
                // Só funciona se o portfolio_size for igual ao número de ações selecionadas
                // A ideia é começar com a proporção mínima e distribuir o restante proporcionalmente. Formulação melhor está no artigo related work 1
                stock_proportions[i] = minimum_stock_proportion + (stock_proportions[i] / unnormalized_weight_sum) * (1 - minimum_new_weight_sum);
            }
        }

        // Restrição de proporção máxima descrita no arquivo upper_bound_reference:
        set<int> over_bound_stocks;
        set<int> initial_chosen_stocks;
        for (int i = 0; i < selected_stocks.size(); ++i) {
            if (selected_stocks[i] == 1) {
                initial_chosen_stocks.insert(i);
            }
        }

        // Enquanto houver alguma ação em initial_chosen_stocks que esteja acima do limite de proporção, ajustar as proporções:
        while (true) {
            int remaining_over_bound = 0;
            for (auto stock : initial_chosen_stocks) {
                if (stock_proportions[stock] > maximum_stock_proportion) {
                    ++remaining_over_bound;
                }
            }
            if (remaining_over_bound == 0) break;

            // calcula a porcentagem total ainda presente em initial_chosen_stocks
            // ações com porcentagem acima do permitido passam pra over_bound_stocks com proporção corrigida para a máxima, 1 por iteração
            double total_percentage = 0.0;
            bool moved_to_upper_bound = false;
            for (auto stock : initial_chosen_stocks) {
                total_percentage += stock_proportions[stock];
                if (!moved_to_upper_bound && stock_proportions[stock] > maximum_stock_proportion) {
                    moved_to_upper_bound = true;
                    over_bound_stocks.insert(stock);
                    stock_proportions[stock] = maximum_stock_proportion;
                    initial_chosen_stocks.erase(stock);
                }
            }

            // Calcula o quanto daria pra distribuir tomando que todas as stocks têm pelo menos a proporção mínima
            // E também que todas as em over_bound_stocks têm a proporção máxima
            double free_percentage = 1 - (initial_chosen_stocks.size() * minimum_stock_proportion + maximum_stock_proportion * over_bound_stocks.size());

            for (auto stock : initial_chosen_stocks) {
                stock_proportions[stock] = minimum_stock_proportion + (stock_proportions[stock] / total_percentage) * free_percentage;
            }
        }
    }

    /**
     * @brief Calcula e retorna o fitness do cromossomo.
     * @return O valor do fitness do cromossomo.
     */
    double compute_fitness() {
        repair(portfolio_size);

        // calcular o fitness com base nos dados históricos de cada ação
        // o fitness será (1 - lambda) * retorno do portfolio - lambda * risco do portfolio
        // onde lambda é o coeficiente de tolerância ao risco
        // o retorno do portfolio é a média dos retornos de cada ação ponderada pela proporção de cada ação no portfolio
        // o risco do portfolio é a soma das covariâncias dos retornos de cada par de ações ponderada pela proporção de cada ação no portfolio
        return 0.0;
    }

    /**
     * @brief Calcula o retorno do portfólio.
     * @return O valor do retorno do portfólio.
     */
    double compute_portfolio_return() {
        // Implementar o cálculo do retorno do portfólio
        return 0.0;
    }

    /**
     * @brief Calcula o risco do portfólio.
     * @return O valor do risco do portfólio.
     */
    double compute_portfolio_risk() {
        // Implementar o cálculo do risco do portfólio
        return 0.0;
    }

private:
    int chromosome_length;
    double minimum_stock_proportion;
    double maximum_stock_proportion;
    int portfolio_size;
    double fitness;
    vector<int> selected_stocks;
    vector<double> stock_proportions;
};

/**
 * @class GeneticAlgorithm
 * @brief Classe que implementa um algoritmo genético para o problema da otimização de portfólio.
 */
class GeneticAlgorithm {
    public:
        /**
        * @brief Construtor que inicializa os parâmetros do algoritmo genético.
        * @param population_size Tamanho da população.
        * @param mutation_rate Taxa de mutação.
        * @param crossover_rate Taxa de crossover.
        * @param elitism_count Número de cromossomos de elitismo.
        * @param tournament_size Tamanho do torneio.
        * @param portfolio_size Tamanho do portfólio.
        * @param num_available_stocks Número de ações disponíveis.
        * @param risk_tolerance_coefficent Coeficiente de tolerância ao risco.
        * @param minimum_stock_proportion Proporção mínima de uma ação no portfólio.
        * @param maximum_stock_proportion Proporção máxima de uma ação no portfólio.
        */
        GeneticAlgorithm(
            int population_size,
            double mutation_rate,
            double crossover_rate,
            int elitism_count,
            int tournament_size,
            int portfolio_size,
            int num_available_stocks,
            double risk_tolerance_coefficent,
            double minimum_stock_proportion,
            double maximum_stock_proportion
        ) : population_size(population_size),
            mutation_rate(mutation_rate),
            crossover_rate(crossover_rate),
            elitism_count(elitism_count),
            tournament_size(tournament_size),
            portfolio_size(portfolio_size),
            num_available_stocks(num_available_stocks),
            risk_tolerance_coefficent(risk_tolerance_coefficent),
            minimum_stock_proportion(minimum_stock_proportion),
            maximum_stock_proportion(maximum_stock_proportion)
        {}

        /**
        * @brief Sobrecarga do operador de chamada de função para executar o algoritmo genético.
        * @param num_generations Número de gerações a serem executadas.
        */
        void operator()(int num_generations) {
            vector<Chromosome> population = init_population();

            for (int i = 0; i < num_generations; ++i) {

                auto parents = select_chromosomes(population);
                population = crossover_population(parents);
                population = mutate_population(population);
                population = repair_population(population);

                double population_fitness = eval_population(population);
                cout << "Generation " << i << " mean fitness: " << population_fitness << endl;

                auto best_chromosome = get_best_chromosome(population);
                cout << "Generation " << i << " best chromosome fitness: " << best_chromosome.compute_fitness() << endl;
            }
        }

    private:
        int population_size;
        double mutation_rate;
        double crossover_rate;
        int elitism_count;
        int tournament_size;
        int portfolio_size;
        int num_available_stocks;
        double risk_tolerance_coefficent;
        double minimum_stock_proportion;
        double maximum_stock_proportion;

        /**
        * @brief Inicializa a população de cromossomos.
        * @return Um vetor de cromossomos representando a população inicial.
        */
        vector<Chromosome> init_population() {
            vector<Chromosome> population;
            for (int i = 0; i < population_size; ++i) {
                population.emplace_back(num_available_stocks, portfolio_size, minimum_stock_proportion, maximum_stock_proportion);
            }
            return population;
        }
        
        /**
        * @brief Repara a população de cromossomos.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos reparados.
        */
        vector<Chromosome> repair_population(vector<Chromosome>& population) {
            for (auto& chromosome : population) {
                chromosome.repair(portfolio_size);
            }
            return population;
        }

        /**
        * @brief Avalia a população de cromossomos.
        * @param population Um vetor de cromossomos representando a população.
        * @return O valor médio de fitness da população.
        */
        double eval_population(vector<Chromosome>& population) {
            double population_fitness = 0.0;
            for (auto& chromosome : population) {
                population_fitness += chromosome.compute_fitness();
            }
            return population_fitness / population_size;
        }

        /**
        * @brief Retorna o melhor cromossomo da população.
        * @param population Um vetor de cromossomos representando a população.
        * @return O cromossomo com o melhor fitness.
        */
        Chromosome get_best_chromosome(vector<Chromosome>& population) {
            auto best_it = max_element(population.begin(), population.end(),
                [](Chromosome& a, Chromosome& b) { return a.compute_fitness() < b.compute_fitness(); });
            return *best_it;
        }

        /**
        * @brief Realiza o crossover na população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos resultantes do crossover.
        */
        vector<Chromosome> crossover_population(vector<Chromosome>& population) {
            // Implementar a lógica de crossover
            return population;
        }

        /**
        * @brief Aplica mutação na população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos resultantes da mutação.
        */
        vector<Chromosome> mutate_population(vector<Chromosome>& population) {
            bernoulli_distribution b(mutation_rate);
            mt19937 generator(random_device());
            for (Chromosome &indiviual : population) {
                for (int i = 0; i < num_available_stocks; i++) {
                    if (b(generator)) {
                        indiviual.mutate_gene(i);
                    }
                }
            }
        }

        /**
        * @brief Seleciona cromossomos da população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos selecionados.
        */
        vector<Chromosome> select_chromosomes(vector<Chromosome>& population) {

            vector<Chromosome> parents;

            sort(
                population.begin(),
                population.end(),
                [](Chromosome& a, Chromosome& b) {return a.compute_fitness() > b.compute_fitness();}
            );

            for (int i = 0; i < elitism_count; ++i) {
                parents.push_back(population[i]);
            }
            
            for (int i = 0; i < population_size - elitism_count; i++) {
                shuffle(population.begin() + elitism_count, population.end(), mt19937{ random_device{}() });
                vector<Chromosome> tournament;
                for (int i = 0; i < tournament_size; i++) {
                    tournament.push_back(population[elitism_count + i]);
                }
                Chromosome winner = get_best_chromosome(tournament);
                parents.push_back(winner);
            }

            return parents;
        }
};
