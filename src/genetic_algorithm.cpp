#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <iomanip>

using namespace std;

const double ALPHA = 0.2;

/**
 * @class Chromosome
 * @brief Classe que representa um cromossomo.
 */
class Chromosome {

    public:

        int chromosome_length;
        double minimum_stock_proportion;
        double maximum_stock_proportion;
        int portfolio_size;
        double fitness;
        vector<int> selected_stocks;
        vector<double> stock_proportions;

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
            double minimum_stock_proportion,
            double maximum_stock_proportion,
            const vector<double> &expected_returns,
            const vector<vector<double>> &cov_matrix,
            double risk_tolerance_coefficient,
            vector<int> selected_stocks = {},
            vector<double> stock_proportions = {}
        ) : 
            chromosome_length(chromosome_length),
            minimum_stock_proportion(minimum_stock_proportion),
            maximum_stock_proportion(maximum_stock_proportion),
            portfolio_size(portfolio_size)
        {   

            random_device rd;
            mt19937 gen(rd());
            bernoulli_distribution dis(0.5);

            if (selected_stocks.empty()) {
                selected_stocks.resize(chromosome_length);
                for (int &stock : selected_stocks) {
                    stock = dis(gen);
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

            fitness = compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient);
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
                mt19937 gen(random_device{}());
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
                vector <pair <double, int> > sorted_proportions;
                for(int i = 0; i < chromosome_length; ++i) {
                    sorted_proportions.push_back({stock_proportions[i], i});
                }
                sort(sorted_proportions.begin(), sorted_proportions.end());

                for (int i = 0; i < selected_stocks_count - portfolio_size; ++i) {
                    if(selected_stocks[sorted_proportions[i].second] == 0) {
                        selected_stocks_count++;
                        continue;
                    }
                    selected_stocks[sorted_proportions[i].second] = 0;
                    stock_proportions[sorted_proportions[i].second] = 0.0;
                }
                //cout << selected_stocks_count << '\n';
                int qt = 0;
                for (int i = 0; i < chromosome_length; ++i) {
                    if (selected_stocks[i] == 1) {
                        qt++;
                    }
                }
                //if(qt != portfolio_size) cout << "oi " << qt << '\n';
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
                set<int> initial_chosen_stocks_copy = initial_chosen_stocks;
                for (auto stock : initial_chosen_stocks_copy) {
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
         * @param expected_returns Vetor de inteiros com os retornos esperados de cada ação.
         * @param cov_matrix Matriz de covariância dos retornos das ações.
         * @param risk_tolerance_coefficient Coeficiente de tolerância ao risco.
         * @return O valor do fitness do cromossomo.
         */
        double compute_fitness(
            vector<double> expected_returns,
            vector<vector<double>> cov_matrix,
            double risk_tolerance_coefficient //lambda na fórmula
        ) {
            repair(portfolio_size);

            // calcular o fitness com base nos dados históricos de cada ação
            // o fitness será (1 - lambda) * retorno do portfolio - lambda * risco do portfolio
            // onde lambda é o coeficiente de tolerância ao risco
            // o retorno do portfolio é a média dos retornos de cada ação ponderada pela proporção de cada ação no portfolio
            // o risco do portfolio é a soma das covariâncias dos retornos de cada par de ações ponderada pela proporção de cada ação no portfolio

            double return_sum = 0.0;
            double risk_sum = 0.0;

            for (int i = 0; i < chromosome_length; i++) {
                if (selected_stocks[i]) {
                    return_sum += stock_proportions[i] * expected_returns[i];
                    for (int j = 0; j < chromosome_length; j++) {
                        if (selected_stocks[j]) {
                            risk_sum += stock_proportions[i] * stock_proportions[j] * cov_matrix[i][j];
                        }
                    }
                }
            }

            return risk_tolerance_coefficient * return_sum - (1 - risk_tolerance_coefficient) * risk_sum;
        }

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
        * @param risk_tolerance_coefficient Coeficiente de tolerância ao risco.
        * @param minimum_stock_proportion Proporção mínima de uma ação no portfólio.
        * @param maximum_stock_proportion Proporção máxima de uma ação no portfólio.
        * @param expected_returns Vetor de inteiros com os retornos esperados de cada ação.
        * @param cov_matrix Matriz de covariância dos retornos das ações.
        */
        GeneticAlgorithm(
            int population_size,
            double mutation_rate,
            double crossover_rate,
            int elitism_count,
            int tournament_size,
            int portfolio_size,
            int num_available_stocks,
            double risk_tolerance_coefficient,
            double minimum_stock_proportion,
            double maximum_stock_proportion,
            vector<double> expected_returns,
            vector<vector<double>> cov_matrix,
            string crossover_type = "uniform",
            string selection_type = "tournament"
        ) : population_size(population_size),
            mutation_rate(mutation_rate),
            crossover_rate(crossover_rate),
            crossover_type(crossover_type),
            elitism_count(elitism_count),
            tournament_size(tournament_size),
            portfolio_size(portfolio_size),
            num_available_stocks(num_available_stocks),
            risk_tolerance_coefficient(risk_tolerance_coefficient),
            minimum_stock_proportion(minimum_stock_proportion),
            maximum_stock_proportion(maximum_stock_proportion),
            expected_returns(expected_returns),
            cov_matrix(cov_matrix),
            selection_type(selection_type)
        {}

        /**
        * @brief Sobrecarga do operador de chamada de função para executar o algoritmo genético.
        * @param num_generations Número de gerações a serem executadas.
        */
        Chromosome operator()(int num_generations) {
            vector<Chromosome> population = init_population();

            //cout << "Starting iterations" << endl;
            Chromosome solution = get_best_chromosome(population);

            for (int i = 0; i < num_generations; ++i) {

                auto parents = select_chromosomes(population);

                population = crossover_population(parents);

                mutate_population(population);

                population = repair_population(population);

                double population_fitness = eval_population(population);

                auto best_chromosome = get_best_chromosome(population);
                if (best_chromosome.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) > solution.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient)) {
                    solution = best_chromosome;
                    cout << "Generation " << i << " best chromosome fitness: " << best_chromosome.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) << endl;
                }
            
                /*if (i % 10 == 0) {
                    cout << "Generation " << i << " best chromosome fitness: " << best_chromosome.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) << endl;
                }*/
            
            }

            return solution;
        }

    private:
        int population_size;
        double mutation_rate;
        double crossover_rate;
        int elitism_count;
        int tournament_size;
        int portfolio_size;
        int num_available_stocks;
        double risk_tolerance_coefficient;
        double minimum_stock_proportion;
        double maximum_stock_proportion;
        vector<double> expected_returns;
        vector<vector<double>> cov_matrix;
        string crossover_type;
        string selection_type;

        /**
        * @brief Inicializa a população de cromossomos.
        * @return Um vetor de cromossomos representando a população inicial.
        */
        vector<Chromosome> init_population() {
            vector<Chromosome> population;
            for (int i = 0; i < population_size; i++) {
                Chromosome chromosome(num_available_stocks, portfolio_size, minimum_stock_proportion, maximum_stock_proportion, expected_returns, cov_matrix, risk_tolerance_coefficient);
                population.push_back(chromosome);
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
                population_fitness += chromosome.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient);
            }
            return population_fitness / population_size;
        }

        /**
        * @brief Retorna o melhor cromossomo da população.
        * @param population Um vetor de cromossomos representando a população.
        * @return O cromossomo com o melhor fitness.
        */
        Chromosome get_best_chromosome(std::vector<Chromosome>& population) {
            Chromosome best = population[0];
            for (int i = 1; i < population.size(); i++) {
                if (population[i].compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) > best.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient)) {
                    best = population[i];
                }
            }
            return best;
        }

        /**
        * @brief Realiza o crossover na população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos resultantes do crossover.
        */
        vector<Chromosome> crossover_population(vector<Chromosome>& population) {
            
            vector<Chromosome> new_population;
            for (int i = 0; i < population_size - 1; i = i + 2) {
                Chromosome parent1 = population[i];
                Chromosome parent2 = population[i+1];
                if (crossover_type == "uniform") {
                    bernoulli_distribution b(0.5);
                    mt19937 generator(random_device{}());
                    int value;
                    for (int stock = 0; stock < num_available_stocks; stock++) {
                        value = b(generator);

                        // Se obtivermos 1, então ocorre a troca. Senão, mantém-se como está.
                        if (value == 1) {
                            swap(parent1.selected_stocks[stock], parent2.selected_stocks[stock]);
                            swap(parent1.stock_proportions[stock], parent2.stock_proportions[stock]);
                        }
                    }
                }
                else if(crossover_type == "flat") {
                    mt19937 generator(random_device{}());
                    for(int stock = 0; stock < num_available_stocks; stock++) {
                        double lower_bound = min(parent1.stock_proportions[stock], parent2.stock_proportions[stock]);
                        double upper_bound = max(parent1.stock_proportions[stock], parent2.stock_proportions[stock]);
                        lower_bound = max(lower_bound, 0.0);
                        uniform_real_distribution<> d(lower_bound, upper_bound);
                        double value1 = d(generator);
                        double value2 = d(generator);
                        if(value1 > 0) parent1.selected_stocks[stock] = 1;
                        if(value2 > 0) parent2.selected_stocks[stock] = 1;
                        parent1.stock_proportions[stock] = value1;
                        parent2.stock_proportions[stock] = value2;
                    }
                }
                else if(crossover_type == "blend") {
                    mt19937 generator(random_device{}());
                    for(int stock = 0; stock < num_available_stocks; stock++) {
                        double blend_interval_size = ALPHA * abs(parent1.stock_proportions[stock] - parent2.stock_proportions[stock]);
                        double lower_bound = min(parent1.stock_proportions[stock], parent2.stock_proportions[stock]) - blend_interval_size;
                        double upper_bound = max(parent1.stock_proportions[stock], parent2.stock_proportions[stock]) + blend_interval_size;
                        lower_bound = max(lower_bound, 0.0);
                        uniform_real_distribution<> d(lower_bound, upper_bound);
                        double value1 = d(generator);
                        double value2 = d(generator);
                        if(value1 > 0) parent1.selected_stocks[stock] = 1;
                        if(value2 > 0) parent2.selected_stocks[stock] = 1;
                        parent1.stock_proportions[stock] = value1;
                        parent2.stock_proportions[stock] = value2;
                    }
                }
                new_population.push_back(parent1);
                new_population.push_back(parent2);
            }

            if (population_size % 2 == 1) {
                new_population.push_back(population[population_size-1]);
            }

            // Implementar outros crossovers
            // Implementados flat crossover e blend crossover

            return new_population;

        }

        /**
        * @brief Aplica mutação na população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos resultantes da mutação.
        */
        void mutate_population(vector<Chromosome>& population) {
            bernoulli_distribution b(mutation_rate);
            mt19937 generator(random_device{}());
            for (Chromosome &indiviual : population) {
                for (int i = 0; i < num_available_stocks; i++) {
                    if (b(generator)) {
                        indiviual.mutate_gene(i);
                    }
                }
            }
        }

        static bool decreasing_compare_func(Chromosome a, Chromosome b, vector<double> &expected_returns, vector<vector<double>> &cov_matrix, double &risk_tolerance_coefficient) {
            return a.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) > b.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient);
        }

        static bool compare_func(Chromosome a, Chromosome b, vector<double> &expected_returns, vector<vector<double>> &cov_matrix, double &risk_tolerance_coefficient) {
            return a.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient) < b.compute_fitness(expected_returns, cov_matrix, risk_tolerance_coefficient);
        }

        /**
        * @brief Seleciona cromossomos da população.
        * @param population Um vetor de cromossomos representando a população.
        * @return Um vetor de cromossomos selecionados.
        */
        vector<Chromosome> select_chromosomes(vector<Chromosome>& population) {

            if(selection_type == "tournament") {
                vector<Chromosome> parents;

                sort(
                    population.begin(),
                    population.end(),
                    [this](Chromosome a, Chromosome b) {
                        return decreasing_compare_func(a, b, expected_returns, cov_matrix, risk_tolerance_coefficient);
                    }
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
            else {
                vector <Chromosome> parents;
                sort(
                    population.begin(), 
                    population.end(),
                    [this](Chromosome a, Chromosome b) {
                        return compare_func(a, b, expected_returns, cov_matrix, risk_tolerance_coefficient);
                    }
                );

                vector <int> probabilities;

                for(int i = 0; i < population_size; i++) {
                    probabilities.push_back(i + 1);
                }

                for (int i = 0; i < elitism_count; ++i) {
                    parents.push_back(population[i]);
                }

                for(int i = 0; i < population_size - elitism_count; i++) {
                    random_device rd;
                    mt19937 gen(rd());
                    uniform_int_distribution<> distr(0, population_size * (population_size + 1) / 2 - 1);
                    int rdnum = distr(gen);
                    int cur = 0;
                    while(rdnum >= probabilities[cur]) {
                        rdnum -= probabilities[cur];
                        cur++;
                    }
                    if(cur == population_size) cur--;
                    parents.push_back(population[cur]);
                }
                return parents;
            }
        }
};

int main() {

    int num_stocks;
    cin >> num_stocks;

    vector<double> expected_returns, final_returns;
    vector<vector<double>> cov_matrix(num_stocks, vector<double>(num_stocks, 0));

    //montar expected_returns e cov_matrix com valores reais
    //para testes, usei valores aleatórios mesmo

    //random number generator between 0 and 1
    /*mt19937 generator(random_device{}());
    uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < NUM_STOCKS; i++) {
        for (int j = 0; j < NUM_STOCKS; j++) {
            cov_matrix[i][j] = distribution(generator);
        }
    }

    for (int i = 0; i < NUM_STOCKS; i++) {
        expected_returns.push_back(distribution(generator));
    }*/

    for(int i = 0; i < num_stocks; i++) {
        double cur;
        cin >> cur;
        expected_returns.push_back(cur);
    }
    for(int i = 0; i < num_stocks; i++) {
        double cur;
        cin >> cur;
        final_returns.push_back(cur);
    }
    for(int i = 0; i < num_stocks; i++) {
        for(int j = 0; j < num_stocks; j++) {
            double cur;
            cin >> cur;
            cov_matrix[i][j] = cur;
        }
    }

    /**
    * @brief Construtor que inicializa os parâmetros do algoritmo genético.
    * @param population_size Tamanho da população.
    * @param mutation_rate Taxa de mutação.
    * @param crossover_rate Taxa de crossover.
    * @param elitism_count Número de cromossomos de elitismo.
    * @param tournament_size Tamanho do torneio.
    * @param portfolio_size Tamanho do portfólio.
    * @param num_available_stocks Número de ações disponíveis.
    * @param risk_tolerance_coefficient Coeficiente de tolerância ao risco.
    * @param minimum_stock_proportion Proporção mínima de uma ação no portfólio.
    * @param maximum_stock_proportion Proporção máxima de uma ação no portfólio.
    * @param expected_returns Vetor de inteiros com os retornos esperados de cada ação.
    * @param cov_matrix Matriz de covariância dos retornos das ações.
    */
    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "uniform",
            "tournament"
        );


        Chromosome solution = ga(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution.stock_proportions[i] * final_returns[i] * solution.selected_stocks[i];
        }
        cout << "Uniform\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) <<solution.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "uniform tournament: "<< real_return << '\n';
    }
    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga_flat(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "flat",
            "tournament"
        );

        Chromosome solution_flat = ga_flat(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution_flat.stock_proportions[i] * final_returns[i] * solution_flat.selected_stocks[i];
        }
        cout << "Flat\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution_flat.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) << solution_flat.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "flat tournament: "<< real_return << '\n';
    }
    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga_blend(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "blend",
            "tournament"
        );


        Chromosome solution_blend = ga_blend(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution_blend.stock_proportions[i] * final_returns[i] * solution_blend.selected_stocks[i];
        }
        cout << "Blend\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution_blend.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) << solution_blend.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "blend tournament: "<< real_return << '\n';
    }

    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "uniform",
            "ranking"
        );


        Chromosome solution = ga(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution.stock_proportions[i] * final_returns[i] * solution.selected_stocks[i];
        }
        cout << "Uniform\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) <<solution.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "uniform ranking: "<< real_return << '\n';
    }
    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga_flat(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "flat",
            "ranking"
        );

        Chromosome solution_flat = ga_flat(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution_flat.stock_proportions[i] * final_returns[i] * solution_flat.selected_stocks[i];
        }
        cout << "Flat\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution_flat.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) << solution_flat.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "flat ranking: "<< real_return << '\n';
    }
    for(int i = 0; i < 50; i++) {
        GeneticAlgorithm ga_blend(
            100,
            0.005,
            0.1,
            20,
            3,
            20,
            num_stocks,
            0.5,
            0.02,
            0.5,
            expected_returns,
            cov_matrix,
            "blend",
            "ranking"
        );


        Chromosome solution_blend = ga_blend(1000);

        double real_return = 0;
        for(int i = 0; i < num_stocks; i++) {
            real_return += solution_blend.stock_proportions[i] * final_returns[i] * solution_blend.selected_stocks[i];
        }
        cout << "Blend\n";
        for(int i = 0; i < num_stocks; i++) {
            cout << solution_blend.selected_stocks[i] << ' ';
        }
        cout << '\n';
        for(int i = 0; i < num_stocks; i++) {
            cout << fixed << setprecision(10) << solution_blend.stock_proportions[i] << ' ';
        }
        cout << '\n';
        cout << "blend ranking: "<< real_return << '\n';
    }

    double standard_return = 0;
    for(int i = 0; i < num_stocks; i++) {
        standard_return += final_returns[i] / num_stocks;
    }


    cout << "Equal\n";
    for(int i = 0; i < num_stocks; i++) {
        cout << 1 << ' ';
    }
    cout << '\n';
    for(int i = 0; i < num_stocks; i++) {
        cout << fixed << setprecision(10) << 1/num_stocks << ' ';
    }
    cout << '\n';
    cout << "equal: " << standard_return << '\n';
    return 0;
}
