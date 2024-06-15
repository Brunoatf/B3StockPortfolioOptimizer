#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <iomanip>

using namespace std;

double target = 1.45;

int main() {
	double melhor_rendimento = 0;
	for(int ab = 0; ab < 6; ab++) {
		double rendimento_medio = 0;
		vector <int> uni_tour_ttt;
		double uni_tour_medio[53];
		for(int i = 0; i < 53; i++) uni_tour_medio[i] = 0;
		string s;
		for(int i = 0; i < 50; i++) {
			bool ok = false;
			while(cin >> s, s == "Generation") {
				int gen;
				cin >> gen;
				cin >> s;
				cin >> s;
				cin >> s;
				double cur;
				cin >> cur;
				if(cur >= target && !ok) {
					ok = true;
					uni_tour_ttt.push_back(gen);
				}
			}
			for(int j = 0; j < 53; j++) {
				int cur;
				cin >> cur;
			}
			for(int j = 0; j < 53; j++) {
				double cur;
				cin >> cur;
				uni_tour_medio[j] += cur;
			}
			cin >> s;
			cin >> s;
			double cur;
			cin >> cur;
			melhor_rendimento = max(melhor_rendimento, cur);
			rendimento_medio += cur;
			cout << cur << ' ';
		}
		/*cout << uni_tour_ttt.size() << '\n';
		for(int cur : uni_tour_ttt) {
			cout << cur << ' ';
		}
		cout << '\n';*/
		//cout << rendimento_medio / 50 << '\n';
		/*cout << 53 << '\n';
		for(int i = 0; i < 53; i++) {
			cout << uni_tour_medio[i] / 50;
		}*/
		cout << '\n';
	}
	cout << melhor_rendimento << '\n';
	return 0;
}