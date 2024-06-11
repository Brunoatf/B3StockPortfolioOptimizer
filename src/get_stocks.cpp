#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <set>
#include <map>
#include <vector>

using namespace std;

int main() {
	ifstream input_file;
	input_file.open("dados_totais.csv");
	map <string, string> stock_first_dates;
	while(!input_file.eof()) {
		string date, open, high, low, close, adj_close, volume, name;
		getline(input_file, date, ',');
		getline(input_file, open, ',');
		getline(input_file, high, ',');
		getline(input_file, low, ',');
		getline(input_file, close, ',');
		getline(input_file, adj_close, ',');
		getline(input_file, volume, ',');
		getline(input_file, name);
		if(stock_first_dates.find(name) == stock_first_dates.end()) stock_first_dates[name] = date;
		else {
			stock_first_dates[name] = min(stock_first_dates[name], date);
		}
	}
	vector <string> stock_names;
	vector <vector <double> > stock_values, after_period;
	for(auto cur : stock_first_dates) {
		if(cur.second <= "2009-05-07") {
			if(cur.second == "") continue;
			//cout << cur.first << ' ' << cur.second << '\n';
			stock_names.push_back(cur.first);
			vector <double> temp;
			stock_values.push_back(temp);
			after_period.push_back(temp);
		}
	}
	ifstream get_data;
	get_data.open("dados_totais.csv");
	while(!get_data.eof()) {
		string date, open, high, low, close, adj_close, volume, name;
		getline(get_data, date, ',');
		getline(get_data, open, ',');
		getline(get_data, high, ',');
		getline(get_data, low, ',');
		getline(get_data, close, ',');
		getline(get_data, adj_close, ',');
		getline(get_data, volume, ',');
		getline(get_data, name);
		if(date <= "2009-05-07") continue;
		if(name == "") continue;
		for(int i = 0; i < stock_names.size(); i++) {
			if(name == stock_names[i]) {
				if(date <= "2019-05-07") {
					if(close == "") {
						stock_values[i].push_back(stock_values[i][stock_values[i].size() - 1]);
					}
					else stock_values[i].push_back(stod(close));
				}
				else {
					if(close == "") continue;
					else after_period[i].push_back(stod(close));
				}
			}
		}
	}
	vector <double> final_return, expected_returns;
	for(int i = 0; i < stock_names.size(); i++) {
		expected_returns.push_back(stock_values[i][stock_values[i].size() - 1] / stock_values[i][0]);
		final_return.push_back(after_period[i][after_period[i].size() - 1] / after_period[i][0]);
	}
	vector <vector <double> > covariance;
	for(int i = 0; i < stock_names.size(); i++) {
		vector <double> temp;
		covariance.push_back(temp);
		for(int j = 0; j < stock_names.size(); j++) {
			while(stock_values[i].size() > stock_values[j].size()) stock_values[i].pop_back();
			while(stock_values[j].size() > stock_values[i].size()) stock_values[j].pop_back();
			double mean1 = 0, mean2 = 0;
			for(double cur : stock_values[i]) mean1 += cur;
			for(double cur : stock_values[j]) mean2 += cur;
			mean1 /= stock_values[i].size(); mean2 /= stock_values[j].size();
			covariance[i].push_back(0);
			for(int k = 0; k < stock_values[i].size(); k++) {
				covariance[i][j] += (stock_values[i][k] - mean1) * (stock_values[j][k] - mean2);
			}
			covariance[i][j] /= (stock_values[i].size() - 1);
		}
	}
	cout << stock_names.size() << '\n';
	for(int i = 0; i < stock_names.size(); i++) cout << expected_returns[i] << ' ';
	cout << '\n';
	for(int i = 0; i < stock_names.size(); i++) cout << final_return[i] << ' ';
	cout << '\n';
	for(int i = 0; i < stock_names.size(); i++) {
		for(int j = 0; j < stock_names.size(); j++) {
			cout << covariance[i][j] << ' ';
		}
		cout << '\n';
	} 
	return 0;
}