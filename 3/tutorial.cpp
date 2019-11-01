#include <vector>
#include <differential_evolution.hpp>
#include <iostream>
#include <random>
#include <iomanip>
#include <fstream>

using namespace std;

using namespace amichel::de;

double random(void) {
	static std::mt19937 generator(1);
	static std::uniform_real_distribution<double> distribution(0, 1);
	return distribution(generator);
}

double normalRandom(void) {
	static std::mt19937 generator(1);
	static normal_distribution<double> distribution(0, 1);
	return distribution(generator);
}

double random(double a, double b) {
	return a + random() * (b - a);
}

double calcPolynom(const vector<double>& coefs, double x) {
	double sum = 0;
	double xpow = 1;
	for (auto& i : coefs) {
		sum += xpow * i;
		xpow *= x;
	}
	return sum;
}

double distanceAbsolute(double a, double b) {
	return pow(a - b, 2);
}

double distanceRelative(double a, double b) {
	return distanceAbsolute(a, b) / fabs(max(a, b));
}

double distance(const vector<double>& coefs1, const vector<double>& coefs2) {
	double sum = 0;
	for (int i = 0; i < coefs1.size(); i++) {
		sum += distanceAbsolute(coefs1[i], coefs2[i]);
	}
	return sqrt(sum/coefs1.size());
}

double calcUniformMetric(double value, double x, double percent) {
	double r = fabs(value * percent * 2);
	if (fabs(x - value) < r) {
		return 1 + (1 - fabs(x - value)/r) * 9;
	} else {
		return r / fabs(x - value);
	}
}

double distance(const vector<pair<double, double>>& a, const vector<pair<double, double>>& b) {
	double percent = 0.2;
	double sum = 0;
	for (int i = 0; i < a.size(); i++) {
		sum += calcUniformMetric(a[i].second, b[i].second, percent);
	}
	return sum / a.size();
}

vector<double> toCoefs(amichel::de::DVectorPtr args) {
	vector<double> coefs;
	for (auto& i : *args) {
		coefs.push_back(i);
	}
	return coefs;
}

class sphere_function {
	vector<pair<double, double>> points;

 public:
  sphere_function(const vector<pair<double, double>>& points) : points(points) {
  }

  virtual double operator()(amichel::de::DVectorPtr args) {
	  auto coefs = toCoefs(args);

		vector<pair<double, double>> points2;
		for (auto& i : points) {
			points2.push_back({ i.first, calcPolynom(coefs, i.first) });
		}

		return distance(points, points2);
  }
};

class my_listener : public listener {
public:
	virtual void start() {
		//cout << setw(10) << "Generation" << setw(10) << "Metric" << endl;
	}
	virtual void end() {
		//cout << endl;
	}
	virtual void error() {}
	virtual void startGeneration(size_t genCount) {
	}
	virtual void endGeneration(size_t genCount, individual_ptr bestIndGen,
		individual_ptr bestInd) {
		if (genCount % 10 == 0) {
			//cout << "\r" << setw(10) << genCount << setw(10) << bestInd->cost();
		}
	}
	virtual void startSelection(size_t genCount) {}
	virtual void endSelection(size_t genCount) {}
	virtual void startProcessors(size_t genCount) {}
	virtual void endProcessors(size_t genCount) {}
};

vector<double> evolve(const vector<pair<double, double>>& points, int coefsCount, int generations, int populationSize) {
	try {
		constraints_ptr constraints(std::make_shared<constraints>(coefsCount, -1.0e3, 1.0e3));
		for (int i = 0; i < coefsCount; i++) {
			(*constraints)[i] = std::make_shared<real_constraint>(-100, 100);
		}

		sphere_function of(points);

		listener_ptr listener(std::make_shared<my_listener>());
		processor_listener_ptr processor_listener(std::make_shared<null_processor_listener>());

		processors<sphere_function>::processors_ptr _processors(std::make_shared<processors<sphere_function>>(32, std::ref(of), processor_listener));

		termination_strategy_ptr terminationStrategy(std::make_shared<max_gen_termination_strategy>(generations));

		selection_strategy_ptr selectionStrategy(std::make_shared<best_parent_child_selection_strategy>());

		mutation_strategy_arguments mutation_arguments(0.5, 0.9);
		mutation_strategy_ptr mutationStrategy(std::make_shared<mutation_strategy_1>(coefsCount, mutation_arguments));
		differential_evolution<sphere_function> de(
			coefsCount, 
			populationSize, 
			_processors, 
			constraints, 
			false /* is minimize? */,
			terminationStrategy, 
			selectionStrategy, 
			mutationStrategy, 
			listener
		);

		de.run();

		individual_ptr best(de.best());
		//std::cout << "Best metric after evolution: " << of(best->vars()) << endl;
		return toCoefs(best->vars());
	} catch (const amichel::de::exception& e) {
		std::cout << "an error occurred: " << e.what();
		return {};
	}
}

void addNoise(double& value, double noisePercent) {
	value *= 1 + random(-noisePercent, noisePercent);
	//value *= 1 + normalRandom()*noisePercent;
}

double calcEvolutionDistance(vector<double>& coefs, int generations, int populationSize, int pointsCount, double noisePercent, bool isPrint) {
	vector<pair<double, double>> points;
	for (int i = 0; i < pointsCount; i++) {
		double x = random(-20, 20);
		double y = calcPolynom(coefs, x);
		addNoise(y, noisePercent);
		points.push_back({ x, y });
	}

	auto result = evolve(points, coefs.size(), generations, populationSize);

	return distance(coefs, result);
}

int main() {
	vector<double> coefs = { -10, 0, 3.1, 88, 5, -15, 10 };

	auto result = calcEvolutionDistance(coefs, 200, 200, 20, 0.0, true);
	cout << result << endl;

	std::ofstream fout("out.txt");

	int counter = 0;
	for (double percent = 0; percent < 0.10; percent += 0.01) {
		for (int points = 8; points < 1500; points *= 1.2) {
			counter++;
			std::cout << "\r" << std::setw(5) << counter * 100 / (28 * 10) << "%";
			auto result = calcEvolutionDistance(coefs, 200, 200, points, percent, false);
			fout << percent << "\t" << points << "\t" << result << std::endl;
		}
	}
	fout.close();
}
