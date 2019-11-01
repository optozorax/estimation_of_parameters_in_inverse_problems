#include <random>
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <functional>
#include <fstream>

class PolynomCreature;
typedef std::shared_ptr<PolynomCreature> PolynomCreature_ptr;

class RandomGenerator {
public:
	RandomGenerator(int seed) : generator(seed) {
	}

	double randomDouble(double a, double b) {
		std::uniform_real_distribution<double> distribution(a, b);
		return distribution(generator);
	}

	int randomInt(int a, int b) {
		std::uniform_int_distribution<int> distribution(a, b);
		return distribution(generator);
	}
private:
	std::mt19937 generator;
};

RandomGenerator rnd(0);

class PolynomCreature {
public:
	void mutateSingle(void) {
		values[rnd.randomInt(0, values.size()-1)] *= rnd.randomDouble(-2, 2);
	}

	void mutate(int count) {
		for (int i = 0; i < count; ++i) {
			mutateSingle();
		}
	}

	PolynomCreature() = default;

	// Конструктор копирования, нужен для создания неизменных особей
	PolynomCreature(const PolynomCreature& other) = default;

	double calcPolynomInPoint(double x) const {
		double sum = 0;
		double powx = 1;
		for (auto& i : values) {
			sum += i * powx;
			powx *= x;
		}
		return sum;
	}

	std::vector<double> values;
};

PolynomCreature_ptr makeRandom(int size, double a, double b) {
	PolynomCreature result;
	result.values.resize(size);
	for (auto& i : result.values) {
		i = rnd.randomDouble(a, b);
	}
	return PolynomCreature_ptr(new PolynomCreature(result));
}

double evalByPosition(double value, double sum, int pos, int size) {
	return 1 / pow(1.1, pos);
}

double evalByValue(double value, double sum, int pos, int size) {
	return value;
}

template<class Creature>
class Evolution {
public:
	typedef std::shared_ptr<Creature> Creature_ptr;
	typedef std::vector<Creature_ptr> Creatures;
	typedef std::pair<Creature_ptr, double> CreatureWithScore;
	typedef std::vector<CreatureWithScore> CreatureScores;
	typedef std::function<Creature_ptr(void)> CreatureGenerator;
	typedef std::function<double(const Creature&)> EvalFunction;

	// Maximizes 
	Evolution(
		int populationSize, 
		EvalFunction max, 
		CreatureGenerator generator
	) : max(max), generator(generator) {
		init(populationSize);
	}

	void init(int populationSize) {
		population.clear();
		population.reserve(populationSize);
		for (int i = 0; i < populationSize; ++i) {
			population.push_back(generator());
		}
	}

	void evolve(int generations, bool isPrint) {
		for (int i = 0; i < generations; ++i) {
			evolveIter();
			if (i % 10 == 0 && isPrint) {
				std::cout << "\r" << std::setw(10) << std::fixed << std::setprecision(1) << i * 100.0 / generations << '%';
			}
		}
	}

	CreatureWithScore getBest(void) const {
		return evalScores(evalByPosition)[0];
	}

private:
	Creatures population;
	EvalFunction max;
	CreatureGenerator generator;

	CreatureScores evalScores(const std::function<double(double value, double sum, int pos, int size)>& normalize) const {
		CreatureScores result;
		double sum = 0;

		// Оцениваем каждое существо
		for (const auto& i : population) {
			double score = max(*i);
			result.push_back({i, score});
			sum += score;
		}

		// Сортируем существ по их оценкам
		// TODO сравнить насколько быстро программа будет работать без этой сортировки
		std::sort(result.begin(), result.end(), [] (const auto& a, const auto& b) -> bool {
			return a.second > b.second;
		});

		// Применяем к числам определённую нормализацию
		int pos = 0;
		sum = 0;
		for (auto& i : result) {
			i.second = normalize(i.second, sum, pos, result.size());
			pos++;
		}

		// Нормализуем все оценки таким образом, чтобы их сумма была равна 1
		for (auto& i : result)
			i.second /= sum;

		return result;
	}

	Creatures makeNewPopulation(const CreatureScores& scores) {
		Creatures result;
		result.reserve(scores.size() + 3);

		// Получаем существо с некоторой вероятностью (у лучших больше вероятность)
		auto getProbablyBestCreature = [&scores] (void) -> const Creature_ptr {
			double p = rnd.randomDouble(0, 1);
			for (int i = 0; i < scores.size(); ++i) {
				if (p < scores[i].second) {
					return scores[i].first;
				}
				p -= scores[i].second;
			}
			throw std::exception();
		};

		int size = scores.size();

		// Выбираем 5% без изменений
		for (int i = 0; i < 0.05 * size; ++i) {
			Creature_ptr creature(new Creature(*getProbablyBestCreature()));
			result.emplace_back(creature);
		}

		// Мутируем 45% одной мутацией
		for (int i = 0; i < 0.45 * size; ++i) {
			result.emplace_back(new Creature(*getProbablyBestCreature()));
			result.back()->mutate(1);
		}

		// Мутируем 10% двумя мутациями
		for (int i = 0; i < 0.10 * size; ++i) {
			result.emplace_back(new Creature(*getProbablyBestCreature()));
			result.back()->mutate(2);
		}

		// Мутируем 10% тремя мутациями
		for (int i = 0; i < 0.10 * size; ++i) {
			result.emplace_back(new Creature(*getProbablyBestCreature()));
			result.back()->mutate(3);
		}

		// Мутируем 10% четырьмя мутациями
		for (int i = 0; i < 0.10 * size; ++i) {
			result.emplace_back(new Creature(*getProbablyBestCreature()));
			result.back()->mutate(4);
		}

		// Мутируем 10% пятью мутациями
		for (int i = 0; i < 0.10 * size; ++i) {
			result.emplace_back(new Creature(*getProbablyBestCreature()));
			result.back()->mutate(5);
		}

		// Создаём новые 10%
		for (int i = 0; i < 0.10 * size; ++i) {
			result.push_back(generator());
		}

		// Нормализуем количество особей
		while (result.size() > size) {
			result.pop_back();
		}
		while (result.size() < size) {
			result.push_back(generator());
		}

		return result;
	}

	void evolveIter(void) {
		auto scores = evalScores(evalByPosition);
		population = makeNewPopulation(scores);
	}
};


typedef std::vector<std::pair<double, double>> Points;

double distanceAbsolute(double a, double b) {
	return pow(a - b, 2);
}

double distance(const PolynomCreature& coefs1, const PolynomCreature& coefs2) {
	double sum = 0;
	for (int i = 0; i < coefs1.values.size(); i++) {
		sum += distanceAbsolute(coefs1.values[i], coefs2.values[i]);
	}
	return sqrt(sum/coefs1.values.size());
}

double calcUniformMetric(double value, double x, double percent) {
	double r = fabs(value * percent * 2);
	if (fabs(x - value) < r) {
		return 1 + (1 - fabs(x - value)/r) * 9;
	} else {
		return r / fabs(x - value);
	}
}

double distance(const Points& a, const PolynomCreature& c) {
	double percent = 0.2;
	double sum = 0;
	for (int i = 0; i < a.size(); i++) {
		sum += calcUniformMetric(a[i].second, c.calcPolynomInPoint(a[i].first), percent);
	}
	return sum / a.size();
}

double calcEvolutionDistance(PolynomCreature& creature, int generations, int populationSize, int pointsCount, double noisePercent, bool isPrint) {
	Points points;
	for (int i = 0; i < pointsCount; ++i) {
		double x = rnd.randomDouble(-20, 20);
		double y = creature.calcPolynomInPoint(x);
		y *= 1 + rnd.randomDouble(-noisePercent, noisePercent);
		points.push_back({x, y});
	}

	Evolution<PolynomCreature> e(populationSize, [&](const PolynomCreature& c) -> double {
		return distance(points, c);
	}, [&]() -> PolynomCreature_ptr {
		return makeRandom(creature.values.size(), -100, 100);
	});

	e.evolve(generations, isPrint);

	auto res = e.getBest();

	using namespace std;

	if (isPrint) {
		cout << "Metric for true coefficients: " << distance(points, creature) << endl;
		cout << endl;

		cout << setprecision(3);
		cout << setw(10) << "Found" << setw(10) << "True" << setw(10) << "|f - t|" << endl;
		cout << setw(10) << "----------" << setw(10) << "----------" << setw(10) << "----------" << endl;
		auto& coefs = creature.values;
		auto& result = res.first->values;
		for (int i = 0; i < coefs.size(); i++) {
			cout << setw(10) << result[i] << setw(10) << coefs[i] << setw(10) << fabs(result[i] - coefs[i]) << endl;
		}
		cout << endl;

		cout << "Metric for result: " << distance(points, *res.first) << endl;
		cout << "Distance: " << distance(creature, *res.first) << endl;
	}

	return distance(creature, *res.first);
}

int main() {
	PolynomCreature creature;
	creature.values = {-10, 0, 3.1, 88, 5, -15, 10};
	auto result = calcEvolutionDistance(creature, 2000, 200, 40, 0.0, true);

	std::ofstream fout("out.txt");

	int counter = 0;
	for (double percent = 0; percent < 0.10; percent += 0.01) {
		for (int points = 8; points < 1500; points *= 1.2) {
			counter++;
			std::cout << "\r" << std::setw(5) << counter * 100 / (28 * 10) << "%";
			auto result = calcEvolutionDistance(creature, 2000, 200, points, percent, false);
			fout << percent << "\t" << points << "\t" << result << std::endl;
		}
	}
	fout.close();
}