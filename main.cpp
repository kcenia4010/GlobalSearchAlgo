#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <vector>
#include <thread>
#include <omp.h>

#define M_PI           3.14159265358979323846 
#define NUM_ITER 10000000

class Point
{
private:
	double x_;
	double z_;
public:
	Point() : x_(0.0), z_(0.0) {}
	Point(double _x, double _z) : x_(_x), z_(_z) {}
	double x() const { return x_; }
	double z() const { return z_; }
	bool operator < (const Point& other) const
	{
		return (this->x_ < other.x_);
	}
	void setx(double _x) { x_ = _x; }
	void setz(double _z) { z_ = _z; }
};

struct Spec
{
	std::multiset<double> Mset;
	std::multimap<double, std::pair<std::set<Point>::iterator, std::set<Point>::iterator>> Rmap;
};

void recalculationR(Spec& spec, std::set<Point>& w, double r, double& M, double& m) {
	spec.Rmap.clear();
	M = *(--spec.Mset.end());
	m = (M > 0) ? r * M : 1;
	std::set<Point>::iterator i = w.begin();
	i++;
	for (i; i != w.end(); i++)
	{
		std::set<Point>::iterator i_1 = --i;
		i++;
		spec.Rmap.insert(std::make_pair(m * (i->x() - i_1->x()) + (i->z() - i_1->z())
			* (i->z() - i_1->z()) / (m * (i->x() - i_1->x())) - 2 * (i->z() + i_1->z()), std::make_pair(i_1, i)));
	}
}

void new_interval(Point& t, Point& t_1, Spec& spec)
{
	std::set<Point>::iterator ti = ((--spec.Rmap.end())->second).second;
	t.setx(ti->x());
	t.setz(ti->z());
	ti--;
	t_1.setx(ti->x());
	t_1.setz(ti->z());
}

void insert_point(Spec& spec, std::set<Point>& w, Point point, double M, double m)
{
	auto insert_return = w.insert(point);

	std::set<Point>::iterator new_elem = insert_return.first;
	std::set<Point>::iterator before_new_elem = --new_elem;
	new_elem++;
	std::set<Point>::iterator after_new_elem = ++new_elem;
	new_elem--;

	spec.Mset.erase(abs((after_new_elem->z() - before_new_elem->z()) / (after_new_elem->x() - before_new_elem->x())));
	spec.Mset.insert(abs((after_new_elem->z() - new_elem->z()) / (after_new_elem->x() - new_elem->x())));
	spec.Mset.insert(abs((new_elem->z() - before_new_elem->z()) / (new_elem->x() - before_new_elem->x())));

	if (M == *(--spec.Mset.end()))
	{
		spec.Rmap.erase(m * (after_new_elem->x() - before_new_elem->x()) + (after_new_elem->z() - before_new_elem->z())
			* (after_new_elem->z() - before_new_elem->z()) / (m * (after_new_elem->x() - before_new_elem->x()))
			- 2 * (after_new_elem->z() + before_new_elem->z()));

		spec.Rmap.insert(std::make_pair(m * (after_new_elem->x() - new_elem->x()) + (after_new_elem->z() - new_elem->z())
			* (after_new_elem->z() - new_elem->z()) / (m * (after_new_elem->x() - new_elem->x()))
			- 2 * (after_new_elem->z() + new_elem->z()), std::make_pair(new_elem, after_new_elem)));

		spec.Rmap.insert(std::make_pair(m * (new_elem->x() - before_new_elem->x()) + (new_elem->z() - before_new_elem->z())
			* (new_elem->z() - before_new_elem->z()) / (m * (new_elem->x() - before_new_elem->x()))
			- 2 * (new_elem->z() + before_new_elem->z()), std::make_pair(before_new_elem, new_elem)));
	}
}

void insert_point_parallel(Spec& spec, std::set<Point>& w, Point point, double M, double m,
	double newM1, double newM2, double Mprev, double Rprev = 0, double newR1 = 0, double newR2 = 0)
{
	auto insert_return = w.insert(point);

	std::set<Point>::iterator new_elem = insert_return.first;
	std::set<Point>::iterator before_new_elem = --new_elem;
	new_elem++;
	std::set<Point>::iterator after_new_elem = ++new_elem;
	new_elem--;

	spec.Mset.erase(Mprev);
	spec.Mset.insert(newM2);
	spec.Mset.insert(newM1);

	if (M == *(--spec.Mset.end()))
	{
		spec.Rmap.erase(m * (after_new_elem->x() - before_new_elem->x()) + (after_new_elem->z() - before_new_elem->z())
			* (after_new_elem->z() - before_new_elem->z()) / (m * (after_new_elem->x() - before_new_elem->x()))
			- 2 * (after_new_elem->z() + before_new_elem->z()));

		spec.Rmap.insert(std::make_pair(newR2, std::make_pair(new_elem, after_new_elem)));
		spec.Rmap.insert(std::make_pair(newR1, std::make_pair(before_new_elem, new_elem)));
	}
}


void search(std::pair<double, std::pair<std::set<Point>::iterator, std::set<Point>::iterator>> R,
	std::function<double(double)> foo, double m, double M,
	double& prevM, double& prevR, double& newM1, double& newM2, double& newR1, double& newR2,
	Point& new_point, int& return_code)
{
	Point t = *R.second.second;
	Point t_1 = *R.second.first;
	double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
	double z = foo(x);
	new_point = Point(x, z);
	prevM = abs((t.z() - t_1.z()) / (t.x() - t_1.x()));
	if (prevM != M)
	{
		newM1 = abs((z - t_1.z()) / (x - t_1.x()));
		newM2 = abs((t.z() - z) / (t.x() - x));
		if (M >= std::max(newM1, newM2))
		{
			prevR = m * (t.x() - t_1.x()) + (t.z() - t_1.z()) * (t.z() - t_1.z()) / (m * (t.x() - t_1.x())) - 2 * (t.z() + t_1.z());
			newR1 = m * (x - t_1.x()) + (z - t_1.z()) * (z - t_1.z()) / (m * (x - t_1.x())) - 2 * (z + t_1.z());
			newR2 = m * (t.x() - x) + (t.z() - z) * (t.z() - z) / (m * (t.x() - x)) - 2 * (t.z() + z);
			return_code = 0;
		}
		else
			return_code = 2;
	}
	else
		return_code = 1;
}
int n = 0;
Point GSA(std::function<double(double)> foo, double a, double b, double r = 2, double epsilon = 0.0001, int N = 10000)
{
	std::set<Point> w;
	w.insert(Point(a, foo(a)));
	w.insert(Point(b, foo(b)));
	Point t(b, foo(b));
	Point t_1(a, foo(a));
	Point res = (foo(a) < foo(b)) ? t_1 : t;
	Spec spec;
	double M = abs((foo(b) - foo(a)) / (b - a));
	spec.Mset.insert(M);
	double m = (M > 0) ? r * M : 1;
	double R = m * (b - a) + (foo(b) - foo(a)) * (foo(b) - foo(a)) / (m * (b - a)) - 2 * (foo(b) + foo(a));
	spec.Rmap.insert(std::make_pair(R, std::make_pair(w.begin(), ++w.begin())));
	int n_thread = 4;
	double* prevM = new double[n_thread];
	double* prevR = new double[n_thread];
	double* newM1 = new double[n_thread];
	double* newM2 = new double[n_thread];
	double* newR1 = new double[n_thread];
	double* newR2 = new double[n_thread];
	Point* new_point = new Point[n_thread];
	int* return_code = new int[n_thread];
	std::vector<std::thread> threads(n_thread);
	int count = 0;
	while (((t.x() - t_1.x()) > epsilon) && (n < N))
	{
		if (n < 4)
		{
			double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
			double z = foo(x);
			res = z < res.z() ? Point(x, z) : res;
			insert_point(spec, w, Point(x, z), M, m);

			if (M != *(--spec.Mset.end()))
			{
				recalculationR(spec, w, r, M, m);
			}

			new_interval(t, t_1, spec);
		}
		else {
			auto iter = spec.Rmap.end();
			std::vector<std::pair<double, std::pair<std::set<Point>::iterator, std::set<Point>::iterator>>> Rmax(n_thread);
			for (int i = 0; i < n_thread; i++)
			{
				iter--;
				Rmax[i] = *iter;
			}
			int j = 0;
			for (auto it = std::begin(threads); it != std::end(threads); ++it)
			{
				*it = std::thread(search, Rmax[j], foo, m, M, std::ref(prevM[j]), std::ref(prevR[j]), std::ref(newM1[j]), std::ref(newM2[j]),
					std::ref(newR1[j]), std::ref(newR2[j]), std::ref(new_point[j]), std::ref(return_code[j]));
				j++;
			}
			for (auto&& i : threads) {
				i.join();
			}
			int idx = -1;
			double maxR_elem = (*(--spec.Rmap.end())).first;
			bool f = false;
			for (int i = 0; i < n_thread; i++)
			{
				if ((return_code[i] == 2) || (return_code[i] == 1))
				{
					f = true;
					insert_point_parallel(spec, w, new_point[i], M, m,
						newM1[i], newM2[i], prevM[i]);
					recalculationR(spec, w, r, M, m);
					res = new_point[i].z() < res.z() ? new_point[i] : res;
					break;
				}
			}
			if (f == false)
			{
				for (int i = 0; i < n_thread; i++)
				{
					insert_point_parallel(spec, w, new_point[i], M, m,
						newM1[i], newM2[i], prevM[i], prevR[i], newR1[i], newR2[i]);
					res = new_point[i].z() < res.z() ? new_point[i] : res;
				}
			}
			new_interval(t, t_1, spec);
		}
		n++;
	}
	delete[] prevM;
	delete[] prevR;
	delete[] newM1;
	delete[] newM2;
	delete[] newR1;
	delete[] newR2;
	delete[] new_point;
	delete[] return_code;
	return res;
}

double foo0(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return sin(x) + sin(10 * x / 3) * a;
}

double foo1(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return 2 * (x - 3) * (x - 3) + exp(x * x / 2) * a;
}

double foo2(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	double res = 0;
	for (int k = 1; k <= 5; k++)
	{
		res += k * sin((k + 1) * x + k);
	}
	return (-1) * res * a;
}

double foo3(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return (3 * x - 1.4) * sin(18 * x) * a;
}

double foo4(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return (-1) * (x + sin(x)) * exp((-1) * x * x) * a;
}

double foo5(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3 * a;
}

double foo6(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return (-1) * sin(2 * M_PI * x) * exp((-1) * x) * a;
}

double foo7(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return (x * x - 5 * x + 6) / (x * x + 1) * a;
}

double foo8(double x)
{
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(x) * sin(x) + cos(x) * cos(x);
	return (-1) * x + sin(3 * x) - 1 * a;
}

std::vector<std::vector<double>> bounds = { {2.7, 7.5}, {-3.0, 3.0},  {0.0, 10.0}, {0.0, 1.2},
	{-10.0, 10.0}, {2.7, 7.5}, {0.0, 4.0}, {-5.0, 5.0}, {0.0, 6.5} };

int main(int argc, char* argv[])
{
	std::vector<std::function<double(double)>> functions({ foo0, foo1, foo2, foo3, foo4, foo5, foo6, foo7, foo8 });
	Point point;
	int f = 0;
	std::function<double(double)> foo = functions[f];
	double a = bounds[f][0];
	double b = bounds[f][1];
	if (argc > 1)
	{
		for (int i = 1; i < argc; i++)
		{
			if (argv[i][0] == '-')
			{
				if (argv[i][1] == 'f')
				{
					int j = 1;
					while (argv[i][j] != '=')
						j++;
					int k = atoi(&argv[i][++j]) - 1;
					foo = functions[k];
					a = bounds[k][0];
					b = bounds[k][1];
				}
			}
		}
	}
	for (f = 0; f < 9; f++) {
		n = 0;
		std::cout << "function" << f << std::endl;
		double st = omp_get_wtime();
		point = GSA(functions[f], bounds[f][0], bounds[f][1], 2.5, 0.01, 10000);
		double en = omp_get_wtime();
		std::cout << "time omp = " << en - st << std::endl;
		std::cout << "x = " << point.x() << std::endl << "z = " << point.z() << std::endl;
		std::cout << "n = " << n << std::endl << std::endl;
	}
	return 0;
}