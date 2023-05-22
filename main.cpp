#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <vector>
#include <thread>
#include "evolvent.h"
#include <omp.h>
#include "grishagin_function.h"

#define M_PI           3.14159265358979323846 
#define NUM_ITER 1000000

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

struct Result {
	std::vector<double> y;
	double f;
	Result(size_t N = 2) : f(0.0) { y.resize(N); }
};

double Calculate(vagrish::GrishaginFunction& foo, const double* y) {
	return foo.Calculate(y);
}

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
		spec.Rmap.erase(Rprev);

		spec.Rmap.insert(std::make_pair(newR2, std::make_pair(new_elem, after_new_elem)));
		spec.Rmap.insert(std::make_pair(newR1, std::make_pair(before_new_elem, new_elem)));
	}
}


void search(std::pair<double, std::pair<std::set<Point>::iterator, std::set<Point>::iterator>> R, TEvolvent evolvent,
	vagrish::GrishaginFunction& foo, double m, double M,
	double& prevM, double& prevR, double& newM1, double& newM2, double& newR1, double& newR2,
	Point& new_point, int& return_code)
{
	Point t = *R.second.second;
	Point t_1 = *R.second.first;
	std::vector<double> newY(2);
	double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
	if ((x < 0) || (x > 1)) {
		x = 0.5 * (t.x() + t_1.x());
	}
	evolvent.GetImage(x, newY.data());
	double z = Calculate(foo, newY.data());
	new_point = Point(x, z);
	prevM = abs((t.z() - t_1.z()) / (t.x() - t_1.x()));
	newM1 = abs((z - t_1.z()) / (x - t_1.x()));
	newM2 = abs((t.z() - z) / (t.x() - x));
	if (prevM != M)
	{
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

void min(TEvolvent& evolvent, vagrish::GrishaginFunction& foo, double x, Result& result) {
	std::vector<double> newY(2);
	evolvent.GetImage(x, newY.data());
	double z = foo.Calculate(newY.data());
	if (z < result.f) {
		result.f = z;
		for (int i = 0; i < 2; ++i) {
			result.y[i] = newY[i];
		}
	}
}

int n = 0;
Result GSA(vagrish::GrishaginFunction& foo, double* a, double* b, double r = 2, double epsilon = 0.0001, int N = 2, int NMax = 10000)
{
	TEvolvent evolvent(N, 10);
	Result result(N);
	evolvent.SetBounds(a, b);
	std::set<Point> w;
	std::vector<double> y1(N);
	std::vector<double> y2(N);
	double A = 0.0, B = 1.0;
	evolvent.GetImage(A, y1.data());
	evolvent.GetImage(B, y2.data());
	w.insert(Point(A, Calculate(foo, y1.data())));
	w.insert(Point(B, Calculate(foo, y2.data())));
	Point t(B, Calculate(foo, y2.data()));
	Point t_1(A, Calculate(foo, y1.data()));
	result.f = (t_1.z() < t.z()) ? t_1.z() : t.z();
	result.y = (t_1.z() < t.z()) ? y1 : y2;
	Spec spec;
	double M = abs((t.z() - t_1.z()) / (B - A));
	spec.Mset.insert(M);
	double m = (M > 0) ? r * M : 1;
	double R = (m * (B - A)) + ((t.z() - t_1.z()) * (t.z() - t_1.z()) / (m * (B - A))) - 2 * (t.z() + t_1.z());
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
	while (((t.x() - t_1.x()) > epsilon) && (n < NMax))
	{
		if (n < 4)
		{
			std::vector<double> newY(N);

			double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
			evolvent.GetImage(x, newY.data());
			double z = Calculate(foo, newY.data());
			if (z < result.f) {
				result.f = z;
				for (int i = 0; i < N; ++i) {
					result.y[i] = newY[i];
				}
			}

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
				*it = std::thread(search, Rmax[j], std::ref(evolvent), std::ref(foo), m, M, std::ref(prevM[j]), std::ref(prevR[j]), std::ref(newM1[j]), std::ref(newM2[j]),
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
					min(evolvent, foo, new_point[i].x(), result);
					break;
				}
			}
			if (f == false)
			{
				for (int i = 0; i < n_thread; i++)
				{
					insert_point_parallel(spec, w, new_point[i], M, m,
						newM1[i], newM2[i], prevM[i], prevR[i], newR1[i], newR2[i]);
					min(evolvent, foo, new_point[i].x(), result);
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
	return result;
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
	std::vector<double> a(2);
	std::vector<double> b(2);
	int N = 2;
	vagrish::GrishaginFunction func;
	for (int i = 1; i < 6; i++) {
		n = 0;
		func.SetFunctionNumber(i);
		func.GetBounds(a.data(), b.data());
		double st = omp_get_wtime();
		Result res = GSA(func, a.data(), b.data(), 2.5, 0.0001, N, 10000);
		double en = omp_get_wtime();
		std::cout << "func " << i << std::endl;
		std::cout << "time omp = " << en - st << std::endl;
		std::cout << " n = " << n << std::endl;
		std::cout << "f = " << res.f << std::endl;
		std::cout << "y = ";
		for (int i = 0; i < N; i++) {
			std::cout << res.y[i] << " ";
		}
		std::cout << std::endl << std::endl;
		double* point = new double[2];
		func.GetOptimumPoint(point);
		if ((abs(func.GetOptimumValue() - res.f) > 0.1))
		{
			std::cout << "fail in " << i << " Optimum= " << func.GetOptimumValue() << std::endl;
			std::cout << std::endl << std::endl;
		}
	}
	return 0;
}