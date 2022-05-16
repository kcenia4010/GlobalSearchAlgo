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

#define _CRT_SECURE_NO_WARNINGS

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
	std::multimap<double, std::set<Point>::iterator> Rmap;
};

struct Result {
	std::vector<double> y;
	double f;
	Result(size_t N = 2) : f(0.0) { y.resize(N); }
};

void GSA(vagrish::GrishaginFunction& foo, double* a, double* b, double r, double epsilon, int N, int NMax, Result& point, int& num_iter)
{
	TEvolvent evolvent(N, 10);
	Result result(N);
	evolvent.SetBounds(a, b);
	std::set<Point> w;
	std::vector<double> y1(N);
	std::vector<double> y2(N);
	double A = a[0], B = b[0];
	evolvent.GetImage(A, y1.data());
	evolvent.GetImage(B, y2.data());
	w.insert(Point(A, foo.Calculate(y1.data())));
	w.insert(Point(B, foo.Calculate(y2.data())));
	Point t(B, foo.Calculate(y2.data()));
	Point t_1(A, foo.Calculate(y1.data()));
	result.f = (t_1.z() < t.z()) ? t_1.z() : t.z();
	Spec spec;
	double M = abs((t.z() - t_1.z()) / (B - A));
	spec.Mset.insert(M);
	double m = (M > 0) ? r * M : 1;
	double R = (m * (B - A)) + ((t.z() - t_1.z()) * (t.z() - t_1.z()) / (m * (B - A))) - 2 * (t.z() + t_1.z());
	spec.Rmap.insert(std::make_pair(R, ++w.begin()));
	int n = 0;
	while (((t.x() - t_1.x()) > epsilon) && (n < NMax))
	{
		n++;

		std::vector<double> newY(N);

		double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
		evolvent.GetImage(x, newY.data());
		double z = foo.Calculate(newY.data());
		if (z < result.f) {
			result.f = z;
			for (int i = 0; i < N; ++i) {
				result.y[i] = newY[i];
			}
		}

		auto insert_return = w.insert(Point(x, z));

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
				- 2 * (after_new_elem->z() + new_elem->z()), after_new_elem));

			spec.Rmap.insert(std::make_pair(m * (new_elem->x() - before_new_elem->x()) + (new_elem->z() - before_new_elem->z())
				* (new_elem->z() - before_new_elem->z()) / (m * (new_elem->x() - before_new_elem->x()))
				- 2 * (new_elem->z() + before_new_elem->z()), new_elem));
		}
		else
		{
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
					* (i->z() - i_1->z()) / (m * (i->x() - i_1->x())) - 2 * (i->z() + i_1->z()), i));
			}
		}

		std::set<Point>::iterator ti = (--spec.Rmap.end())->second;
		t.setx(ti->x());
		t.setz(ti->z());
		ti--;
		t_1.setx(ti->x());
		t_1.setz(ti->z());
	}
	point = result;
	num_iter = n;
}

double num_iter;
Result GlobalSearchAlgo(vagrish::GrishaginFunction& foo, double* a, double* b, double r = 2, double epsilon = 0.00001, int N = 2, int NMax = 10000)
{
	int n_thread = 4;
	std::vector<int> num_iters(n_thread);
	std::vector<std::thread> threads(n_thread);
	double distance = abs(b[0] - a[0]);
	double start = a[0], end = a[0] + distance / n_thread;
	std::vector<Result> points(n_thread);
	int i = 0;
	std::vector<double> s(N);
	std::vector<double> e(N);
	for (auto it = std::begin(threads); it != std::end(threads); ++it) {
		s[0] = s[1] = start;
		e[0] = e[1] = end;
		*it = std::thread(GSA, std::ref(foo), s.data(), e.data(), r, epsilon, N, NMax / n_thread, std::ref(points[i]), std::ref(num_iters[i]));
		start = end;
		if (i == n_thread - 2)
			end = b[0];
		else
			end += distance / n_thread;
		i++;
	}
	for (auto&& i : threads) {
		i.join();
	}
	Result res = points[0];
	num_iter = num_iters[0];
	for (int j = 1; j < n_thread; j++) {
		if (points[j].f < res.f)
			res = points[j];
		if (num_iters[j] > num_iter)
			num_iter = num_iters[j];
	}

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
	std::vector<double> a(2);
	std::vector<double> b(2);
	int N = 2;
	vagrish::GrishaginFunction func;
	for (int i = 1; i < 101; i++) {
		num_iter = 0;
		func.SetFunctionNumber(i);
		func.GetBounds(a.data(), b.data());
		double st = omp_get_wtime();
		Result res = GlobalSearchAlgo(func, a.data(), b.data(), 2.5, 0.01, N, 10000);
		double en = omp_get_wtime();
		std::cout << "func " << i << std::endl;
		std::cout << "time omp = " << en - st << std::endl;
		std::cout << " n = " << num_iter << std::endl;
		std::cout << "f = " << res.f << std::endl;
		std::cout << "y = ";
		for (int i = 0; i < N; i++) {
			std::cout << res.y[i] << " ";
		}
		std::cout << std::endl << std::endl;
		double* point = new double[2];
		func.GetOptimumPoint(point);
		if ((func.GetOptimumValue() - res.f > 0.1)) 
		{
			std::cout << "fail in " << i << " res.f = " << res.f << std::endl;
			std::cout << "res.y = " << res.y[0] << " " << res.y[1] << std::endl;
			std::cout << std::endl << std::endl;
		}
	}
	return 0;
}