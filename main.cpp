#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <vector>
#include <chrono>
#include "evolvent.h"
#include "grishagin_function.h"
#include <omp.h>

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
	std::multimap<double, std::set<Point>::iterator> Rmap;
};
int n = 0;

struct Result {
	std::vector<double> y;
	double f;
	Result(size_t N) : f(0.0) { y.resize(N); }
};

double Calculate(vagrish::GrishaginFunction& foo, const double* y) {
	volatile double a = 1;
	for (size_t i = 0; i < NUM_ITER; i++)
		a *= sin(y[0]) * sin(y[0]) + cos(y[0]) * cos(y[0]);
	return foo.Calculate(y) * a;
}

Result GSA(vagrish::GrishaginFunction& foo, double* a, double* b, double r = 2.0, double epsilon = 0.001, int N = 2, int NMax = 1000)
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
	Spec spec;
	double M = abs((t.z() - t_1.z()) / (B - A));
	spec.Mset.insert(M);
	double m = (M > 0) ? r * M : 1;
	double R = (m * (B - A)) + ((t.z() - t_1.z()) * (t.z() - t_1.z()) / (m * (B - A))) - 2 * (t.z() + t_1.z());
	spec.Rmap.insert(std::make_pair(R, ++w.begin()));
	while (((t.x() - t_1.x()) > epsilon) && (n < NMax))
	{
		n++;

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

double Paraboloid(std::vector<double> x) {
	return 3.0 * x[0] * x[0] + 1.6 * x[1] * x[1];
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
		Result res = GSA(func, a.data(), b.data(), 11.5, 0.01, 2, 10000);
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
	}
	return 0;
}