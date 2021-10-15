#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <vector>

#define M_PI           3.14159265358979323846 

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

Point GSA(std::function<double(double)> foo, double a, double b, double r = 2, double epsilon = 0.00001, int N = 10000)
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
	double R = m * (b - a) + (foo(b) - foo(a)) * (foo(b) - foo(a)) / (m * (b - a)) - 2 * (foo(b) - foo(a));
	spec.Rmap.insert(std::make_pair(R, ++w.begin()));
	int n = 0;
	while (((t.x() - t_1.x()) > epsilon) && (n < N))
	{
		n++;

		// Проведение испытания 
		double x = 0.5 * (t.x() + t_1.x()) - (t.z() - t_1.z()) / (2 * m);
		double z = foo(x);
		res = z < res.z() ? Point(x, z) : res;

		// Вставка в w новой точки
		auto insert_return = w.insert(Point(x, z));

		std::set<Point>::iterator new_elem = insert_return.first;
		std::set<Point>::iterator before_new_elem = --new_elem;
		new_elem++;
		std::set<Point>::iterator after_new_elem = ++new_elem;
		new_elem--;

		// Вычисление M
		spec.Mset.erase(abs((after_new_elem->z() - before_new_elem->z()) / (after_new_elem->x() - before_new_elem->x())));
		spec.Mset.insert(abs((after_new_elem->z() - new_elem->z()) / (after_new_elem->x() - new_elem->x())));
		spec.Mset.insert(abs((new_elem->z() - before_new_elem->z()) / (new_elem->x() - before_new_elem->x())));
		if (M == *(--spec.Mset.end()))
		{
			// Вычисление R
			spec.Rmap.erase(m * (after_new_elem->x() - before_new_elem->x()) + (after_new_elem->z() - before_new_elem->z())
				* (after_new_elem->z() - before_new_elem->z()) / (m * (after_new_elem->x() - before_new_elem->x()))
				- 2 * (after_new_elem->z() - before_new_elem->z()));

			spec.Rmap.insert(std::make_pair(m * (after_new_elem->x() - new_elem->x()) + (after_new_elem->z() - new_elem->z())
				* (after_new_elem->z() - new_elem->z()) / (m * (after_new_elem->x() - new_elem->x()))
				- 2 * (after_new_elem->z() - new_elem->z()), after_new_elem));

			spec.Rmap.insert(std::make_pair(m * (new_elem->x() - before_new_elem->x()) + (new_elem->z() - before_new_elem->z())
				* (new_elem->z() - before_new_elem->z()) / (m * (new_elem->x() - before_new_elem->x()))
				- 2 * (new_elem->z() - before_new_elem->z()), new_elem));
		}
		else
		{
			// пересчет R
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
					* (i->z() - i_1->z()) / (m * (i->x() - i_1->x())) - 2 * (i->z() - i_1->z()), i));
			}
		}

		// Нахождение интервала с максимальным R 
		std::set<Point>::iterator ti = (--spec.Rmap.end())->second; // если R одинаковые?
		t.setx(ti->x());
		t.setz(ti->z());
		ti--;
		t_1.setx(ti->x());
		t_1.setz(ti->z());
	}
	return res;
}

double foo1(double x)
{
	return sin(x) + sin(10 * x / 3);
}

double foo2(double x)
{
	return 2 * (x - 3) * (x - 3) + exp(x * x / 2);
}

double foo3(double x)
{
	double res = 0;
	for (int k = 1; k <= 5; k++)
	{
		res += k * sin((k + 1) * x + k);
	}
	return (-1) * res;
}

double foo4(double x)
{
	return (3 * x - 1.4) * sin(18 * x);
}

double foo5(double x)
{
	return (-1) * (x + sin(x)) * exp((-1) * x * x);
}

double foo6(double x)
{
	return sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3;
}

double foo7(double x)
{
	return (-1) * sin(2 * M_PI * x) * exp((-1) * x);
}

double foo8(double x)
{
	return (x * x - 5 * x + 6) / (x * x + 1);
}

double foo9(double x)
{
	return (-1) * x + sin(3 * x) - 1;
}

int main(int argc, char* argv[])
{
	std::vector<std::function<double(double)>> functions({ foo1, foo2, foo3, foo4, foo5, foo6, foo7, foo8, foo9 });
	Point point;
	std::function<double(double)> foo = foo1;
	double a = -3;
	double b = 3;
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
				}
				if (argv[i][1] == 'a')
				{
					int j = 1;
					while (argv[i][j] != '=')
						j++;
					a = atof(&argv[i][++j]);
				}
				if (argv[i][1] == 'b')
				{
					int j = 1;
					while (argv[i][j] != '=')
						j++;
					b = atof(&argv[i][++j]);
				}
			}
		}
	}
	point = GSA(foo2, a, b, 2);
	std::cout << "x = " << point.x() << std::endl << "z = " << point.z();
	return 0;
}

