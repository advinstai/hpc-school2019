#include <iostream>
#include <cstdlib>

static const int val1 = -5;
const int val2 = 10;
static int val3 = -20;
int val4 = -15;
extern int errno;

class Test {
public:
int add_abs(const int v1, const int v2)
{
	return abs(v1)+abs(v2);
};
private:
int add_abs(const double v1, double v2)
{
	return abs(v1)+abs(v2);
};
};

namespace func {
static int add_abs(const int v1, const int v2)
{
	return abs(v1)+abs(v2);
}
}

static int add_abs(const int v1, int v2)
{
	return abs(v1)+abs(v2);
}

extern "C" int add_abs(double v1, double v2)
{
	return abs(v1)+abs(v2);
}


int main(int argc, char **argv)
{
     class Test t;
     int val5 = 20;

     std::cout << add_abs(val1,val2) << " / "
               << add_abs(val3,val4) << " / "
               << t.add_abs(errno,val5) << std::endl;
     return 0;
}
