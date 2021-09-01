#include <math.h>

#include <iostream>

using namespace std;

int main() {
    int a, b, c, d, e;
    cin >> a >> b >> c;
    d = (b * b) - (4 * a * c);
    e = sqrt(d);
    if (d > 0) {
        cout << "Two different roots x1=" << (-b + e) / (2 * a) << " , x2=" << (-b - e) / (2 * a) << endl;
    } else if (d == 0) {
        cout << "Two same roots x=" << (-b + e) / (2 * a) << endl;
    } else {
        cout << "No real root" << endl;
    }

    return 0;
}
