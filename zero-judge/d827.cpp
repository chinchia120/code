#include <iostream>
using namespace std;

int main() {
    int a, b, c;
    while (cin >> a) {
        b = a / 12;
        c = a % 12;
        cout << 50 * b + 5 * c << endl;
    }
    return 0;
}
