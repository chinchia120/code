#include <iostream>
using namespace std;

int main() {
    int a;
    while (cin >> a) {
        a = (a * a * a + 5 * a + 6) / 6;
        cout << a << endl;
    }
    return 0;
}