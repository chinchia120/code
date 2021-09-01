#include <iostream>
using namespace std;

int main() {
    int a;
    while (cin >> a) {
        a = a * a - a + 2;
        cout << a << endl;
    }
    return 0;
}
