#include <iostream>
using namespace std;

int main() {
    int a, b, c;
    cin >> a >> b;
    do {
        c = a % b;
        a = b;
        b = c;
    } while (c != 0);
    cout << a << endl;
    return 0;
}
