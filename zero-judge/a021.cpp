#include <iostream>
using namespace std;

int main() {
    long long int a, b, d;
    char c;
    cin >> a >> c >> b;
    if (c == '+') d = a + b;
    if (c == '-') d = a - b;
    if (c == '*') d = a * b;
    if (c == '/') d = a / b;
    cout << d << endl;

    return 0;
}
