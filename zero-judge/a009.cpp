#include <iostream>
using namespace std;

int main() {
    char ch;
    while (cin >> ch) {
        int k = int(ch) - 7;
        cout << char(k);
    }
    return 0;
}
