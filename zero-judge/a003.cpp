#include <iostream>
using namespace std;

int main() {
    int a, b, s;
    cin >> a >> b;
    s = (a * 2 + b) % 3;
    if (s == 0) {
        cout << "´¶³q" << endl;
    }
    if (s == 1) {
        cout << "¦N" << endl;
    }
    if (s == 2) {
        cout << "¤j¦N" << endl;
    }
  
    return 0;
}
