#include <iostream>
using namespace std;

int main(void) {
    int a;
    while (cin >> a) {
        if (a == 0) {
            return 0;
        } else {
            for (int i = 1; i < a; i++) {
                if (i % 7 != 0) {
                    cout << i << " ";
                }
            }
            cout << endl;
        }
    }
    return 0;
}
