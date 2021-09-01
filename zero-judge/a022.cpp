#include <iostream>
using namespace std;

int main() {
    string str;
    while (cin >> str) {
        for (int i = 0; i < str.length(); i++) {
            if (str.at(i) == str.at(str.length() - i - 1)) {
                cout << "yes" << endl;
                break;
            } else {
                cout << "no" << endl;
                break;
            }
        }
    }
    return 0;
}