#include <iostream>
using namespace std;

int main() {
    string st;
    char ch;
    int first, sum = 0, num1, num2;
    while (cin >> st) {
        ch = (st.at(0));
        if (ch == 'A') {
            first = 10;
        } else if (ch == 'B') {
            first = 11;
        } else if (ch == 'C') {
            first = 12;
        } else if (ch == 'D') {
            first = 13;
        } else if (ch == 'E') {
            first = 14;
        } else if (ch == 'F') {
            first = 15;
        } else if (ch == 'G') {
            first = 16;
        } else if (ch == 'H') {
            first = 17;
        } else if (ch == 'I') {
            first = 34;
        } else if (ch == 'J') {
            first = 18;
        } else if (ch == 'K') {
            first = 19;
        } else if (ch == 'L') {
            first = 20;
        } else if (ch == 'M') {
            first = 21;
        } else if (ch == 'N') {
            first = 22;
        } else if (ch == 'O') {
            first = 35;
        } else if (ch == 'P') {
            first = 23;
        } else if (ch == 'Q') {
            first = 24;
        } else if (ch == 'R') {
            first = 25;
        } else if (ch == 'S') {
            first = 26;
        } else if (ch == 'T') {
            first = 27;
        } else if (ch == 'U') {
            first = 28;
        } else if (ch == 'V') {
            first = 29;
        } else if (ch == 'W') {
            first = 32;
        } else if (ch == 'X') {
            first = 30;
        } else if (ch == 'Y') {
            first = 31;
        } else if (ch == 'Z') {
            first = 33;
        }
        for (int i = 1; i < 9; i++) {
            sum = sum + (st.at(i) - '0') * (9 - i);
        }
        sum = sum + (st.at(9) - '0') + (first / 10) + (first % 10) * 9;
        if (sum % 10 == 0) {
            cout << "real" << endl;
        } else {
            cout << "fake" << endl;
        }
    }
    return 0;
}