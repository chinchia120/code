#include <iostream>
#include <ctype.h>
using namespace std;

int main()
{
	// qustion 1
	int a = 3, b = 5;
	cout << "if A=3 and B=5,we know that A+B=" << a + b << endl;
	// question 2
	float c = 1.2, d = 3.0;
	cout << "if A=1.2 and B=3.0,we know that A*B=" << c * d << endl;
	// qustion 3
	int e = 4;
	cout << "if diameter R=4, the area of circle=" << e / 2 * e / 2 * 3.14 << endl;
	// qustion 4
	char CH = 'A';
	cout << "if charcter is 'A' the ASCII code number is " << (int)'A' << endl;
	// qustion 5
	int p = 3, q;
	q = p + 96;
	cout << "if the number is 3,it corrresponds the lowercase " << (char)q << endl;
	// qustion 6
	float g, h;
	float i;
	cout << "input two integer A and B." << endl;
	cin >> g >> h;
	i = g / h;
	cout << "A/B=" << i << endl;
	// qustion 7
	float j, k, m;
	int n;
	cout << "input two integer A and B." << endl;
	cin >> j >> k;
	m = j / k;
	n = j / k;
	if (m - n > 0.4)
	{
		cout << "after rounding A/B=" << n + 1 << endl;
	}
	else if (m - n <= 0.4)
	{
		cout << "after rounding A/B=" << n << endl;
	}
	// qustion 8
	int l;
	cout << "input a integer A." << endl;
	cin >> l;
	if (l > 0)
	{
		cout << "the absolute value of A=" << l << endl;
	}
	if (l < 0)
	{
		cout << "the absolute value of A=" << -l << endl;
	}
	else if (l == 0)
	{
		cout << "the absolute value of A=" << 0 << endl;
	}
	// qustion 9
	char ch;
	cout << "input a charcter." << endl;
	cin >> ch;
	if (isalpha(ch))
	{
		if (isupper(ch))
		{
			cout << "it is uppercase" << endl;
		}
		else if (islower(ch))
		{
			cout << "it is lowercase" << endl;
		}
	}
	else if (isdigit(ch))
	{
		cout << "it is number" << endl;
	}
	else
	{
		cout << "error" << endl;
	}
	// qustion 10
	int total, hr, min, sec;
	cout << "input a integer A." << endl;
	cin >> total;
	hr = total / 3600;
	min = total % 3600 / 60;
	sec = total % 3600 % 60;
	cout << hr << "hour " << min << "min " << sec << "sec " << endl;

	return 0;
}