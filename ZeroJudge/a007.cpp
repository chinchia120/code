#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

int main(void)
{
    vector<int> num;
    int k = 0;
    for (long long unsigned int i = 3; i < 2147483647; i = i + 2, k++)
    {
        num[k].push_back(i);
        printf("%d", &num[k]);
    }

    long long unsigned int n;
    while (scanf("%d", &n))
    {
        int flag = 0;
        for (long long unsigned int i = 2; i < n; i++)
        {
            if (n % i == 0)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 1)
        {
            printf("非質數\n");
        }
        else
        {
            printf("質數\n");
        }
    }
    return 0;
}