#include <iostream>
using namespace std;
int main(void) {
    long long int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        long long int start, end, sum = 0;
        scanf("%d %d", &start, &end);
        for (long long int j = start; j <= end; j++) {
            sum = sum + j;
            if (sum < end) {
                continue;
            } else {
                printf("%d\n", j);
                break;
            }
        }
    }
    return 0;
}
