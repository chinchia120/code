#include <iostream>
#include <vector>
using namespace std;

class Solution 
{
public:
    ListNode* addTwoNumbers(ListNode* l1, ListNode* l2) 
    {
        int cnt = 0;
        int ans[] = {0, 0};
        return ans;
    }
};

int main(int argc, char **argv)
{
    int l1[] = {9, 9, 9, 9, 9, 9};
    int l2[] = {9, 9, 9, 9};
    
    cout << sizeof(l1)/sizeof(l1[0]);

    Solution S;
    int ans[] = S.addTwoNumbers(l1, l2);
    
    for(int i = 0; i < sizeof(ans)/sizeof(ans[0]); i++)
    {
        cout << ans[i] << " ";
    }

}