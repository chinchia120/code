#include <iostream>
#include <queue>
using namespace std;

class MyStack 
{
public:
    queue<int> my_stack;

    MyStack() 
    {
        
    }
    
    void push(int x) 
    {   
        queue<int> tmp_stack = my_stack;
        int len = my_stack.size();

        queue<int> my_empty;
        swap(my_empty, my_stack);
        
        my_stack.push(x);

        for(int i = 0; i < len; i++)
        {   
            my_stack.push(tmp_stack.front());
            tmp_stack.pop();
        }
    }
    
    int pop() 
    {
        int tmp = my_stack.front();
        my_stack.pop();

        return tmp;
    }
    
    int top() 
    {
        return my_stack.front();
    }
    
    bool empty() 
    {
        return my_stack.empty();
    }
};

int main(int argc, char **argv)
{   
    MyStack myStack;
    myStack.push(1);
    myStack.push(2);
    myStack.top();
    myStack.pop();
    myStack.empty();

    return 0;
}