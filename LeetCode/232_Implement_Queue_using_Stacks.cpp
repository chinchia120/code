#include <iostream>
#include <queue>
using namespace std;

class MyQueue 
{
public:
    queue<int> my_stack;

    MyQueue() 
    {
        
    }
    
    void push(int x) 
    {
        my_stack.push(x);
    }
    
    int pop() 
    {
        int tmp = my_stack.front();
        my_stack.pop();

        return tmp;
    }
    
    int peek() 
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
    MyQueue myQueue;
    myQueue.push(1);
    myQueue.push(2);
    myQueue.peek();
    myQueue.pop();
    myQueue.empty();

    return 0;
}