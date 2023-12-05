#include <stdio.h>

struct time
{
    int hours;
    int minutes;
};

void showTime(struct time *);

int main()
{
    struct time now, *ptr;

    ptr = &now;

    ptr->hours = 18;
    // now.hours = 18;
    //(*ptr).hours = 18;

    (*ptr).minutes = 35;
    // now.minutes = 35;
    // ptr->minutes = 35;

    printf("   %2d : %2d\n", now.hours, now.minutes);
    showTime(ptr);

    return 0;
}

void showTime(struct time *ptr)
{
    if (ptr->hours >= 12)
    {
        printf("PM %2d : %2d\n", ptr->hours - 12, ptr->minutes);
    }
    else
    {
        printf("PM %2d : %2d\n", ptr->hours, ptr->minutes);
    }
}