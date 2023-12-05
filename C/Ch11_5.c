#include <stdio.h>
#include <string.h>

int main(){
    typedef int onHand;
    struct item{
        char name[30];
        float cost;
        onHand quantity;
    };
    typedef struct item inventory;

    inventory phone;

    strcpy(phone.name, "iPhone");
    phone.cost = 27500.00;
    phone.quantity = 100;

    printf("item   name   = %s\n", phone.name);
    printf("item   cost   = %.2f\n", phone.cost);
    printf("item quantity = %d\n", phone.quantity);

    return 0;
}