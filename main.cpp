#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    parallelor d1=parallelor();
    dijkstor d2=dijkstor();
    ERGraph graph(1000,1,d1,d2);
    graph.prepush(0,9,1);
}

