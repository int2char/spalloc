#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    parallelor d1=parallelor();
    dijkstor d2=dijkstor();
    ERGraph graph(100,1,d2,d1);
    cout<<"graph finied"<<endl;
    graph.prepush(0,3,1);
}

