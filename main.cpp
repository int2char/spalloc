#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    //parallelor d=parallelor();
    dijkstor d=dijkstor();
    ERGraph graph(10000,1,d);
    cout<<"graph init success"<<endl;
    cout<<"prepushing "<<endl;
    graph.prepush(0,99,1);
}

