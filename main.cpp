#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    struct timeval starttime, endtime;
    parallelor d=parallelor();
    //dijkstor d=dijkstor();
    ERGraph graph(800,1,d);
    gettimeofday(&starttime,NULL);
    cout<<"prepushing "<<endl;
    graph.prepush(33,88,1);
    gettimeofday(&endtime,NULL);
    cout<<"time is:"<<endtime.tv_sec- starttime.tv_sec<<endl;
}

//ERGraph graph(100,1,d);
// gettimeofday(&starttime,NULL);
//graph.prepush(66,33,1);
