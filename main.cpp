#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    struct timeval starttime, endtime;
    //parallelor d=parallelor();
    dijkstor d=dijkstor();
    ERGraph graph(10000,1,d);
    cout<<"graph init success"<<endl;
    gettimeofday(&starttime,NULL);
    cout<<"prepushing "<<endl;
    graph.prepush(44,33,1);
    gettimeofday(&endtime,NULL);
    cout<<"time is:"<<endtime.tv_sec- starttime.tv_sec<<endl;
}
//ERGraph graph(100,1,d);
// gettimeofday(&starttime,NULL);
//graph.prepush(66,33,1);
