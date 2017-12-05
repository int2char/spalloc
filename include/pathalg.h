//
// Created by root on 17-5-9.
//
#ifndef CSPALLOC_PATHALG_H
#include"limits.h"
#define CSPALLOC_PATHALG_H
#define INFCOST INT_MAX/2
#include<bits/stdc++.h>
#include <unistd.h>
#include"edge.h"
#define ML 50
#define BS 5
#define WD 5
#define inf INT_MAX/2
using namespace std;
class pairless {
    public:
        bool operator()(pair<int,int>a,pair<int,int>b){
            return a.second>b.second;
        }
};
class algbase {
    protected:
        vector<int> getrout(int &s, int &t, vector<edge> &edges, vector<int> &pre) {
            vector<int> rout;
            int pp = pre[t];
            while (pp >= 0) {
                rout.push_back(pp);
                pp = pre[edges[pp].s];
            }
            reverse(rout.begin(), rout.end());
            return rout;
        }
    public:
        algbase(){};
        virtual bool cutcake(int)=0;
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo)=0;
	 	virtual vector<int> routalg(int s,int t,int bw)=0;
	 	virtual void prepush(int s,int t,int bw)=0;
};
class dijkstor:public algbase{
    public:
		vector<vector<int>>neibour;
		vector<int>ancestor;
		int edgesize,nodenum,pesize,pnodesize,maxbw;
		vector<vector<vector<int>>>mask;
		vector<edge>edges;
		vector<int>dist;
		vector<int>pre;
		vector<int>leveln;
		vector<int>exn2n;
		vector<vector<int>>rela;
		vector<int>order;
        dijkstor(){};
        void topsort()
        {
        	cout<<" in top sort "<<endl;
        	queue<int>zero;
        	vector<int>orz(nodenum,-1);
        	order=orz;
        	for(int i=0;i<pnodesize;i++)
        		zero.push(i);
        	int biao=0;
			while(!zero.empty())
			{
				int node=zero.front();
				zero.pop();
				order[node]=biao++;
				for(int i=0;i<neibour[node].size();i++)
				{
					if((--ancestor[edges[neibour[node][i]].t])==0)
								zero.push(edges[neibour[node][i]].t);
				}
				cout<<endl;
			}
        }
        virtual bool cutcake(int index){
        	cout<<"cut "<<index<<endl;
        	if(maxbw-(index+1)*10>=0)
        			maxbw-=(index+1)*10;
			else
				{
					cout<<"failure"<<endl;
					return false;
				}
        	vector<int>tmp;
        	for(int i=0;i<edgesize;i++)
        		tmp.push_back(i);
        	mask[index].push_back(tmp);
        	leveln[index]++;
        	return true;
        }
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
        	maxbw=500;
        	rela=relate;
        	edgesize=extenedges.size(),nodenum=ginf.enodesize,pesize=ginf.pesize,pnodesize=ginf.pnodesize;
        	exn2n=ginf.exn2n;
        	edges=extenedges;
        	vector<vector<int>>nd(nodenum,vector<int>());
        	neibour=nd;
        	vector<int>ad(nodenum,0);
        	ancestor=ad;
        	for(int i=0;i<edgesize;i++)
        		{
        			neibour[extenedges[i].s].push_back(extenedges[i].t);
        			ancestor[extenedges[i].t]++;
        		}
    		vector<int>bl(BS,0);
    		leveln=bl;
    		vector<vector<vector<int>>>nm(BS,vector<vector<int>>());
    		mask=nm;
    		vector<int>dd(nodenum*ML,inf);
    		dist=dd;
    		vector<int>pp(nodenum*ML,-1);
    		pre=pp;
    		//topsort();
        }
        virtual vector<int> routalg(int s,int t,int bw){
        	int index=bw/10-1;
        	if(leveln[index]==0)if(!cutcake(index))return vector<int>();
        	int tnode=-1;
        	while(true)
        	{
				vector<int>visited(nodenum,0);
				vector<int>pre(nodenum,-1);
				int vflag=1;
				queue<int>que;
				que.push(s);
	        	visited[s]=1;
				while(!que.empty()&&vflag)
				{
					int node=que.front();
					que.pop();
					for(int i=0;i<neibour[node].size();i++)
					{
						int to=neibour[node][i];
						if(visited[to]==0){pre[to]=node;que.push(to);visited[to]=1;}
						else{continue;}
						if(exn2n[to]==t){tnode=to;vflag=0;break;}
					}

				}
				vector<int>rout;
				int pp=pre[tnode];
				cout<<"tnode is "<<tnode<<endl;
				while(pp>=0)
				{
					rout.push_back(pp);
					pp=pre[pp];
				}
				cout<<endl;
				if(rout.size()>0)
				{
					for(int i=0;i<rout.size();i++)
						cout<<rout[i]<<" ";
					cout<<endl;
					return rout;
				}
				else
				{
					if(!cutcake(index))return vector<int>();
				}
        	}
	 	}
        static bool compare(pair<int,int>&a,pair<int,int>&b)
        {
        	if(a.second<b.second)
        		return true;
        	return false;
        }
        virtual void prepush(int s,int t,int bw)
        {
        	vector<vector<int>>nei(nodenum+2,vector<int>());
        	vector<int>height(nodenum+2,0);
        	vector<int>value(nodenum+2,0);
        	vector<int>weight(pesize,1);
        	for(int i=0;i<edges.size();i++)
        	{
        		nei[edges[i].s].push_back(i);
        		nei[edges[i].t].push_back(i);
        	}
        	/*vector<pair<int,int>>odernode;
        	for(int i=0;i<nodenum;i++)
        		odernode.push_back(make_pair(i,order[i]));
        	sort(odernode.begin(),odernode.end(),compare);*/
        	for(int i=0;i<=WD;i++)
        		height[s+pnodesize*i]=nodenum;
			for(int j=0;j<nei[s].size();j++)
				if(edges[nei[s][j]].s==s)
				{
					weight[nei[s][j]%pesize]*=-1;
					value[edges[nei[s][j]].t]=1;

				}
        	int mark=1;
        	int cc=1;
        	cout<<"before mark"<<endl;
        	while(mark==1)
        	{
        		mark=0;
        		for(int i=0;i<nodenum;i++)
					{
        				if(value[i]>0&&i%pnodesize!=s&&i%pnodesize!=t)
						{
							int flag=0;
							int minheight=INT_MAX;
							for(int j=0;j<nei[i].size();j++)
							{
								int v=(edges[nei[i][j]].s==i)?weight[nei[i][j]%pesize]:-weight[nei[i][j]%pesize];
								int to=(edges[nei[i][j]].s==i)?edges[nei[i][j]].t:edges[nei[i][j]].s;
								if(v>0&&value[i]>0)
								{
									flag=(flag==2)?2:1;
									minheight=min(minheight,height[to]);
									if(height[i]==height[to]+1)
									{
										value[i]--;
										value[to]++;
										weight[nei[i][j]%pesize]*=-1;
										flag=2;
										mark=1;
									}
								}
							}
							if(flag==1)
								{
									height[i]=minheight+1,mark=1;
								}
						}
    				}

        	}
        	cout<<"after mark"<<endl;
        	for(int i=0;i<=WD;i++)
        	{
        		cout<<value[pnodesize*i+t]<<endl;
        	}
        }
};
class parallelor:public algbase
{
	private:
		edge *dev_edges,*aedges;
		int*dev_m,*m,*dev_pre,*pre,*pred,*dev_pred,*dev_d,*d,*dev_mask,*mask,*dev_leveln,*leveln;
		int*dev_rela,*rela;
		int *dev_chan,*chan;
		int presize,dsize,masksize,levelnsize;
		int edgesize,nodenum,pesize,pnodesize;
		int neisize,duansize;
		int *choosel,*dev_choosel;
		int *rout,*dev_rout;
		int *routn,*dev_routn,*order,*dev_order;
		vector<int>hleveln,ancestor;
		int *ordernode,*dev_ordernode;
		int maxbw;
		int *dev_qian,*qian,*dev_qbeg,*qbeg,*dev_qsize,*qsize;
		epair *dev_nei,*nei;
		int *dev_duan,*duan;
		int *dev_beg,*beg;
		int *dev_value,*value;
		int *dev_height,*height;
		vector<vector<int>>neibour;
		void allocate(int maxn,int maxedges);
		void copydata(int s,vector<edge>&edges,int nodenum);
		void dellocate();
	public:
	 	 parallelor();
	 	 void topsort();
	 	 virtual bool cutcake(int index);
	     virtual void prepush(int s,int t,int bw);
	 	 virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo);
	 	 virtual vector<int> routalg(int s,int t,int bw);
	 	 virtual ~parallelor(){dellocate();}
};
#endif //CSPALLOC_PATHALG_H
