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
#define WD 8
#define LY 100
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
            vector<edge> exedges;
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
        	pesize=edges.size();
        	int W=WD+1;
        	nodenum=pnodesize*W;
        	vector<vector<int>>neie(nodenum*LY,vector<int>());
        	vector<vector<int>>nein(nodenum*LY,vector<int>());
        	vector<int>height(nodenum*LY,0);
        	vector<int>value(nodenum*LY,0);
        	vector<int>esign(pesize*LY,1);
        	for(int k=0;k<LY;k++)
        	{
        		int startn=k*nodenum;
        		int starte=k*pesize;
        		for(int i=0;i<edges.size();i++)
        			for(int j=0;j<W-1;j++)
        			{
        				int s=edges[i].s*W+j+startn;
        				int t=edges[i].t*W+j+1+startn;
        				neie[s].push_back(i+1+starte);
        				neie[t].push_back(-(i+1+starte));
        				nein[s].push_back(t);
        				nein[t].push_back(s);
        			}
        	}
        	for(int k=0;k<LY;k++)
        		{
        		int startn=k*nodenum;
        		for(int i=0;i<W;i++)
        			height[startn+s+pnodesize*i]=W+1;
        		}
        	for(int k=0;k<LY;k++)
        	{
        		int startn=k*nodenum;
        		for(int i=1;i<W;i++)
        			height[startn+s+pnodesize*i]=INT_MAX;
        	}
        	for(int k=0;k<LY;k++)
        			{
        			for(int i=0;i<edges.size();i++)
        				if(edges[i].s==s)
        					{
        					value[k*W*pnodesize+W*edges[i].t+1]=1;
        					esign[k*edges.size()+i]*=-1;
        					}
        			}
        	int mark=1;
        	int cc=1;
        	cout<<"before mark here"<<endl;
        	while(mark==1)
        	{
        		mark=0;
        		for(int i=0;i<nodenum*LY;i++)
					{
        				if(value[i]>0&&i%pnodesize!=s&&i%pnodesize!=t)
						{
							int minheight=INT_MAX;
							for(int j=0;j<nein[i].size();j++)
							{
								if(value[i]>0&&(esign[abs(neie[i][j])-1]*neie[i][j])>0)
								{
									int to=nein[i][j];
									if(height[i]==height[to]+1)
									{
										value[i]--;
										value[to]++;
										esign[abs(neie[i][j])-1]*=-1;
										mark=1;
									}
									minheight=min(minheight,height[to]);
								}
							}
							if(value[i]>0&&minheight<INT_MAX)
								height[i]=minheight+1,mark=1;
						}
    				}
        	}
        	cout<<"after mark"<<endl;
        	for(int i=0;i<LY*W*pnodesize;i++)
        		if(value[i]!=0)
        			{
        				int bi=i%nodenum;
        				cout<<i/nodenum<<" "<<bi<<" "<<bi/W<<" "<<bi%W<<" "<<height[i]<<" "<<value[i]<<endl;
        			}
        	int count=0;
        	for(int i=0;i<edges.size()*LY;i++)
        		if(esign[i]<0)
        			count++;
        	cout<<"count is: "<<count<<endl;
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
		vector<edge>edges;
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
