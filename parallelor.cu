
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include"pathalg.h"
static const int WORK_SIZE = 256;
__global__ void BFShigh(int t,int *m,int index,epair*nei,int *d,int *chan,int edgesize,int tedgesize,int round,int pnodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=tedgesize)return;
	int from=nei[i].f;
	if (chan[from]<0)return;
	chan[from]=-1;
	int to=nei[i].t;
	d[to]=round;
	if((to%pnodenum)==t)*m=1;
}
__global__ void BFShighN(int t,int *m,int index,epair*nei,int* duan,int*beg,int *d,int *chan,int round,int pnodenum,int nodenum)
{
	int from=threadIdx.x + blockIdx.x*blockDim.x;
	if(from>=nodenum)return;
	if (chan[from]<0)return;
	for(int k=beg[from];k<(beg[from]+duan[from]);k++)
	{
		int to=nei[k].t;
		d[to]=round;
		if((to%pnodenum)==t)*m=1;
	}
}
__global__ void initchan(int s,int *chan,int *d,int *pred,int nodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nodenum)return;
	chan[i]=(i==s)?1:-1;
	d[i]=(i==s)?0:inf;
	pred[i]=d[i];
}
__global__ void chanchan(int *m,int *pred,int *d,int *chan,int totalsize,int nodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=totalsize)return;
	chan[i]=-1;
	if(d[i]<pred[i])
	{
		chan[i]=1;
		pred[i]=d[i];
	}
}
void parallelor::copydata(int s,vector<edge>&edges,int nodenum){
	memset(pre,-1,sizeof(int)*nodenum);
	*m=0;
	for(int i=0;i<nodenum;i++)
		d[i]=INT_MAX/2;
	d[s]=0;
	for(int i=0;i<edges.size();i++)
		aedges[i]=edges[i];
	cudaMemcpy(dev_edges,aedges,edges.size()* sizeof(edge),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_m,m,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d,d,sizeof(int)*nodenum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pre,pre,sizeof(int)*nodenum,cudaMemcpyHostToDevice);
};
void parallelor::dellocate(){
	/*delete[]d;
	delete[]pre;
	delete[]aedges;
	delete m;
	cudaFree(dev_edges);
	cudaFree(dev_m);
	cudaFree(dev_d);
	cudaFree(dev_pre);*/
};
void parallelor::allocate(int maxn,int maxedge){
	m=new int;
	d=new int[maxn],pre=new int[maxn];
	aedges=new edge[maxedge];
	cudaMalloc(&dev_edges, sizeof(edge)*maxedge);
	cudaMalloc((void**)&dev_d,maxn*sizeof(int));
	cudaMalloc((void**)&dev_pre,maxn*sizeof(int));
	cudaMemcpy(duan,dev_duan,duansize*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMalloc((void**)&dev_m,sizeof(int));
}
bool parallelor::cutcake(int index){

	cout<<"cut "<<index<<endl;
	if(maxbw-(index+1)*10>=0)
		maxbw-=(index+1)*10;
	else
		{
			cout<<"failure"<<endl;
			return false;
		}
	hleveln[index]++;
	return true;

};
void parallelor::topsort()
{
       	cout<<" in top sort "<<endl;
       	cout<<"node num is "<<nodenum<<endl;
       	queue<int>zero;
       	order=new int[nodenum];
       	ordernode=new int[nodenum];
       	for(int i=0;i<nodenum;i++)
       		if(ancestor[i]==0)
       			zero.push(i);
       	int biao=0;
			while(!zero.empty())
			{
				int node=zero.front();
				zero.pop();
				order[node]=biao++;
				ordernode[biao-1]=node;
				for(int i=0;i<neibour[node].size();i++)
				{
					if((--ancestor[aedges[neibour[node][i]].t])==0)
							zero.push(aedges[neibour[node][i]].t);
				}
			}
		cout<<biao<<" "<<nodenum<<endl;
		cout<<"out top"<<endl;
};
void parallelor::init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
	cout<<"in cuda init"<<endl;
	maxbw=500;
	//allocate in cuda
	edgesize=extenedges.size();nodenum=ginf.enodesize;
	edges=extenedges;
	pesize=ginf.pesize;pnodesize=ginf.pnodesize;
	/*cout<<"es "<<edgesize<<"  pes"<<pesize<<endl;
	dsize=ML*nodenum,presize=ML*nodenum;
	neisize=BS*ML*edgesize;
	duansize=nodenum;
	vector<vector<int>>nd(nodenum,vector<int>());
    neibour=nd;
    vector<int>as(nodenum,0);
    ancestor=as;
	for(int i=0;i<edgesize;i++)
		{
			neibour[extenedges[i].s].push_back(i);
			ancestor[extenedges[i].t]++;
		}
	levelnsize=BS;
	cudaMalloc(&dev_edges, sizeof(edge)*edgesize);
	cudaMalloc((void**)&dev_d,dsize*sizeof(int));
	cudaMalloc((void**)&dev_pred,dsize*sizeof(int));
	cudaMalloc((void**)&dev_pre,presize*sizeof(int));
	cudaMalloc((void**)&dev_chan,presize*sizeof(int));
	cudaMalloc((void**)&dev_m,sizeof(int));
	cudaMalloc((void**)&dev_choosel,sizeof(int));
	cudaMalloc((void**)&dev_nei,neisize*sizeof(epair));
	cudaMalloc((void**)&dev_rela,(WD+1)*edgesize*sizeof(int));
	cudaMalloc((void**)&dev_rout,WD*sizeof(int));
	cudaMalloc((void**)&dev_routn,sizeof(int));
	cudaMalloc((void**)&dev_duan,duansize*sizeof(int));
	cudaMalloc((void**)&dev_beg,duansize*sizeof(int));
	cudaMalloc((void**)&dev_order,nodenum*sizeof(int));
	cudaMalloc((void**)&dev_ordernode,nodenum*sizeof(int));
	cudaMalloc((void**)&dev_qian,edgesize*sizeof(int));
	cudaMalloc((void**)&dev_qsize,nodenum*sizeof(int));
	cudaMalloc((void**)&dev_qbeg,nodenum*sizeof(int));
	//new in host ;
	aedges=new edge[edgesize];
	choosel=new int;
	m=new int;
	pred=new int[dsize];
	d=new int[dsize],pre=new int[presize],chan=new int[presize];
	leveln=new int[levelnsize];
	rela=new int[(WD+1)*edgesize];
	nei=new epair[neisize];
	duan=new int[duansize];
	beg=new int[duansize];
	rout=new int[WD];
	for(int i=0;i<WD;i++)
		rout[i]=-1;
	routn=new int;
	//init in host ;
	*m=0;
	*choosel=0;
	memset(pre,-1,sizeof(int)*presize);
	memset(chan,0,sizeof(int)*presize);
	for(int i=0;i<dsize;i++)
		d[i]=inf,pred[i]=inf;
	for(int i=0;i<edgesize;i++)
		aedges[i]=extenedges[i];
	for(int i=0;i<relate.size();i++)
		for(int j=0;j<WD+1;j++)
			if(j<relate[i].size())
				rela[i*WD+j]=relate[i][j];
			else
				rela[i*WD+j]=-1;
	memset(leveln,0,sizeof(int)*levelnsize);
	int h=0;
	topsort();
	vector<vector<int>>vqian(nodenum,vector<int>());
	int g=0;
	for(int i=0;i<BS*ML;i++)
		{
			int t=0;
			for(int j=0;j<nodenum;j++)
				{
					beg[j]=t;
					duan[j]=neibour[ordernode[j]].size();
					t+=neibour[ordernode[j]].size();
					for(int k=0;k<neibour[ordernode[j]].size();k++)
						{
							nei[h].f=ordernode[j];
							nei[h].t=aedges[neibour[ordernode[j]][k]].t;
							if(i==0)
								vqian[nei[h].t].push_back(nei[h].f);
							h++;
						}
				}
		}
	qian=new int[edgesize];
	qsize=new int[nodenum];
	qbeg=new int[nodenum];
	int y=0;
	int ss=0;
	for(int i=0;i<vqian.size();i++)
		{
			qsize[i]=vqian[i].size();
			qbeg[i]=ss;
			for(int j=0;j<vqian[i].size();j++)
				qian[y++]=vqian[i][j];
			ss+=vqian[i].size();
		}
	cudaMemcpy(dev_edges,aedges,edgesize* sizeof(edge),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_m,m,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_choosel,choosel,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d,d,dsize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pred,pred,dsize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pre,pre,presize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_chan,chan,presize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_rela,rela,(WD+1)*edgesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_routn,m,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_rout,rout,WD*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nei,nei,neisize*sizeof(epair),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_duan,duan,duansize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_beg,beg,duansize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_order,order,nodenum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_ordernode,ordernode,nodenum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_qian,qian,edgesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_qsize,qsize,nodenum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_qbeg,qbeg,nodenum*sizeof(int),cudaMemcpyHostToDevice);
	//for(int i=0;i<nodenum;i++)
		//cout<<order[i]<<endl;
	vector<int>tmp(levelnsize,0);
	hleveln=tmp;
	cout<<"out init"<<endl;*/
};
parallelor::parallelor()
{

};
vector<int> parallelor:: routalg(int s,int t,int bw)
{
	int zero=0;
	int index=bw/10-1;
	if(hleveln[index]<=0)cutcake(index);
	int max=0;
	int h=10;
	t=ordernode[t];
	s=ordernode[s];
	cout<<"blasting "<<endl;
	//while(true)
	{
		initchan<< <(nodenum/WORK_SIZE)+1, WORK_SIZE >> >(s,dev_chan,dev_d,dev_pred,nodenum);
		int kk=1,gg=8;
		cudaMemcpy(dev_m, &zero, sizeof(int), cudaMemcpyHostToDevice);
		do{
			/*cudaMemcpy(chan,dev_chan,nodenum*sizeof(int), cudaMemcpyDeviceToHost);
			int cc=0;
			for(int i=0;i<nodenum;i++)
				if(chan[i]>=0)
					cc++;
			cout<<cc<<endl;*/
			BFShigh << <(edgesize/WORK_SIZE)+1, WORK_SIZE >> >(t,dev_m,index,dev_nei,dev_d,dev_chan,edgesize,edgesize,kk,pnodesize);
			//BFShighN<< <(nodenum/WORK_SIZE)+1, WORK_SIZE >> >(t,dev_m,index,dev_nei,dev_duan,dev_beg,dev_d,dev_chan,kk,pnodesize,nodenum);
			chanchan<< <(nodenum/WORK_SIZE)+1, WORK_SIZE >> >(dev_m,dev_pred,dev_d,dev_chan,nodenum,nodenum);
			cudaMemcpy(m, dev_m, sizeof(int), cudaMemcpyDeviceToHost);
			kk++;
		}
		while(*m==0);
		cout<<"kk is: "<<kk<<endl;
		/*cudaMemcpy(d, dev_d, sizeof(int)*nodenum, cudaMemcpyDeviceToHost);
		int k=0;
		while(t<nodenum)
		{
			k++;81
135
1449
2673
4896
			cout<<d[t]<<endl;
			t+=pnodesize;
		}
		cout<<"over "<<endl;
		/*cudagetrout<< <1,1>> >(dev_qian,dev_qsize,dev_qbeg,dev_d,s,t,dev_rout,dev_routn,dev_choosel,1,nodenum,pnodesize);
		cudaMemcpy(routn,dev_routn,sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(rout,dev_rout,WD*sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(choosel,dev_choosel,sizeof(int),cudaMemcpyDeviceToHost);
		cout<<"size is "<<*routn<<endl;
		if(*routn==0)
			{
				if(!cutcake(index))
					return vector<int>();
			}
		else
		{
			cout<<(index+1)*10<<"/"<<*choosel<<": ";
			for(int i=0;i<*routn;i++)
				cout<<rout[i]<<" ";
			cout<<endl;
			return vector<int>();
		}*/
	}
	return vector<int>();
};
int fls(int x)
{
	int position;
	int i;
	if(x!=0)
		for(i=(x>>1),position=0;i!=0;++position)
			i>>=1;
	else
		position=-1;
	return pow(2,position+1);
}
/*__global__ void push(int*dev_h,int*dev_v,int*dev_ev,int*dev_s,int*dev_t,int E,int W,int *mark,int sorce,int end)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int id=threadIdx.x;
	__shared__  int flows[WORK_SIZE];
	__shared__  int biao[WORK_SIZE];
	int eid=i/W;
	if(eid>E)return;
	int offset=i%W;
	biao[id]=offset;
	int s=dev_s[eid]*(W+1)+offset;
	int t=dev_t[eid]*(W+1)+offset;
	if(dev_ev[eid]>0&&dev_v[s]>0&&s!=sorce&&dev_s[eid]!=end)
	{
		if(dev_h[s]==dev_h[t]+1)
			flows[id]=1,dev_ev[eid]*=-1;
	}
	if(dev_ev[eid]<0&&dev_v[t]>0&&t!=sorce&&dev_t[eid]!=end)
	{
		if(dev_h[s]+1==dev_h[t])
			flows[id]=2,dev_ev[eid]*=-1;
	}
	int start=(id/W)*W;
	for(int d=W;d>1;d=d/2)
	{
		if(id-start<d/2)
			if(flows[id]<flows[id+d/2])
				flows[id]=flows[id+d/2],biao[id]=biao[id+d/2];
	}
	if(i%W==0)
	{
		if(flows[id]==1)atomicAdd(&dev_v[t+biao[id]],1),atomicSub(&dev_v[s+biao[id]],1),*mark=1;
		if(flows[id]==2)atomicAdd(&dev_v[s+biao[id]],1),atomicSub(&dev_v[t+biao[id]],1),*mark=1;
	}
};*/

__global__ void push(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=N||dev_v[i]==0||i/W==s||i/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int offset=i%W;
	int value=dev_v[i];
	for(int j=0;j<max;j++)
	{
		if(dev_nein[b+j]<INT_MAX&&value>0)
		{
			if(h==dev_h[dev_nein[b+j]]+1&&dev_neie[b+j]*dev_esign[abs(dev_neie[b+j])]>0)
				{
					if(dev_neie[b+j]>0)
						dev_emark[abs(dev_neie[b+j])]=dev_nein[b+j];
					else
						dev_emark[abs(dev_neie[b+j])]=-i;
					value--;
				}
		}
		else
			break;
	}
};
__global__ void aggregate(int* dev_esign,int*dev_v,int* dev_emark,int W,int E,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E)return;
	int id=threadIdx.x;
	int eiid=id/W;
	int start=eiid*W;
	int eid=i/W;
	__shared__ int emids[500];
	emids[eiid]=INT_MAX;
	__syncthreads();
	if(dev_emark[i]<INT_MAX)
		emids[eiid]=dev_emark[i];
	__syncthreads();
	int emid=emids[eiid];
	if(id==start&&emid<INT_MAX)
		{
			int s=abs(dev_esign[eid])*(W+1)+abs(emid)%(W+1);
			if(emid>0)
			{	atomicSub(&dev_v[s],1);
				atomicAdd(&dev_v[abs(emid)],1);
				*mark=1;
				dev_esign[eid]*=-1;
			}
			else
			{
				atomicAdd(&dev_v[s],1);
				atomicSub(&dev_v[abs(emid)],1);
				*mark=1;
				dev_esign[eid]*=-1;
			}
		}
	dev_emark[i]=INT_MAX;
};
__global__ void aggregate1(int* dev_esign,int*dev_v,int* dev_emark,int W,int E,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E)return;
	__shared__  int emark[WORK_SIZE];
	int id=threadIdx.x;
	emark[id]=dev_emark[i];
	int start=(id/W)*W;
	int eid=i/W;
	for(int d=W;d>1;d=d/2)
		if(id-start<d/2)
			if(abs(emark[id])>abs(emark[id+d/2]))
				emark[id]=emark[id+d/2];
	if(id==start&&emark[id]<INT_MAX)
		{
			int s=abs(dev_esign[eid])*(W+1)+abs(emark[id])%(W+1);
			if(emark[id]>0)
			{	atomicSub(&dev_v[s],1);
				atomicAdd(&dev_v[abs(emark[id])],1);
				*mark=1;
				dev_esign[eid]*=-1;
			}
			else
			{
				atomicAdd(&dev_v[s],1);
				atomicSub(&dev_v[abs(emark[id])],1);
				*mark=1;
				dev_esign[eid]*=-1;
			}
		}
	dev_emark[i]=INT_MAX;
};
__global__ void aggregate2(int* dev_esign,int*dev_v,int* dev_emark,int W,int E,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E)return;
	int emid=dev_emark[i];
	if(emid<INT_MAX)
	{
		int s=abs(dev_esign[i])*(W+1)+abs(emid)%(W+1);
		if(emid>=0)
		{	atomicSub(&dev_v[s],1);
			atomicAdd(&dev_v[abs(emid)],1);
			*mark=1;
			dev_esign[i]*=-1;
		}
		else
		{
			atomicAdd(&dev_v[s],1);
			atomicSub(&dev_v[abs(emid)],1);
			*mark=1;
			dev_esign[i]*=-1;
		}
	}
	dev_emark[i]=INT_MAX;
};
__global__ void relable(int*dev_h,int*dev_v,int N,int*mark,int*dev_nein,int*dev_neie,int *dev_esign,int max,int W,int s,int t)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=N||dev_v[i]==0||i/W==s||i/W==t)return;
	int b=i*max;
	int mini=INT_MAX;
	for(int j=0;j<max;j++)
	{
		if(dev_nein[b+j]<INT_MAX)
		{
			if(dev_neie[b+j]*dev_esign[abs(dev_neie[b+j])]>0)
				mini=min(mini,dev_h[abs(dev_nein[b+j])]);
		}
		else
			break;
	}
	if(mini!=INT_MAX)
		dev_h[i]=mini+1,*mark=1;

};
void parallelor::prepush(int s,int t,int bw)
{
	cout<<"begin"<<endl;
	int W=WD+1;
	int*dev_h,*h,*dev_v,*v,*dev_neie,*neie,*dev_nein,*nein;
	int*dev_esign,*esign;
	int *dev_emark,*emark,*mark,*dev_mark;
	int *minarray;
	h=new int[W*pnodesize];
	v=new int[W*pnodesize];
	minarray=new int[pnodesize];
	mark=new int;
	vector<vector<int>>rawneie(W*pnodesize,vector<int>());
	vector<vector<int>>rawnein(W*pnodesize,vector<int>());
	for(int i=0;i<edges.size();i++)
		for(int j=0;j<W-1;j++)
			{
				int s=edges[i].s*W+j;
				int t=edges[i].t*W+j;
				rawneie[s].push_back(i);
				rawneie[t].push_back(-i);
				rawnein[s].push_back(t);
				rawnein[t].push_back(s);
			}
	int max=0;
	for(int i=0;i<rawnein.size();i++)
		if(rawnein[i].size()>max)max=rawnein[i].size();
	max++;
	neie=new int[W*pnodesize*max];
	nein=new int[W*pnodesize*max];
	for(int i=0;i<W*pnodesize;i++)
	{
		for(int j=0;j<max;j++)
		{
			if(j<rawneie[i].size())
				{
					neie[i*max+j]=rawneie[i][j];
					nein[i*max+j]=rawnein[i][j];
				}
			else
				{
					neie[i*max+j]=INT_MAX;
					nein[i*max+j]=INT_MAX;
				}
		}
	}
	emark=new int[edges.size()];
	esign=new int[edges.size()];
	for(int i=0;i<edges.size();i++)
		emark[i]=INT_MAX;
	for(int i=0;i<edges.size();i++)
		esign[i]=edges[i].s;
	for(int i=0;i<W*pnodesize;i++)
		{
			h[i]=0;
			v[i]=0;
		}
	for(int i=0;i<edges.size();i++)
		if(edges[i].s==s)
			{
				v[W*edges[i].t+1]=1;
				esign[i]*=-1;
			}
	for(int i=s*W;i<s*W+W;i++)
		h[i]=WD+1;
	cudaMalloc((void**)&dev_h,W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_mark,sizeof(int));
	cudaMalloc((void**)&dev_v,W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_neie,W*max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_nein,W*max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_esign,edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_emark,edges.size()*sizeof(int));
	cudaMemcpy(dev_h,h,W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v,v,W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neie,neie,W*max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nein,nein,W*max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_esign,esign,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_emark,emark,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	*mark=1;
	int time=0;
	while(*mark>0)
	{
		*mark=0;
		cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
		push<<<W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,dev_esign,dev_emark,dev_neie,dev_nein,W*pnodesize,max,W,s,t);
		//aggregate<<<edges.size()*(W-1)/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W-1,edges.size()*(W-1),dev_mark);
		aggregate2<<<edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W-1,edges.size(),dev_mark);
		relable<<<W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,W*pnodesize,dev_mark,dev_nein,dev_neie,dev_esign,max,W,s,t);
		cudaMemcpy(mark,dev_mark,sizeof(int),cudaMemcpyDeviceToHost);
		time++;
	}
	cout<<"times is :"<<time<<endl;
	cudaMemcpy(v,dev_v,W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(h,dev_h,W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	for(int i=0;i<W*pnodesize;i++)
		if(v[i]!=0)
			cout<<i<<" "<<i/W<<" "<<h[i]<<" "<<v[i]<<endl;
	delete[] h;
	delete[] minarray;
	delete[] v;
	delete[] mark;
	delete[] neie;
	delete[] nein;
	delete[]emark;
	delete[]esign;
	cudaFree(dev_h);
	cudaFree(dev_mark);
	cudaFree(dev_v);
	cudaFree(dev_neie);
	cudaFree(dev_nein);
	cudaFree(dev_esign);
	cudaFree(dev_emark);
};
