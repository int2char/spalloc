#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include"pathalg.h"
static const int WORK_SIZE = 512;
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
	cout<<"out cuda init"<<endl;
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
	while(true)
	for(int i=0;i<1;i++)
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
			k++;
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
__global__ void push(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]>INT_MAX/2&&h==hnbj+1)
				{
					dev_emark[eid]=(ebj>0)?nbj:i;
					value--;
					*mark=1;
				}
				minheight=min(minheight,hnbj);
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX){dev_h[i]=minheight+1;*mark=1;}
};
__global__ void push1(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]==0&&h==hnbj+1)
				{
					dev_emark[eid]++;
					value--;
					*mark=1;
				}
				minheight=min(minheight,hnbj);
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX){dev_h[i]=minheight+1;*mark=1;}
};
__global__ void pushrelable(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]==i)
				{
					dev_emark[eid]++;
					atomicSub(&dev_v[i],1);
					atomicAdd(&dev_v[nbj],1);
					value--;
					dev_esign[eid]*=-1;
					*mark=1;
				}
				minheight=min(minheight,hnbj);
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX){dev_h[i]=minheight+1;*mark=1;}
};
__global__ void aggregate3(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
		int s,t;
	if(dev_emark[i]>0)
	{
		if(dev_esign[i]>0)
		{
			s=dev_st[i];
			t=dev_te[i]+1;
		}
		else
		{
			t=dev_st[i];
			s=dev_te[i]+1;
		}
		for(int k=0;k<W;k++)
			{
				int h1=dev_h[s+k];
				int h2=dev_h[t+k];
				if(dev_v[s+k]>0&&h1==h2+1)
				{
					atomicSub(&dev_v[s+k],1);
					atomicAdd(&dev_v[t+k],1);
					dev_esign[i]*=-1;
					break;
				}
			}
	}
	dev_emark[i]=0;
};
__global__ void aggregate4(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
	int s,t;
	dev_emark[i]=INT_MAX;
	if(dev_esign[i]>0)
	{
		s=dev_st[i];
		t=dev_te[i]+1;
	}
	else
	{
		t=dev_st[i];
		s=dev_te[i]+1;
	}
	for(int k=0;k<W;k++)
		{
			int h1=dev_h[s+k];
			int h2=dev_h[t+k];
			if(dev_v[s+k]>0&&h1==h2+1)
			{
				dev_emark[i]=s+k;
				break;
			}
		}
};
__global__ void aggregate2(int* dev_esign,int*dev_v,int* dev_emark,int W,int E)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
	int emid=dev_emark[i];
	if(emid<INT_MAX)
	{
		int s=abs(dev_esign[i])-2+emid%W;
		if(dev_esign[i]>0)
		{	atomicSub(&dev_v[s],1);
			atomicAdd(&dev_v[emid],1);
		}
		else
		{	atomicAdd(&dev_v[s],1);
			atomicSub(&dev_v[emid],1);
		}
		dev_esign[i]*=-1;
	}
	dev_emark[i]=INT_MAX;
};
__global__ void aggregate5(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E)
{
        int i = threadIdx.x + blockIdx.x*blockDim.x;
        if(i>=E*LY*W)return;
        int s,t;
        int eid=i/W;
        int k=i%W;
        if(dev_esign[eid]>0)
        {
                s=dev_st[eid];
                t=dev_te[eid]+1;
        }
        else
        {
                t=dev_st[eid];
                s=dev_te[eid]+1;
        }
        int h1=dev_h[s+k];
        int h2=dev_h[t+k];
        if(dev_v[s+k]>0&&h1==h2+1)
                dev_emark[eid]=s+k;
};

__global__ void relable(int*dev_h,int*dev_v,int N,int*mark,int*dev_nein,int*dev_neie,int *dev_esign,int max,int W,int s,int t)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	if(i>=N*LY||dev_v[i]==0||bi/W==s||bi/W==t)return;
	int b=i*max;
	int mini=INT_MAX;
	for(int j=0;j<max;j++)
	{
		if(dev_nein[b+j]<INT_MAX)
		{
			if((dev_neie[b+j]^dev_esign[abs(dev_neie[b+j])-1])>0)
				mini=min(mini,dev_h[dev_nein[b+j]]);
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
	int *dev_st,*dev_te;
	int* st,*te;
	h=new int[W*pnodesize*LY];
	v=new int[W*pnodesize*LY];
	mark=new int;
	vector<vector<int>>rawneie(W*pnodesize,vector<int>());
	vector<vector<int>>rawnein(W*pnodesize,vector<int>());
	for(int i=0;i<edges.size();i++)
		for(int j=0;j<W-1;j++)
		{
			int s=edges[i].s*W+j;
			int t=edges[i].t*W+j+1;
			rawneie[s].push_back(i+1);
			rawneie[t].push_back(-(i+1));
			rawnein[s].push_back(t);
			rawnein[t].push_back(s);
		}
	int max=0;
	for(int i=0;i<rawnein.size();i++)
		if(rawnein[i].size()>max)max=rawnein[i].size();
	max++;
	neie=new int[W*pnodesize*max*LY];
	nein=new int[W*pnodesize*max*LY];
	for(int k=0;k<LY;k++)
		{
		int startn=k*W*pnodesize;
		int starte=k*edges.size();
		for(int i=0;i<W*pnodesize;i++)
			{
				for(int j=0;j<max;j++)
				{
					if(j<rawneie[i].size())
						{
							if(rawneie[i][j]<0)
								neie[(startn+i)*max+j]=rawneie[i][j]-starte;
							else
								neie[(startn+i)*max+j]=rawneie[i][j]+starte;
							nein[(startn+i)*max+j]=rawnein[i][j]+startn;
						}
					else
						{
							neie[(startn+i)*max+j]=INT_MAX;
							nein[(startn+i)*max+j]=INT_MAX;
						}
				}
			}
		}
	emark=new int[LY*edges.size()];
	esign=new int[LY*edges.size()];
	st=new int[LY*edges.size()];
	te=new int[LY*edges.size()];
	for(int i=0;i<LY*edges.size();i++)
		emark[i]=INT_MAX;//0
	for(int k=0;k<LY;k++)
		for(int i=0;i<edges.size();i++)
			{
				st[i+k*edges.size()]=edges[i].s*W+k*W*pnodesize;
				te[i+k*edges.size()]=edges[i].t*W+k*W*pnodesize;
				esign[i+k*edges.size()]=edges[i].s*W+k*W*pnodesize+1;
			}
	for(int i=0;i<W*LY*pnodesize;i++)
		{
			h[i]=0;
			v[i]=0;
		}
	for(int k=0;k<LY;k++)
		{
		for(int i=0;i<edges.size();i++)
			if(edges[i].s==s)
				{
				v[k*W*pnodesize+W*edges[i].t+1]=1;
				esign[k*edges.size()+i]*=-1;
				}
		}
	for(int k=0;k<LY;k++)
		{
		int start=k*W*pnodesize;
		for(int i=s*W;i<s*W+W;i++)
			h[i+start]=WD;
		}
	cudaMalloc((void**)&dev_h,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_mark,sizeof(int));
	cudaMalloc((void**)&dev_v,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_neie,LY*W*max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_nein,LY*W*max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_esign,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_emark,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_st,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_te,LY*edges.size()*sizeof(int));
	cudaMemcpy(dev_h,h,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v,v,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neie,neie,LY*W*max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nein,nein,LY*W*max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_esign,esign,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_emark,emark,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_st,st,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_te,te,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	*mark=1;
	int time=0;
	cout<<"max is "<<max<<endl;
	time_t start,end;
	start=clock();
	while(*mark>0)
	//for(int i=0;i<3;i++)
	{
		//cout<<"************"<<endl;
		*mark=0;
		cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
		push<<<LY*W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,dev_esign,dev_emark,dev_neie,dev_nein,W*pnodesize,max,W,s,t,dev_mark);
		aggregate2<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W,edges.size());
		//aggregate3<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,dev_st,dev_te,dev_h,W,edges.size());
		/*cudaMemcpy(emark,dev_emark,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
		for(int i=0;i<LY*edges.size();i++)
			if(emark[i]>0)
				cout<<"gota... "<<i<<"s:"<<st[i]<<" "<<te[i]<<" "<<emark[i]<<endl;*/
		//relable<<<LY*W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,W*pnodesize,dev_mark,dev_nein,dev_neie,dev_esign,max,W,s,t);
		//aggregate2<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W,edges.size(),W*pnodesize,dev_mark);
		cudaMemcpy(mark,dev_mark,sizeof(int),cudaMemcpyDeviceToHost);
		/*cudaMemcpy(v,dev_v,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(h,dev_h,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(esign,dev_esign,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
		int flow=0;
		for(int i=0;i<LY*W*pnodesize;i++)
			if(v[i]!=0)
				{
					int bi=i%(W*pnodesize);
					if(bi/W==t)flow+=v[i];
					cout<<i/(W*pnodesize)<<" "<<bi<<" "<<bi/W<<" "<<bi%W<<" "<<h[i]<<" "<<v[i]<<endl;
					/*if(i==319)
					{
						for(int j=0;j<max;j++)
							if(nein[i*max+j]<INT_MAX)
								cout<<neie[i*max+j]<<" "<<esign[abs(neie[i*max+j])-1]<<" "<<h[nein[i*max+j]]<<endl;
					}
				}*/
		time++;
	}
	end=clock();
	cout<<"GPU time is: "<<end-start<<endl;
	cudaMemcpy(v,dev_v,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(h,dev_h,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	int flow=0;
	for(int i=0;i<LY*W*pnodesize;i++)
		if(v[i]!=0)
			{
				int bi=i%(W*pnodesize);
				if(bi/W==t)flow+=v[i];
				//cout<<i/(W*pnodesize)<<" "<<bi<<" "<<bi/W<<" "<<bi%W<<" "<<h[i]<<" "<<v[i]<<endl;
			}
	cudaMemcpy(esign,dev_esign,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
	int count=0;
	for(int i=0;i<edges.size()*LY;i++)
		if(esign[i]<0)
			count++;
	cout<<"resort"<<endl;
	/*for(int i=0;i<edges.size()*LY;i++)
		{
			if(esign[i]<0)
			{
				int sorce=edges[i].t*W;
				if(sorce/W==t)
				{
					int pre=edges[i].s*W;
					cout<<pre<<" ";
					while((pre/W)!=s)
					{
						int flag=0;
						for(int h=0;h<W;h++)
						{
							pre++;
							for(int k=0;k<max;k++)
								{
									if(nein[pre*max+k]<INT_MAX)
										if(esign[abs(neie[pre*max+k])-1]<0&&neie[pre*max+k]<0)
										{
											esign[abs(neie[pre*max+k])-1]*=-1;
											pre=edges[abs(neie[pre*max+k])-1].s*W;
											cout<<pre<<" ";
											flag=1;
										}
										if(flag==1)break;
								}
							if(flag==1)break;
						}
					}
					cout<<endl;
				}
			}
		}*/
	cout<<"flow is"<<flow<<endl;
	cout<<"count is "<<count<<endl;
	cout<<"die is "<<time<<endl;
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
