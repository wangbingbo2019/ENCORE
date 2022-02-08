#coding:gbk
'''
Created on 2019-8-15

@author: Bingbo Wang
'''

# coding: utf-8

# In[ ]:

import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
G2 = nx.Graph()
def Pattern_detection(filename):
    global G2
    DEfilename=filename
    openPPI(r'newHnet-2015name')
    (genenumber,fc)=readdata(DEfilename)

    (zlist,fclist,size)=lcc_zscore(genenumber,fc,100,50,max(fc))
    max_fcvalue=refine(zlist,fclist,size) 
    
    (zlist,fclist,size)=lcc_zscore(genenumber,fc,1000,50,max_fcvalue)
    
    (p1,fccutoff1,fccutoff2,minpick,sig)=leastsqfit(fclist,zlist)
   
    print('FC cutoffs for disease neighborhood and core:',pow(2,fccutoff1),pow(2,fccutoff2))
    
    DrawPicture(DEfilename,pow(2,max_fcvalue),max(zlist)+2,p1,zlist,fclist,fc)
    
    (module,coremodule)=getlist(fccutoff1,fccutoff2,genenumber,fc)
   
    core=getLCC(G2,coremodule)
    Dmodule=getLCC(G2,module)
    periphery=list(set(Dmodule)-set(core))
    g=nx.subgraph(G2,Dmodule)
    gg=nx.subgraph(G2,core)
    f=open(r'disease neighborhood.txt','w')
    for (x,y) in g.edges():
        f.write('%s'%G2.node[x]['label'])
        f.write('\t')
        f.write('%s'%G2.node[y]['label'])
        f.write('\n')
    f.close()
    f=open(r'core and periphery genes.txt','w')
    f.write('GeneSymbol')
    f.write('\t')
    f.write('Categories')
    f.write('\n')
    for x in core:
        f.write('%s'%G2.node[x]['label']) 
        f.write('\t')
        f.write('core')
        f.write('\n')
    for x in periphery:
        f.write('%s'%G2.node[x]['label']) 
        f.write('\t')
        f.write('periphery')
        f.write('\n')
    f.close()
    return(g,gg)
    
def openPPI(filename):
    global G2
    a=open(filename,"r")    
    for i in a:
        n=i.strip().split("\t")
        G2.add_edge(int(n[0]),int(n[1]))
        G2.add_node(int(n[0]),label=n[2])
        G2.add_node(int(n[1]),label=n[3])
    a.close()
def readdata(filename):
    global G2
    p=0.05
    
    fcalfa=0.01
    genename=[]
    genenumber=[]
    fc=[]
    pvalue=[]
    number=list(G2.nodes())
    first=open(filename,"r")
    for i in first:
        n=i.strip().split("\t")
        if n[0]=='Gene':
            continue
        
        if(float(n[1])<=p):
            genename.append(n[0])
            fc.append(abs(float(n[2])))
            pvalue.append(float(n[1]))
            flag=0
            for x in range(len(G2.node)):
                if G2.node[number[x]]['label']==n[0]:
                    genenumber.append(number[x])
                    flag=1
                    break
            if flag==0:
                genenumber.append(0)
    print('Number of DE genes with p-value<0.05:',len(fc),'Found in PPI:',len(list(set([x for x in genenumber if x>0]))))
    return (genenumber,fc)
def lcc_zscore(genelist,fclist,ran,freq,maxfc):
    global G2
    if ran<=100:
        print("Preprocessing...\n")
    import random
    zscore=[]
    fccutoff=[]
    sizeoflcc=[]
    alfa=(maxfc-min(fclist))/freq
    for u in range(freq):
        source=[]
        for x in range(len(genelist)):
        #if((min(fc)+fcalfa*u)<=fc[x]<=(min(fc)+fcalfa*(u+1))):
            cutoff=min(fclist)+alfa*(u+1)
            if(fclist[x]>=cutoff):
                if(int(genelist[x])!=0):
                    source.append(genelist[x])
                
        if len(source)<=0:
            print("no nodes\n")
            continue
        genes=set()
        for j in source:
            genes.add(int(j))

################################################## Get LCC
        g=nx.subgraph(G2,genes)
        l=nx.connected_components(g)

        length=[]
        lcc=[]
        for i in l:
            length.append(len(i))
            lcc.append(i)
        
        list.sort(length,reverse=True)
        largest=length[0]
###LCC size
        if ran>100:
            print('size of LCC:',largest)
        sizeoflcc.append(largest)

        LC=[]
        for j in lcc:
            if len(j)==largest:
              LC.append(j)    
###LCC nodes
        #print(LC)
        p = {}
        for i in range(len(length)):
            size = length[i]
    #if not p.has_key(size):
            if not size in p:
                p[size] = 0
            p[size] += 1
   

####################################################zscore
        number_of_sims = ran
        all_genes = G2.nodes()
        all_p = {}
        all_largest = {}
        l_list = []

# simulations with randomly distributed seed nodes
        for i in range(number_of_sims):
            black_nodes = random.sample(all_genes,len(genes))
    
            g2=nx.subgraph(G2,black_nodes)
            l2=nx.connected_components(g2)
            length2=[]
            for i in l2:
                length2.append(len (i))
            list.sort(length2,reverse=True)
            largest2=length2[0]
    
    #iterkeys()
            for s in p.keys():
        #if not all_p.has_key(s):
                if not s in all_p:
                    all_p[s] = 0
                all_p[s] += p[s]

    #if not all_largest.has_key(largest2):
            if not largest2 in all_largest:
                all_largest[largest2] = 0
            all_largest[largest2] += 1
            l_list.append(largest2)

        l_mean = np.mean(l_list)
        l_std  = np.std(l_list)

        if l_std == 0:
            z_score = 0
        else:
            z_score = (1.*largest - l_mean)/l_std

###Z score
        if ran>100:
            print('z_score:',z_score)
            print('fc_cutoff:',pow(2,cutoff))
        zscore.append(z_score)
        fccutoff.append(cutoff)
    return (zscore,fccutoff,sizeoflcc)

def refine(zlist,fclist,size):
    flag=0
    max_fcvalue=max(fclist)
    for i in range(len(zlist)):
        #if zlist[i]>0.5 and size[i]>3:
        if size[i]>5:
            flag=0
        else:
            flag=flag+1
        if flag>=3:
            max_fcvalue=fclist[i]
            break
    return max_fcvalue

def leastsqfit(fclist,zlist):
    from scipy.optimize import leastsq
    z=[]
    for i in range(len(zlist)):
        z.append(pow(zlist[i],3))
    z1 = np.polyfit(fclist,z, 8)
    
    p1 = np.poly1d(z1)
    z2 = np.polyfit(fclist,zlist, 8)
    p2 = np.poly1d(z2)
    #print(z1)
    #print(p1)
    f=[]
    maximum=[]
    maximumfclist=[]
    for i in range(len(fclist)):
        f.append(p1(fclist[i]))
    #print(f)
    if f[0]>f[1]:
        maximum.append(f[0])
        maximumfclist.append(fclist[0])
    for i in range(1,len(fclist)-1):        
        if f[i]>f[i-1] and f[i]>f[i+1]:
            maximum.append(f[i])
            maximumfclist.append(fclist[i])
    
    #fccutoff1=
    
    s=maximum.index(max(maximum))
    #print(s)
    fccutoff1=maximumfclist[s]
    del maximum[s]
    del maximumfclist[s]
    if len(maximum)<1:
        fccutoff2=fccutoff1
    else:
        s=maximum.index(max(maximum))
        fccutoff2=maximumfclist[s]
        
    if fccutoff1>fccutoff2:
        t=fccutoff1
        fccutoff1=fccutoff2
        fccutoff2=t
    #print(fccutoff1,fccutoff2)    
    p=[]    
    s=fclist.index(fccutoff1)
    r=s+5
    if r>=50:
        r=50
    if s>=5:
        for i in range(s-5,r):
            p.append(zlist[i])
    if s<5:
        for i in range(10):
            p.append(zlist[i])
    s=zlist.index(max(p))
    fccutoff1=fclist[s]
    
    p=[]    
    s=fclist.index(fccutoff2)
    r=s+5
    if r>=50:
        r=50
    if s>=5:
        for i in range(s-5,r):
            p.append(zlist[i])
    if s<5:
        for i in range(10):
            p.append(zlist[i])
    s=zlist.index(max(p))
    fccutoff2=fclist[s]
        
        
    minfc=[]
    minzscore=[]
    #minfc.append()
    for i in range(len(fclist)):
        if (fccutoff1<=fclist[i])&(fclist[i]<=fccutoff2):
            minfc.append(fclist[i])
            minzscore.append(zlist[i])
    e=minzscore.index(min(minzscore))
    minpick=minfc[e]
    sig=abs(p2(fccutoff2)-p2(minpick))/p2(fccutoff1)
    #print(p1(fccutoff1),p1(fccutoff2),p1(minpick))
    #print(p2(fccutoff1),p2(fccutoff2),p2(minpick))
    return(p2,fccutoff1,fccutoff2,minpick,sig)
    #return(p2,maximumfclist)
    
    
def DrawPicture(titlename,xlim,ylim,p1,zlist,fclist,fc):
    
    size=11
    
    ### %pylab inline
    matplotlib.rc('xtick', labelsize=size-1) 
    matplotlib.rc('ytick', labelsize=size-1) 
    #%pylab inline
    newfc=[]
    newfclist=[]
    for i in range(len(fc)):
        newfc.append(pow(2,abs(fc[i])))
        #newfc.append(abs(fc[i]))
        
        #print(fc[i])
    for i in range(len(fclist)):
        newfclist.append(pow(2,fclist[i]))
        #newfclist.append(fclist[i])
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(4, 4),dpi=600)
    #axes.boxplot(newfc,vert=0,positions=[-1],widths=1)

    #axes[0][0].plot(fclistall[0],zlist1)
    #axes[0][0].plot(fclistall[0],zlist1,'ro')
    #axes[0].hist(newfc, bins=100, normed=True, color='seagreen')
    axes.plot(newfclist,zlist,'blue')
    axes.plot(newfclist,zlist,'ro')
    
    results, edges = np.histogram(newfc,bins=100, normed=True)
    binWidth = edges[1] - edges[0]
    axes.bar(edges[:-1], 16*results*binWidth, binWidth, color='seagreen')
    #weights = np.ones_like(newfc)/float(len(newfc))
    #axes[0].hist(newfc, weights=weights, bins=100, histtype='stepfilled', normed=True, color='seagreen')
    
    axes.plot(newfclist,p1(fclist),'gold')
    
    axes.set_title(titlename, fontsize=size)

    # adding horizontal grid lines

    axes.set_xlabel('Fold Change Cutoff', fontsize=size)
    axes.set_ylabel('Z-score', fontsize=size)
    axes.set_ylim([-2, ylim])
    axes.set_xlim([1, xlim])
    #axes[0].set_yticks(fontsize=size)

   
#plt.yaxis()
#plt.grid(True)
# add x-tick labels
#plt.axis([0, 10])
    plt.show()
    #fig.savefig(r'titlename+'.png',dpi=600, bbox_inches='tight')   
                
def getLCC(G,genes):
    g=nx.subgraph(G,genes)
    l=nx.connected_components(g)
    #print(l)
    length=[]
    lcc=[]
    flag=0
    for i in l:
        length.append(len(i))
        lcc.append(i)
        #print(len(i))
    #if flag==0:
        #return [[]]
    list.sort(length,reverse=True)
    largest=length[0]
###LCC size
    #print(largest)


    LC=[]
    for j in lcc:
        if len(j)==largest:
            for x in j:
                LC.append(x)    
###LCC nodes
    #print(LC)
    return LC

def getlist(cutoff1,cutoff2,genenumber,fc):
    module=[]
    coremodule=[]
   
    for x in range(len(genenumber)):
        if(fc[x]>=cutoff1):
    #if(fc[x]>=2.046):
            if(int(genenumber[x])!=0):
                module.append(genenumber[x])
        if(fc[x]>=cutoff2):
    #if(fc[x]>=2.046):
            if(int(genenumber[x])!=0):
                coremodule.append(genenumber[x])       
    return (module,coremodule)
