#!/usr/bin/env python2.7

import math,operator

def mean(seq,w):
    return sum(map(operator.mul,seq,w))/float(sum(w))

def var(seq,w):
    avg=mean(seq,w)
    return sum(map(operator.mul,[(x-avg)**2 for x in seq],w))/float(sum(w))

def std(seq,w):
    return var(seq,w)**0.5
    
def stderr(seq,w):
    if len(seq)>1:
        return (var(seq,w)/(len(seq)-1))**0.5
    else:
        return 0.0

class Stat:
    def __init__(self,*args,**kargs):
        self.dict=dict((x,()) for x in args)
        self.dict.update(kargs)
        for k,v in self.dict.items():
            assert not isinstance(v, (str, unicode));
            if isinstance(v,int):
                v=(v,)
                self.dict[k]=v
        self.operations={}
    def clear(self):
        for k,v in self.dict.items():
            if len(v)==0:
                self.__dict__[k]=[]
            elif len(v)==1:
                self.__dict__[k]=[[] for x in range(v[0])]
            elif len(v)==2:
                self.__dict__[k]=[[[] for y in range(v[1])] for x in range(v[0])]
            elif len(v)==3:
                self.__dict__[k]=[[[[] for z in range(v[2])] for y in range(v[1])] for x in range(v[0])]
    def register(self,prefix,operation):
        self.operations[prefix]=operation
        for k in self.dict:
            if prefix+k in self.__dict__:
                raise IndexError
    def update(self,weight):
        for k,v in self.dict.items():
            a=getattr(self,k)
            for pf,op in self.operations.items():
                avgName=pf+k
                #print avgName,':',a,';',v
                if len(v)==0:
                    setattr(self,avgName,op(a,weight))
                elif len(v)==1:
                    setattr(self,avgName,[op(a[x],weight) for x in range(v[0])])
                elif len(v)==2:
                    setattr(self,avgName,[[op(a[x][y],weight) for y in range(v[1])] for x in range(v[0])])
                elif len(v)==3:
                    setattr(self,avgName,[[[op(a[x][y][z],weight) for z in range(v[2])] for y in range(v[1])]])

if __name__=='__main__':
    a=Stat('Apple',Pearl=(2),Banana=(2,3))
    a.register('Avg',mean)
    a.register('SD',std)
    a.register('SE',stderr)
    a.Apple=[1,2,3,4,5]
    for i in range(2):
        a.Pearl[i]=[1,2,3,4,5]
    for i in range(2):
        for j in range(3):
            a.Banana[i][j]=[1,2,3,4,5]
    weight=[1.0 for x in a.Apple]
    a.update(weight)
    print 'AvgApple,SDApple,SEApple:'
    print a.AvgApple,a.SDApple,a.SEApple
    print 'AvgPearl,SDPearl,SEPearl:'
    for i in range(2):
        print a.AvgPearl[i],a.SDPearl[i],a.SEPearl[i]
    print 'AvgBanana,SDBanana,SEBanana:'
    for i in range(2):
        for j in range(3):
            print a.AvgBanana[i][j],a.SDBanana[i][j],a.SEBanana[i][j]
