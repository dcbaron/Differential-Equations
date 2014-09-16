x,y,z,t='x','y','z','t'
from sys import float_info; eps=float_info.epsilon
from random import choice
import matplotlib.pyplot as plt

def sign(x):
    if x>0:return 1
    elif x<0: return -1
    else: return 0

def deriv(expr):
    return lambda **args: eval(expr, dict(args,sign=sign))#,sign=sign)

def getdecimalplaces(x):
    n=0
    while round(x,n) != x and n<100:
        n+=1
    return n

reverse = lambda l : type(l)(reversed(l))


def eulers(derivs, vars=[t,x,y], init_conds=[0,1,1],
           interval=[0,1], step=.01, bound=8, include_t=False):
    '''derivs:[x', y',...], dep_vars=[x, y, ...]'''
    
    n=getdecimalplaces(step)
    
    if len(vars)==len(derivs): vars=['dep_var']+list(vars)
    tname=vars[0]
    

    derivs=[deriv(d) for d in derivs] 
    dtdt = lambda **args : 1
    derivs= dict(zip(vars, [dtdt]+derivs))

    init_conds=dict(zip(vars, init_conds))
    t=init_conds[tname]
    
    def getNewVals(vals, step):
        newvals={}
        for var in vals:
            loop=True
            while loop:
                try:
                    val = vals[var]+step*derivs[var](**vals)
                    loop=False
                except ZeroDivisionError:
                    vals[choice(vals.keys())]+=eps
            newvals[var]=val
        newvals[tname]=round(vals[tname]+step,n)
        return newvals, newvals[tname]

    def addpoint(vals, points, t):
        return points + [(t,tuple(vals[var] for var in vars[1:]))]
        
    def get_points(a,b,t,init_conds,step):
        vals = init_conds.copy()
        
        points=[]
        stepsign=sign(step)
        while sign(b-t)==stepsign:
            oldvals = vals
            
            if (sign(t-a) in {0,stepsign}) \
               or (abs(a-t) >= abs(step)):
                vals, t = getNewVals(vals,step)
            else:
                vals, t = getNewVals(vals,a-t)
                t=vals[tname]=a

            valset=set()
            for v in vars[1:]: valset.add(abs(vals[v]))
            if max(valset) > bound: return points

            if sign(t-a) in {0,stepsign}:
                points = addpoint(vals, points, t)

        points.pop()        
        t=oldvals[tname]
        vals, t = getNewVals(oldvals, b-t)
        points = addpoint(vals, points, b)
        return points
    
    if len(interval)!=2:raise
    a,b=interval
    if a>=b: raise

    if (a<=t) and (b>=t):
        highpoints=addpoint(init_conds, [], t)
    else: highpoints=[]
    lowpoints=[]
    
    if b>t:
        highpoints += get_points(a,b,t,init_conds,step)        
    
    if t!=init_conds[tname]:raise
    
    if a<t:
        lowpoints = get_points(b,a,t,init_conds,-step)
        
    points=reverse(lowpoints)+highpoints
    return points if include_t else [x[1] for x in points]


p=eulers([y,'-x**3-y**3'],[t,x,y],[0,1,1],[-1,1],.01)

def tups_to_lists(p):
    return [q[0] for q in p],[q[1] for q in p]

def eulerslists(derivs, vars=[t,x,y], init_conds=[0,1,1],
           interval=[0,1], step=.01, bound=8, include_t=False):
    return tups_to_lists(eulers(
        derivs, vars, init_conds,
           interval, step, bound, include_t))

def combineplots(*tups):
    '''takes arbitrarily many 2-tuples (X-list, Y-list);
    returns a 2-tuple'''
    def loop(tups=tups,X=[],Y=[]):
        empty = True
        x,y=[],[]
        for tup in tups:
            x.append(tup[0].pop(0))
            y.append(tup[1].pop(0))
            if tup[0]==[]:
                tup[0].append(x[-1])
                tup[1].append(y[-1])
            else:
                empty=False
        if empty:
            return X+[x], Y+[y]
        else:
            return loop(tups, X+[x], Y+[y])
    return loop()

def lotsapaths(xmin=-1,xmax=1,ymin=-1,ymax=1,spacing=.30,**args):
    tups = []
    x,y=xmin,ymin
    while x<=xmax:
        while y<=ymax:
            args['init_conds']=[0,x,y]
            tups.append(eulerslists(**args))
            y+=spacing#-spacing*(1-abs(y))*.9
        y=ymin
        x+=spacing#-spacing*(1-abs(x))*.9
    return combineplots(*tups)

def makecoords(expr, BL,UR,gridspacing=.1):
    from pylab import arange, array
    X=arange(BL[0],UR[0],gridspacing)
    Y=arange(BL[1],UR[1],gridspacing)
    fn=deriv(expr)
    Z=array([[fn(x=x,y=y) for x in X] for y in Y])
    return X,Y,Z

#X,Y,UV =makecoords('[(-x-y**2),(-y-x**2)]',[-8,-8],[10,10],gridspacing=1)
#U=[[x[0] for x in row] for row in UV]
#V=[[x[1] for x in row] for row in UV]
#plt.quiver(X,Y,U,V,units='inches',scale_units='width')


def phaseport(derivs,BL,UR,gridspacing=.01,linedensity=1):
    from pylab import arange, array
    X=arange(BL[0],UR[0],gridspacing)
    Y=arange(BL[1],UR[1],gridspacing)
    dx, dy = map(deriv,derivs)
    U=array([[dx(x=x,y=y) for x in X] for y in Y])
    V=array([[dy(x=x,y=y) for x in X] for y in Y])
    return plt.streamplot(X,Y,U,V,density=linedensity)


def thing(n):
    if n==1:
        p=lotsapaths(xmin=-3,
             ymin=-3,
             xmax=3,
             ymax=3,
             spacing=.25,
             derivs=[y,'-x**3-y**3'],#'-x-y**2','-y-x**2'],
             vars=[t,x,y],
             interval=[-4,4],
             step=.01,
             bound=8)
        plt.plot(*p)
    elif n==2:
        p=phaseport(['-x-y**2','-y-x**2'],#y,'-x**3-y**3'],
            [-2,-2],
            [2,2],
            linedensity=1)
    if n==3:
        p=plt.contour(*makecoords('x**4+2*y**2',[-5,-3],[5,3]),levels=range(10))
    #plt.show()
    return p

def make(init_conds):
    p=plt.plot( *tups_to_lists( eulers(
        [y,'-x**3-y**3'],
        init_conds=init_conds,
        interval=[0,100],
        step=.01)))
    return p

def transform(x,y,tup):
    '''tup= (U-list, V-list),
    U-list = [(x1(t0),x2(t0),...),(x1(t1),x2(t1),...),...]
    x,y = ['f(u,v)','g(u,v)']'''
    U=tup[0]
    V=tup[1]
    UV=zip(U,V)
    x,y= deriv(x),deriv(y)

    UV=(zip(utup, vtup) for utup,vtup in UV)
    X=[[x(u=u,v=v) for u,v in tlist] for tlist in UV]
    Y=[[y(u=u,v=v) for u,v in tlist] for tlist in UV]
    return X,Y

##make([0,2,2])
##make([0,-2,2])
##thing(2)
##plt.annotate('x-nullcline',xy=[-1.6,0])
##plt.plot([-2,2],[0,0],'k',[-2,2],[2,-2],'k')
##plt.annotate('y-nullcline',xy=[-1.6,1])
##plt.plot([0,0],[-2,2],'--')
##plt.xlabel("x-axis")
##plt.ylabel('y-axis')
##thing(2)


##system1=['x/(x**2+y**2)**.5 - x-y*sign(1-x**2-y**2)',
##        'y/(x**2+y**2)**.5 + x*sign(1-x**2-y**2)-y']
##system2=['x*(x**2+y**2) - x-y',
##        'y*(x**2+y**2) + x-y']
##system=['-x+y*(1/10 + x**2)','1/2-y*(1/10 + x**2)']
##p=lotsapaths(xmin=0,
##             ymin=0,
##             xmax=5,
##             ymax=5,
##             spacing=.25,
##             derivs=system,#'-x-y**2','-y-x**2'],
##             vars=[t,x,y],
##             interval=[-4,4],
##             step=.01,
##             bound=6)
##p=transform('u','v',p)
##
##plt.plot(*p)
##
##plt.show()

## Main Loop
if __name__ == '__main__':
    print "Enter a two-dimensional system of differential equations:"
    xprime = raw_input("x' = ")
    yprime = raw_input("y' = ")
    print "Enter a list of initial conditions x(0), y(0)."
    print "Separate each pair of values with a semicolon."
    print "E.g. 1,1 ; 1,2 ; 0,0 ; -1,3"
    print
    print "OR type 'skein' to view many trajectories at once."
    conds_input = raw_input("?: ")
    if 'skein' in conds_input.lower():
        pass
    else:
        conds_input = conds_input.split(';')
        init_conds = []
        for cond in conds_input:
            init_conds.append([0]+[int(x) for x in cond.split(',')])
        plt.plot(*eulerslists([xprime,yprime], init_conds=init_conds[0]))
        plt.show()
    

    
    
