# %% markdown
# # Solution of the heat equation with the FDM
#
# ## $u_t=c^2\nabla^2u:0<x<a,\ 0<y<b$
#
# ---
# Let's import the required packages first of all
# %% import
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# %% visualize
def visualize(type, u, name, inv=200):
    if (type is steady):
        f=plot(0, u)
        f.savefig(name)
    else:
        animate(plot, u, u.shape[2], name, inv)
# %% markdown
# ---
# ### Before we get on with the problem at hand, let's make our life easier by defining function for frequent tasks
#
# ---
# The function to plot thermal map contours- takes the time co-ordinate in the data set and the data set(expected a 3d array-like structure)
# %% plotter
def plot(t, u):
    plt.clf()
    plt.contourf(u[:,:,t], cmap=plt.cm.jet)
    plt.colorbar()
    return plt
# %% markdown
# ---
# The function to animate our plots- takes the function responsible for creating each frame, the data set(expected a 3d array-like structure), number of frames to add in the movie, name to save the movie file with, and the timespan of each interval(optional, default value is 200 ms)
# %% animator
def animate(plot, u, f, name, inv):
    anim=ani.FuncAnimation(plt.figure(), plot, frames=f, fargs=(u,), interval=inv)
    anim.save(name)
# %% markdown
# ---
# The toplevel function to be accessed by the main code
# %% calculate
def calculate(bound, type, u, c=0.005, resolution=[1,1,1], accuracy=1):
    if(type is steady):
        std(bound, type, u, c, resolution, accuracy)
    if(type is transient):
        trans(bound, type, u, c, resolution, accuracy)
# %% markdown
# ---
# The function to manage the steady state computation
# %% steady
def std(bound, type, u, c_2, res, acc):
    k=1
    while(k):
        grid_fill(bound, type, 0, u, c_2, res)
        if(np.all(abs(u[:,:,1]-u[:,:,0])<=acc)):
            k=0
        u[:,:,0]=u[:,:,1]
# %% markdown
# ---
# The function to manage transient state computation
# %% transient
def trans(bound, type, u, c_2, res, acc):
    for t in range(u.shape[2]-1):
        grid_fill(bound, type, t, u, c_2, res)
# %% markdown
# ---
# A function to loop through our data and do the calculation - takes the type of equation to calculate (one of the two strings "steady" or "transient"), the data-set, the value of $c^2$ as in the original differential equation at the top of this documentation, and the resolution list, explained more later.
#
# _The steady state area is under construction and may be full of rubbish_
# %% edge
def edge(x, y):
    if (y== u.shape[0]-1 and x != 0):
        return 0
    elif (x== u.shape[1]-1 and y!=u.shape[0]-1):
        return 2
    elif (y==0 and x!=u.shape[1]-1):
        return 4
    elif (x==0 and y != 0):
        return 6
    else:
        return -1
# %% corner
def corner(x, y):
    if (y== u.shape[0]-1 and x == 0):
        return 3
    elif (x== u.shape[1]-1 and y==u.shape[0]-1):
        return 0
    elif (y==0 and x==u.shape[1]-1):
        return 1
    elif (x==0 and y == 0):
        return 2
    else:
        return -1
# %% grid fill
def grid_fill(bound, type, t, u, c_2, res):
        for y in range(0, u.shape[0]):
            for x in range(0, u.shape[1]):
                ev=edge(x,y)
                if(ev!=-1):
                    bound[ev](bound[ev+1], type, u, x, y, t, c_2, res)
                else:
                    type(u, x, y, t, c_2, res, xp=u[y][x+1][t], xm=u[y][x-1][t], yp=u[y+1][x][t], ym=u[y-1][x][t])
# %% markdown
# ---
# A function to represent the eqn- takes the type of the equation to acess(same as the calculate() function defined above), the dataset, the value of the grid point in the dataset to calculate for in the sequence, x, y, t, then the values of $c^2$, and the resolution list, explained in detail later.
#
# The transient state equation is:
#
# $u(x,y,t+\Delta t)=u(x,y,t)+\dfrac{\Delta t.c^2}{(\Delta x\Delta y)^2}((\Delta y)^2(u(x+\Delta x, y, t)+u(x-\Delta x,y,t)-2u(x,y,t))+(\Delta x)^2(u(x,y+\Delta y,t)+u(x, y-\Delta y, t)-2u(x,y,t)))$
#
# _The steady state area is under construction and may be full of rubbish_
# %% equation
def steady(u, x, y, t, c_2, res, xp, xm, yp, ym):
        u[y][x][t+1]=(1/4)*(yp+ym+xp+xm)
def transient(u, x, y, t, c_2, res, xp, xm, yp, ym):
        u[y][x][t+1]=u[y][x][t]+(res[2]**-1*c_2*(res[1]*res[0])**2)*(res[0]**-2*(xp+xm-2*u[y][x][t])+res[1]**-2*(yp+ym-2*u[y][x][t]))
# %% markdown
# ---
#
# ---
# ### Now let's define our boundary conditions
#
# ---
# For Dirichlet's boundary condition we have the data values given at boundaries at all times.
#
# The function takes the data set, and a dictionary with strings denoting sides("top", "bottom", "left", "right") and the corresponding value of the function
# %% Dirichlet
def dirichlet(value, type, u, x, y, t, c_2, res):
    u[y, x,:]=value
# %% markdown
# ---
# For Neumann's boundary condition we have the value of first derivative at the boundaries.
#
# This works the same as the dirichlet() above except this time the sides are given as a list of strings, rather than a dictionary.
#
# _This area is broken_
# %% Neumann
def neumann(value, type, u, x, y, t, c_2, res):
    ev=edge(x,y)
    cv=corner(x,y)
    if (ev==0 and cv==-1):
        type(u, x, y, t, c_2, res, xp=u[y][x+1][t], xm=u[y][x-1][t], ym=u[y-1][x][t], yp=u[y-1][x][t]+2*res[0]**-1*value)
    elif(ev==0 and cv==0):
        type(u, x, y, t, c_2, res, yp=u[y-1][x][t]+2*res[0]**-1*value, xp=u[y][x-1][t]+2*res[1]**-1*value, xm=u[y][x-1][t], ym=u[y-1][x][t])
    elif (ev==1 and cv==-1):
        type(u, x, y, t, c_2, res, xp=u[y][x-1][t]+2*res[1]**-1*value, xm=u[y][x-1][t], yp=u[y+1][x][t], ym=u[y-1][x][t])
    elif (ev==1 and cv==1):
        type(u, x, y, t, c_2, res, xp=u[y][x-1][t]+2*res[1]**-1*value, ym=u[y+1][x][t]-2*res[0]**-1*value, xm=u[y][x-1][t], yp=u[y+1][x][t])
    elif (ev==2 and cv==-1):
        type(u, x, y, t, c_2, res, ym=u[y+1][x][t]-2*res[0]**-1*value, xp=u[y][x+1][t], xm=u[y][x-1][t], yp=u[y+1][x][t])
    elif (ev==2 and cv==2):
        type(u, x, y, t, c_2, res, ym=u[y+1][x][t]-2*res[0]**-1*value, xm=u[y][x+1][t]-2*res[1]**-1*value, xp=u[y][x+1][t], yp=u[y+1][x][t])
    elif (ev==3 and cv==-1):
        type(u, x, y, t, c_2, res, xm=u[y][x+1][t]-2*res[1]**-1*value, xp=u[y][x+1][t], yp=u[y+1][x][t], ym=u[y-1][x][t])
    elif (ev==3 and cv==3):
        type(u, x, y, t, c_2, res, xm=u[y][x+1][t]-2*res[1]**-1*value, yp=u[y-1][x][t]+2*res[0]**-1*value, xp=u[y][x+1][t], ym=u[y-1][x][t])
# %% markdown
# ---
#
# ---
# ### Let's put some data that represents the system
#
# in the (x, y, t) format
#
# we have a list sys_scale that measures the dimensions of the physical system
#
# we have list res that measures resolution at which the program is supposed to function - here, resolution is defined as the number of steps a unit of system dimension has to be broken through all calculations. for example for a sys_scale value of 5(say cm), and a resolution value of 10, will create a grid with 1cm is represented as 10 steps, so we will get 50 steps
#
# then we have the value of $c^2$ as it appears in the original differential equation at the top of this documentation
#
# then we have a list of boundary value of the function in the clockwise order, i.e. top, right, bottom, left
# %% data
sys_scale=[5,5,60]
res=[10,10,2]
c_2=0.005
acc=.01
bound=[neumann, 100, dirichlet, 400, neumann, 100, dirichlet, 800]
analysis=transient
name="analysis.gif"
inv=50
# %% markdown
# ---
#
# ---
# ### Let's create a mesh of our surface first. I'll do that using a 2-D Numpy array.
#
# To add the time variation into the picture, I'll add another dimension to the array.
#
# Let's start with all grid points having a value of 0 intitally.
# %% mesh
arr=np.zeros((sys_scale[1]*res[1]+1,sys_scale[0]*res[0]+1,sys_scale[2]*res[2]+1))
# %% markdown
# ---
#
# ---
# ### Let's try the steady state analysis now
#
# First let's go with Dirichlet's conditions on all four boundaries
#
# Let's set the boundary conditions
#
# Then we calculate the points
#
# _The steady state area is under construction and maybe full of rubbish_
# %% Dirichlet's steady
u=np.copy(arr)
calculate(bound, analysis, u, accuracy=acc, c=c_2, resolution=res)
# %% markdown
# ---
# Now let's try to plot it
# %% plot Dirichlet's steady state
visualize(analysis, u, name , inv)
