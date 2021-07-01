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
def animate(plot, u, f, name, inv=200):
    anim=ani.FuncAnimation(plt.figure(), plot, frames=f, fargs=(u,), interval=inv)
    anim.save(name)
# %% markdown
# ---
# A function to loop through our data and do the calculation - takes the type of equation to calculate (one of the two strings "steady" or "transient"), the data-set, the value of $c^2$ as in the original differential equation at the top of this documentation, and the resolution list, explained more later.
#
# _The steady state area is under construction and may be full of rubbish_
# %% looper
def calculate(type, u, c_2, res):
    for t in range(0, u.shape[2]-1):
        for y in range(1, u.shape[0]-1):
            for x in range(1, u.shape[1]-1):
                eqn(type, u, x, y, t, c_2, res)
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
def eqn(type, u, x, y, t, c_2, res):
    if (type=="steady"):
        u[y][x][0]=(1/4)*(u[y+1][x][0]+u[y-1][x][0]+u[y][x+1][0]+u[y][x-1][0])
    if(type=="transient"):
        u[y][x][t+1]=u[y][x][t]+(res[2]**-1*c_2*(res[1]*res[0])**2)*(res[0]**-2*(u[y][x+1][t]+u[y][x-1][t]-2*u[y][x][t])+res[1]**-2*(u[y+1][x][t]+u[y-1][x][t]-2*u[y][x][t]))
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
def dirichlet(u, edge):
    if("top" in edge):
        u[u.shape[1]-1,:,:]=edge["top"]
    if("bottom" in edge):
        u[0,:,:]=edge["bottom"]
    if("right" in edge):
        u[:,u.shape[0]-1,:]=edge["right"]
    if("left" in edge):
        u[:,0,:]=edge["left"]
# %% markdown
# ---
# For Neumann's boundary condition we have the value of first derivative at the boundaries.
#
# This works the same as the dirichlet() above except this time the sides are given as a list of strings, rather than a dictionary.
#
# _This area is broken_
# %% Neumann
def neumann(u, edge):
    if("top" in edge):
        u[u.shape[1]-1,:,:]=0
    if("bottom" in edge):
        u[0,:,:]=0
    if("right" in edge):
        u[:,u.shape[0]-1,:]=0
    if("left" in edge):
        u[:,0,:]=0
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
res=[10,10,5]
c_2=0.005
bound=[100,10,100,500]
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
dirichlet(u, {"top":bound[0], "right":bound[1], "bottom":bound[2], "left":bound[3]})
calculate("steady", u)
# %% markdown
# ---
# Now let's try to plot it
# %% plot Dirichlet's steady state
plot(0)
plt.savefig("dirichlet_steady.svg")
# %% markdown
# ---
#
# ---
# Let's try the Nuemann's condition on two ends and Dirichlet's condition on two ends.
#
# Same process as above.
#
# _The steady state area is under construction and may be full of rubbish_
# %% Neumann's steady
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "bottom":bound[2]})
neumann(u, ["right", "left"])
calculate("steady", u)
# %% markdown
# ---
# Now let's try to plot it
# %% plot Neumann's steady state
f=plot(0)
f.savefig("neumann_steady.svg")
# %% markdown
# ---
#
# ---
# ### Let's try the transient state analysis
#
# First let's try the Dirichlet's boundary condition
# %% Dirichlet's transient
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "right":bound[1], "bottom":bound[2], "left":bound[3]})
calculate("transient", u, c_2, res)
# %% markdown
# ---
# Let's try to animate it
# %% Animate Dirichlet's transient state
animate(plot, u, u.shape[2], "dirichlet_trans.gif", 50)
# %% markdown
# ---
#
# ---
# Let's try the Neumann's boundary conditions
# %% Neumann's transient
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "bottom":bound[2]})
neumann(u, ["right", "left"])
calculate("transient", u, c_2, res)
# %% markdown
# ---
# Let's get it animated
# %% Animate Neumann's transient state
animate(plot, u, u.shape[2], "neumann_trans.gif",50)
# %% markdown
# ---
#
# ---
#
# ---
# below is the experimentation zone
#
# %% trial
