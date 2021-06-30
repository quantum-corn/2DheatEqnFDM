# %% markdown
# # Solution of the heat equation with the FDM
#
# $u_t=c^2\nabla^2u:0<x<a,\ 0<y<b$
#
# Let's import the required packages first of all
# %% import
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# %% markdown
# Before we get on with the problem at hand, let's make our life easier by defining function for frequent tasks
# the function to plot heatmap contours
# %% plotter
def plot(t, u):
    plt.clf()
    plt.contourf(u[:,:,t], cmap=plt.cm.jet)
    plt.colorbar()
    return plt
# %% markdown
# the function to animate our plots
# %% animator
def animate(plot, u, f, name, inv):
    anim=ani.FuncAnimation(plt.figure(), plot, frames=f, fargs=(u,), interval=inv)
    anim.save(name)
# %% markdown
# a function loop through our array and do the calculation
# %% looper
def calculate(type, u, c_2=None, res=None):
    for t in range(0, u.shape[2]-1):
        for x in range(1, u.shape[0]-1):
            for y in range(1, u.shape[1]-1):
                eqn(type, u, x, y, t, c_2, res)
# %% markdown
# a function to represent the eqn
# %% equation
def eqn(type, u, x, y, t, c_2, res):
    if (type=="steady"):
        u[x][y][0]=(1/4)*(u[x+1][y][0]+u[x-1][y][0]+u[x][y+1][0]+u[x][y-1][0])
    if(type=="transient"):
        u[x][y][t+1]=u[x][y][t]+(res[2]**-1*c_2*(res[0]res[1])**2)*(res[1]**-2*(u[x+1][y][t]+u[x-1][y][t]-2*u[x][y][t])+res[0]**-2*(u[x][y+1][t]+u[x][y-1][t]-2*u[x][y][t]))
# %% markdown
# Now let's define our boundary conditions
# For Dirichlet's boundary condition we have the value of boundaries.
# %% Dirichlet
def dirichlet(u, edge):
    if("top" in edge):
        u[:,u.shape[1]-1,:]=edge["top"]
    if("bottom" in edge):
        u[:,0,:]=edge["bottom"]
    if("right" in edge):
        u[u.shape[0]-1,:,:]=edge["right"]
    if("left" in edge):
        u[0,:,:]=edge["left"]
# %% markdown
# For Neumann's boundary condition we have the value of first derivative at the boundaries.
# %% Neumann
def neumann(u, edge):
    if("top" in edge):
        u[:,u.shape[1]-1,:]=0
    if("bottom" in edge):
        u[:,0,:]=0
    if("right" in edge):
        u[u.shape[0]-1,:,:]=0
    if("left" in edge):
        u[0,:,:]=0
# %% markdown
# Let's put some data that represents the system
# in the (y, x, t) format
# %% data
sys_scale=[1,1,10]
res=[25,25,10]
c_2=1
bound=[100,100,100,100]
# %% markdown
# Let's create a mesh of our surface first. I'll do that using a 2-D Numpy array.
# To add the time variation into the picture, I'll add another dimension to the array.
# %% mesh
arr=np.zeros((sys_scale[0]*res[0]+1,sys_scale[1]*res[1]+1,sys_scale[2]*res[2]+1))
# %% markdown
# ### Let's try the steady state analysis now
# First let's go with Dirichlet's conditions on all four boundaries
# Let's set the boundary conditions
# Then we calculate the points
# %% Dirichlet's steady
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "right":bound[1], "bottom":bound[2], "left":bound[3]})
calculate("steady", u)
# %% markdown
# Now let's try to plot it
# %% plot Dirichlet's steady state
plot(0)
plt.savefig("dirichlet_steady.svg")
# %% markdown
# Let's try the Nuemann's conditions on two ends and D-irichlet's conditions on two ends.
# Same process as above.
# %% Neumann's steady
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "bottom":bound[2]})
neumann(u, ["right", "left"])
calculate("steady", u)
# %% markdown
# Now let's try to plot it
# %% plot Neumann's steady state
f=plot(0)
f.savefig("neumann_steady.svg")
# %% markdown
# ### Let's try the transient state analysis
# First let's try the Dirichlet's boundary condition
# %% Dirichlet's transient
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "right":bound[1], "bottom":bound[2], "left":bound[3]})
calculate("transient", u, c_2, res)
# %% markdown
# Let's try to animate it
# %% Animate Dirichlet's transient state
animate(plot, u, u.shape[2], "dirichlet_trans.gif", 100)
# %% markdown
# Let's try the Neumann's boundary conditions
# %% Neumann's transient
u=np.copy(arr)
dirichlet(u, {"top":bound[0], "bottom":bound[2]})
neumann(u, ["right", "left"])
calculate("transient", u, c_2, res)
# %% markdown
# Let's get it animated
# %% Animate Neumann's transient state
animate(plot, u, u.shape[2], "neumann_trans.gif")
# %% markdown
# below is the experimentation zone
# %% trial
