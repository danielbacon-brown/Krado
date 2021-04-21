#!/bin/env python3

from pylab import *
from mpl_toolkits.mplot3d import *

from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
from matplotlib.patches import FancyArrowPatch

rc('text', usetex=True)
rc('font', size=20)

fig = figure(figsize=(11,6))
ax = fig.gca(projection='3d')

ax.set_axis_off()

def f(x,t):
    return t/2 * 0.55*(sin(2*x)+0.4*x**2-0.65)

c_plane   = colorConverter.to_rgba('b', alpha=0.15)

N = 50
y = linspace(-1,1,N)
t = linspace(0,2,N)
yy, tt = meshgrid(y, t)
zz = f(yy,tt)

ax.plot(0*ones(y.shape), y, f(y,0), '-g', linewidth=3)
ax.plot(2*ones(y.shape), y, f(y,2), '-g', linewidth=3)

yt = 0.7*y
zt = f(yt, t) + 0.2*t

ax.plot(t, yt, zt, '-r', linewidth=3, zorder = 1)

ax.plot([2,2,2], [-1,yt[-1],yt[-1]], [zt[-1],zt[-1],-1], 'k--')
ax.plot(2*ones(y.shape), yt, f(yt,2)+0.1*(y+1), 'g:', linewidth=2)
ax.plot((2,2),
        (yt[0], yt[-1]),
        (f(yt[0], 2), f(yt[-1], 2) + 0.1*(y[-1]+1)), 'og')
ax.plot((0,2,2),
        (-1,-1,zt[-1]),
        (0,yt[-1],-1), 'ok')

ax.text(0, -1.1, 0, r'$p(0)=0$', ha='right', va='center')
ax.text(2, -1.05, zt[-1], r'$p(T)$', ha='right', va='center')
ax.text(0, -1.0, 1, r'$p$', ha='right', va='bottom')
ax.text(0, 1, -1.1, r'$q$', ha='center', va='top')
ax.text(0, -1, -1.1, r'$t=0$', ha='right', va='top')
ax.text(2, -1, -1.1, r'$t=T$', ha='right', va='top')
ax.text(2, yt[-1]-0.05, -1.05, r'$q(T)=q^*$', ha='left', va='top')
ax.text(0, 0.5, 0.05, r'$\mathcal{M}(0)$', ha='center', va='bottom')
ax.text(2, 0.1, -0.8, r'$\mathcal{M}(T)$', ha='center', va='bottom')

arrowprops = dict(mutation_scale=20,
                  linewidth=2,
                  arrowstyle='-|>',
                  color='k')

# For arrows, see
# https://stackoverflow.com/questions/29188612/arrows-in-matplotlib-using-mplot3d
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

a = Arrow3D([0,2], [-1,-1], [-1,-1], **arrowprops)
ax.add_artist(a)
a = Arrow3D([0,0], [-1,-1], [-1,1], **arrowprops)
ax.add_artist(a)
a = Arrow3D([0,0], [-1,1], [-1,-1], **arrowprops)
ax.add_artist(a)

# For surface illumination, see
# http://physicalmodelingwithpython.blogspot.de/2015/08/illuminating-surface-plots.html

# Get lighting object for shading surface plots.
from matplotlib.colors import LightSource

# Get colormaps to use with lighting object.
from matplotlib import cm

# Create an instance of a LightSource and use it to illuminate the surface.
light = LightSource(70, -120)
white = ones((zz.shape[0], zz.shape[1], 3))
illuminated_surface = light.shade_rgb(white*(0,1,0), zz)

ax.plot_surface(tt, yy, zz,
                cstride=1, rstride=1,
                alpha=0.3, facecolors=illuminated_surface,
                linewidth=0,
                zorder=10)

verts = [array([(-1,-1), (-1,1), (1,1), (1,-1), (-1,-1)])]

ax.plot_surface(((0,0),(0,0)), ((-1,-1),(1,1)), ((-1,1),(-1,1)),
                color=c_plane)

ax.plot_surface(((2,2),(2,2)), ((-1,-1),(1,1)), ((-1,1),(-1,1)),
                color=c_plane)

ax.plot((0,2), (yt[0], yt[-1]), (zt[0], zt[-1]), 'or')

ax.set_xlim3d(0, 2)
ax.view_init(elev=18, azim=-54)

show()
