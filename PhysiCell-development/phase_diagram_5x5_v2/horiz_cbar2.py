import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# 1. Define the discrete values and corresponding colors
bounds = [0, 1, 2, 3,4] 
# colors = ['#2c7bb6', '#0a793a', '#77a353', 'red'] 
colors = ['darkblue', 'lightblue', 'yellow', 'red'] 
                    # Lutz: light-green(uninhibited), light-blue, yellow, red
                    # cbar_name = from_list(None, [[0.5, 1, 0.5],[0,0.5,1],[1,1,0],[1,0,0]], len(self.discrete_variable))
                    # Roman's cartography app: dark blue, light-blue, yellow, red
                    # Red: 215,25,28
                    # yellow: 253,174,97
                    # light blue: 171,217,233
                    # blue: 44,123,182

                    # cbar_name = from_list(None, [[44./255, 123./255, 182./255],[171./255,217./255,233./255],[253./255,174./255,97./255],[215./255,25./255,28./255]], len(self.discrete_variable))
                    # cbar_name = from_list(None, [[44/255.,123/255.,182/255.], [253/255.,174/255.,97/255.], 
                        # [171/255.,217/255.,233/255.],[215/255.,25/255.,28/255.]], len(self.discrete_variable))

cmap = mpl.colors.ListedColormap(colors)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# 2. Create a figure and an axes for the colorbar
# fig.subplots_adjust is used to make space for the colorbar and its label
fig, ax = plt.subplots(figsize=(5, 2))
fig.subplots_adjust(bottom=0.5)

# 3. Create a ScalarMappable object to pass to colorbar (since there's no main plot)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([]) # Set an empty array

# 4. Add the colorbar to the axes
cbar = fig.colorbar(sm, cax=ax, orientation='horizontal', ticks=bounds)

# 5. Customize the colorbar
cbar.set_label('Contact Inhibition', loc='center', labelpad=15)
cbar.ax.xaxis.set_ticks_position('bottom') 

# To center the tick labels within the color segments, an alternative approach to setting ticks is needed.
# By default, ticks are at the boundaries. You can use the midpoint of bounds for centering labels.
midpoints = [(b1 + b2) / 2 for b1, b2 in zip(bounds[:-1], bounds[1:])]
cbar.set_ticks(midpoints)
cbar.set_ticklabels(['None','Pressure', 'Surface', 'Both'])

plt.show()

