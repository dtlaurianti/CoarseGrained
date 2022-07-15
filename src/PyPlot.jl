using PyPlot
#function plot_smoothed_surface
#
#Input: A csv file which has been smoothed by the R function
#Output: No output, just a plot
function plot_smoothed_surface(data)
  #Grab data from file
  df = DataFrame(CSV.File(data))
  x = df.x
  y = df.y
  z = df.z

  # plot surface
  pyplot()
  pygui(true)
  plt = Plots.plot(x, y, z, st=:surface, extra_kwargs=Dict(:subplot=>Dict("3d_colorbar_axis" => [0.85, 0.05, 0.05, 0.9])))
  display(plt)
  return plt
end

#function plot_surface
#
#Input: x,y,z vectors
#Output: a 3d surface plot of those points
function plot_surface(x, y, z)
  pyplot()
  pygui(true)
  plt = Plots.plot(x, y, z, st=:surface, extra_kwargs=Dict(:subplot=>Dict("3d_colorbar_axis" => [0.85, 0.05, 0.05, 0.9])))
  display(plt)
  return plt
end




# function subplot_smoothed_surface
# a subplot version of plot_smoothed_surface
# Input: A csv file which has been smoothed by the R function, (optional) ax
# Output: A subplot
function subplot_smoothed_surface(data, fig, ax=None)
  #Grab data from file
  df = DataFrame(CSV.File(data))
  x = df.x
  y = df.y
  z = df.z

  #Plot surface
  triang = mtri.Triangulation(x,y)
  if ax != None
      ax = fig.add_subplot(1,2,ax,projection="3d")
  else
      ax = fig.add_subplot(1,2,1,projection="3d")
  end
  surf = ax.plot_trisurf(triang,z,cmap=plt.cm.CMRmap,antialiased=True)
  return surf
end

#Test Examples
#G = nx.fast_gnp_random_graph(10, 0.3)
#A = nx.to_numpy_array(G)
#P = RN.generateRandomPartitions(10,7,1000)

#surfaceplots(P,A,"n10_p7_s1000")

# plot_smoothed_surface("data/visualization_data/test.csv")
