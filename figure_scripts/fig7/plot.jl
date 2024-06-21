using DelimitedFiles
using CairoMakie

function plot_m_H()
  
  fig = Figure(size = (400, 300))
  ax = Axis(fig[1, 1],
      xlabel = "T (K)",
      ylabel = "m"
  )
	
  data = readdlm("assets/M_H.txt", skipstart=1)
  sc1 = scatter!(ax, data[:, 1], data[:, 2], markersize = 10, label="M-T curve")
  sc1.color = :transparent
  sc1.strokewidth = 1
  sc1.strokecolor = :purple
  lines!(ax, data[:, 1], data[:, 2])
  
  vlines!(ax, 431, linewidth=1, color=:black, linestyle=:dash)
  
  axislegend()

  save("M_T.pdf", fig)

  return fig

end

plot_m_H()