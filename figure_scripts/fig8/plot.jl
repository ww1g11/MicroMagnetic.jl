using DelimitedFiles
using CairoMakie

function plot()

    data1 = readdlm("assets/cuda_time.txt", skipstart=1)
    data2 = readdlm("assets/benchmark.txt", skipstart=0)

    fig = Figure(size = (500, 360), fontsize = 18)
    ax = Axis(fig[1, 1],
        xlabel = "N",
        ylabel = "Throughput (M cells/s)",
        #xscale = log10,
        #yscale = log10,
        limits = ((nothing), (0,1200))
    )

    scatterlines!(ax, sqrt.(data2[:,1]), data2[:,3]/1e6, markersize = 12, color=:slateblue1, strokecolor=:transparent, label="MuMax3")
    scatterlines!(ax, sqrt.(data1[:,1]), data1[:,3]/1e6, markersize = 12, marker = :rect, color = :sienna1, label="MicroMagnetic.jl")
    
    axislegend(position=(0.95, 0.05), labelsize=16)

    save("throughput.pdf", fig)

end

plot()