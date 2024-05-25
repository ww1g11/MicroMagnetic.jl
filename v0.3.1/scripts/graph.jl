
gvs = ["sim.gv",
        "mesh.gv",
        "energy.gv",
        "driver.gv"]

for gv in gvs
    png = gv[1:end-2]*"png"
    cmd = `dot -Tpng -Gdpi=200 $gv -o $png`
    println(cmd)
    run(cmd)
end
