struct GivenConfigurations <: AbstractUpdate
    loadU_dir::String
    loadU_format::String
    filename_load::Vector{String}
    current_fileindex::Vector{Int64}

    function GivenConfigurations(U,loadU_dir,loadU_format)
        ildg = nothing
        println_verbose_level1(U[1],"load U from ",loadU_dir)
        if loadU_format == "JLD"
            datatype = "JLD"
            filename_load =  filter(f -> contains(f,".jld2"),readdir("./$(loadU_dir)"))
            #filename_load =  filter(f -> contains(f,".jld"),readdir("./$(parameters.loadU_dir)"))
        elseif loadU_format == "ILDG"
            datatype = "ILDG"
            filename_load =  filter(f -> contains(f,"ildg"),readdir("./$(loadU_dir)"))
        elseif loadU_format == "BridgeText"
            datatype = "BridgeText"
            filename_load =  filter(f -> contains(f,"txt"),readdir("./$(loadU_dir)"))
        else
            error("loadU_format should be JLD, ILDG or BridgeText. now $(loadU_format)")
        end


        numfiles = length(filename_load)
        println_verbose_level1(U[1],"Num of files = $numfiles")
        for file in filename_load
            println_verbose_level1(U[1],"$file")
        end

        current_fileindex = Vector{Int64}(undef,1)
        current_fileindex[1] = 0

        return new(loadU_dir,loadU_format,filename_load,current_fileindex)

    end
end

function update!(updatemethod::T,U) where T <: GivenConfigurations
    updatemethod.current_fileindex[1] += 1
    itrj = updatemethod.current_fileindex[1]
    fileformat = updatemethod.loadU_format
    loadU_dir = updatemethod.loadU_dir
    filename_i = updatemethod.filename_load[itrj]
    _,_,L... = size(U[1])
    NC = U[1].NC

    if fileformat == "JLD"
        error("not supported now")
    elseif fileformat == "ILDG"
        ildg = ILDG(loadU_dir*"/"*filename_i)
        i = 1
        load_gaugefield!(U,i,ildg,L,NC)
    elseif fileformat == "BridgeText"
        filename = loadU_dir*"/"*filename_i
        load_BridgeText!(filename,U,L,NC)
    end
end

