
module Verbose_print
    abstract type Verbose_level end

    struct Verbose_3  <: Verbose_level
        fp::Union{Nothing,IOStream}
        Verbose_3() = new(nothing)
        Verbose_3(filename::String) = new(open(filename,"w"))
        Verbose_3(fp::IOStream) = new(fp)
    end

    struct Verbose_2  <: Verbose_level
        fp::Union{Nothing,IOStream}
        Verbose_2() = new(nothing)
        Verbose_2(filename::String) = new(open(filename,"w"))
        Verbose_2(fp::IOStream) = new(fp)
    end

    struct Verbose_1  <: Verbose_level
        fp::Union{Nothing,IOStream}
        Verbose_1() = new(nothing)
        Verbose_1(filename::String) = new(open(filename,"w"))
        Verbose_1(fp::IOStream) = new(fp)
    end

    function Base.flush(v::Verbose_level)
        if v.fp != nothing
            flush(v.fp)
        end
    end

    function Base.versioninfo(v::Verbose_level)
        versioninfo(verbose=true)
        if v.fp != nothing
            versioninfo(v.fp,verbose=true)
        end
    end



    function println_verbose1(v::Verbose_3,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose2(v::Verbose_3,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose3(v::Verbose_3,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose1(v::Verbose_2,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose2(v::Verbose_2,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose3(v::Verbose_2,val...)
        return 
    end

    function println_verbose1(v::Verbose_1,val...)
        println(val...)
        if v.fp != nothing
            println(v.fp,val...)
        end
    end

    function println_verbose2(v::Verbose_1,val...)
        return
    end

    function println_verbose3(v::Verbose_1,val...)
        return 
    end

    function print_verbose1(v::Verbose_3,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose2(v::Verbose_3,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose3(v::Verbose_3,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose1(v::Verbose_2,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose2(v::Verbose_2,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose3(v::Verbose_2,val...)
        return 
    end

    function print_verbose1(v::Verbose_1,val...)
        print(val...)
        if v.fp != nothing
            print(v.fp,val...)
        end
    end

    function print_verbose2(v::Verbose_1,val...)
        return
    end

    function print_verbose3(v::Verbose_1,val...)
        return 
    end
end
