
module Verbose_print
    abstract type Verbose_level end

    struct Verbose_3  <: Verbose_level
    end

    struct Verbose_2  <: Verbose_level
    end

    struct Verbose_1  <: Verbose_level
    end



    function println_verbose1(v::Verbose_3,val...)
        println(val...)
    end

    function println_verbose2(v::Verbose_3,val...)
        println(val...)
    end

    function println_verbose3(v::Verbose_3,val...)
        println(val...)
    end

    function println_verbose1(v::Verbose_2,val...)
        println(val...)
    end

    function println_verbose2(v::Verbose_2,val...)
        println(val...)
    end

    function println_verbose3(v::Verbose_2,val...)
        return 
    end

    function println_verbose1(v::Verbose_1,val...)
        println(val...)
    end

    function println_verbose2(v::Verbose_1,val...)
        return
    end

    function println_verbose3(v::Verbose_1,val...)
        return 
    end
end
