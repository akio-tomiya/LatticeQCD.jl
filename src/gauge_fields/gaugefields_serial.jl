module Gaugefields_serial

    

    struct GaugeFields{NC,NX,NY,NZ,NT,NDW}
        g::Array{ComplexF64,6}
        NV::Int64

        function GaugeFields(NC,NDW,NX,NY,NZ,NT)
            NV = NX*NY*NZ*NT
            g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        end
    end
end