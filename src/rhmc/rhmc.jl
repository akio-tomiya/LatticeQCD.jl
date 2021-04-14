module Rhmc
    import ..AlgRemez:AlgRemez_coeffs,calc_coefficients

    const coeffs_12 = AlgRemez_coeffs(35.76877880344718, [-2.611603580931762e-6, -2.8381078304162547e-5, -0.00024575851782184877, -0.0020470080999742007, -0.016923799105426088, -0.14034361960343106, -1.1877898203825932, -10.955212949978659, -146.07789284928893, -14714.341373866839], [0.00023220051347978252, 0.0014679843480651033, 0.006577084315568452, 0.027426964178005095, 0.11250019513122003, 0.46049080638938455, 1.8988299965694169, 8.104664281763002, 40.29762890763353, 497.70654517555414], 10)
    const coeffs_m12 = AlgRemez_coeffs(0.027957342505180133, [0.009504168034266293, 0.014392806670277538, 0.02668521749757119, 0.05270932006455869, 0.10589384745339038, 0.2139495082299973, 0.4353956294494188, 0.9089955204708765, 2.107201294600078, 7.749985438494519], [5.143593197266517e-5, 0.0006352731089632578, 0.003158674944452018, 0.013481986300116959, 0.055592857978478315, 0.22755516086118974, 0.9333880277034011, 3.892302237847687, 17.438878032822647, 110.24954086601954], 10)
    const coeffs_14 = AlgRemez_coeffs(6.9461071718499525, [-3.7385105377001425e-6, -1.7764883979837546e-5, -6.334438533952787e-5, -0.00021317868898123004, -0.0007061292823812292, -0.0023271447494456277, -0.007658442853216887, -0.025215367727471134, -0.08325492911813415, -0.27706279011926743, -0.9413396741768589, -3.376769295570277, -14.000699219260555, -87.48968373506541, -2438.159453743938], [7.533204879983615e-5, 0.0004155648337491737, 0.001355472748322587, 0.0038120901090773723, 0.010182414902576235, 0.02668386113490887, 0.06943317053167551, 0.18025949646037842, 0.46811809693900475, 1.219485561990294, 3.205950452609407, 8.636487104246942, 24.866247747939997, 87.21468047228943, 671.5034010759199], 15)
    const coeffs_m14 = AlgRemez_coeffs(0.14396552993778106, [0.00034605541183977854, 0.0007361352806089118, 0.0014491361385218552, 0.002897390997668544, 0.005861553385743302, 0.011923530079141716, 0.0243152683972108, 0.04966483224890658, 0.10166850937810712, 0.20917319384688998, 0.4358756590644266, 0.9388526693245018, 2.2065131459490663, 6.583620229613641, 42.16178846519866], [3.8123410780916776e-5, 0.0002935285649316097, 0.0010295079603283045, 0.002964168149734324, 0.007985151479544388, 0.020992458457826137, 0.054687054756901776, 0.14201748314340246, 0.3686998563362629, 0.959381398013239, 2.5141383694277666, 6.715476095132465, 18.886399620855748, 61.602902654298305, 339.8288033824946], 15)
    const coeffs_18 =AlgRemez_coeffs(2.6226563721874356, [-5.960771208125207e-6, -2.365926748832181e-5, -7.312974957271435e-5, -0.00021647258378734637, -0.0006342614065525487, -0.0018530026882076763, -0.005410134459579222, -0.01580639454625644, -0.04630025630383474, -0.13658135694721724, -0.4103890831099445, -1.2939417075231405, -4.636609634886856, -23.781760392176484, -426.7735702750801], [6.511406585080085e-5, 0.0003827736362029256, 0.001267715019499123, 0.003583090268588116, 0.009587354346137916, 0.02513865988838697, 0.06541997845298224, 0.169826364132841, 0.4409229237078518, 1.1481048390885387, 3.015181840776192, 8.102022022873221, 23.168853584928378, 79.54719981908272, 550.2916347706022], 15)
    const coeffs_m18 = AlgRemez_coeffs(0.3812928032069815, [5.704553583093206e-5, 0.00015212625565764492, 0.000349624141431336, 0.0007978804843345692, 0.0018271691389798897, 0.004194095333537474, 0.009639799649923194, 0.022183627514153276, 0.05116785778935707, 0.11868662341141359, 0.2793062781945737, 0.6824916959670321, 1.8418776768771283, 6.536228531781246, 56.90662030012766], [4.652078712894077e-5, 0.0003218215104770888, 0.0011049316663925537, 0.003159705062233523, 0.00849036686736275, 0.02229761527729769, 0.05806003413186596, 0.15074220148748674, 0.39131776875161906, 1.018351818022971, 2.6701839814977193, 7.144670683970085, 20.193812967613667, 66.8802591890845, 393.15622001947446], 15)

    const coneffs_12_n15 = AlgRemez_coeffs(52.801530615011785, [-6.603122250408879e-7, -4.385050932934162e-6, -2.0643603407379996e-5, -8.956905281940202e-5, -0.00037893462099888793, -0.0015892755651349718, -0.006647330899525596, -0.027809825699088738, -0.11674098552484793, -0.49492401946519393, -2.1537889443543996, -10.04004845875417, -56.31414203731785, -541.0537188024954, -47896.82148291968], [9.76789330722235e-5, 0.0004861285762563459, 0.0015447925188979874, 0.0043078364853888105, 0.011474686157071166, 0.03004933652397656, 0.07819887457469368, 0.20311391670879247, 0.5278790951523255, 1.3769550228941574, 3.6292670454376412, 9.835734487685084, 28.774699349047975, 106.22558098021766, 1108.8405448985493], 15)
    const coneffs_m12_n15 =  AlgRemez_coeffs(0.018938844922720747, [0.006232875919569935, 0.007671885966613157, 0.010882179000141433, 0.016605126881131888, 0.02616286350978532, 0.04176562692239934, 0.06703080125283097, 0.10785465346597474, 0.17392670732092233, 0.28161086985411626, 0.46047137387370596, 0.772252892886211, 1.3840926887668725, 2.9688778021493913, 11.07305624880225], [2.308717887146002e-5, 0.0002409965637633696, 0.000889670459783518, 0.0026027542764653425, 0.007053765864978664, 0.018591747424104358, 0.04849595340124775, 0.1260376463356921, 0.327370440293837, 0.8519322874091946, 2.2309978372893666, 5.942658243141148, 16.571804748421712, 52.66096512396873, 262.08312473142433], 15)

    const coeffs_14_n10 = AlgRemez_coeffs(5.7165557622098175, [-1.1472527436322903e-5, -8.22437195900441e-5, -0.0004911917249478209, -0.0028595740550169564, -0.01657546771377983, -0.09632764566721617, -0.5684080282956644, -3.574637513226843, -29.402805614736042, -910.454747964212], [0.00017585352463504747, 0.0012116676043259563, 0.005511731061457466, 0.023034355065820862, 0.09439559579370578, 0.38561054626048247, 1.5842943313284588, 6.695935386170518, 32.015219356384215, 297.1675517580863], 10)
    const coeffs_m14_n10 = AlgRemez_coeffs(0.17493050738884694, [0.000659835331355923, 0.0018359285824308825, 0.0051025759090143645, 0.014493338675151415, 0.04146036428067062, 0.1190535368426458, 0.3449286017982327, 1.0347965655277906, 3.585219290905995, 23.743072811520367], [8.614668677164346e-5, 0.0007996196969643771, 0.0038232149092825908, 0.01615861364506301, 0.06638822575850153, 0.2711990933978193, 1.1113834065181245, 4.644638810303384, 21.127906621091135, 145.57570030585526], 10)

    const coeffs_m18_n10 = AlgRemez_coeffs(2.3792132499634393, [-1.6178331673230056e-5, -9.324673262980408e-5, -0.0004616164445402523, -0.0022470840473836217, -0.01091397032788275, -0.05314702608352868, -0.26225649448725546, -1.3670670699683314, -8.966025770950466, -176.1213306168453], [0.00015069584886010342, 0.0010969629890608236, 0.005038625333107273, 0.021097078414511326, 0.0864633632847519, 0.35303065630397845, 1.4486843270648686, 6.09968682387749, 28.716462343100822, 241.67613686838126], 10)
    const coeffs_18_n10 = AlgRemez_coeffs(0.4203070069550793, [0.00012205466614334875, 0.0004400966659764438, 0.0014872525305126727, 0.005058114989269569, 0.017260926270113545, 0.05909197271336461, 0.20435471009883857, 0.7359818587930066, 3.136603644754106, 28.83640562962562], [0.0001059268835215699, 0.0008914747121053526, 0.004196936783670219, 0.017671206571184022, 0.07251494889428817, 0.29607916032240045, 1.2134381593989532, 5.080750861110906, 23.33715928002066, 169.8786011269988], 10)
    """
    f(x) = x^(y/z) = coeffs.α0 + sum_i^n coeffs.α[i]/(x + coeffs.β[i])
    f(x) = x^(-y/z) = coeffs_inverse.α0 + sum_i^n coeffs_inverse.α[i]/(x + coeffs_inverse.β[i])
    """
    struct RHMC #type for the rational Hybrid Monte Carlo
        y::Int64
        z::Int64
        coeffs::AlgRemez_coeffs
        coeffs_inverse::AlgRemez_coeffs

        function RHMC(order::Rational;n=10,lambda_low=0.0004,lambda_high=64,precision=42)
            num = numerator(order)
            den = denominator(order)
            return RHMC(num,den,n=n,lambda_low=lambda_low,lambda_high=lambda_high,precision=precision)
        end

        function RHMC(y,z;n=10,lambda_low=0.0004,lambda_high=64,precision=42)
            println("-------------------------------------------------------------")
            println("RHMC mode!")
            
            order = y // z # y/z
            num = numerator(order)
            den = denominator(order)
            @assert num != 0 "numerator should not be zero!"
            @assert num*den != 1 "$(num ÷ den) should not be 1!"
             
            if n == 10 && num == 1 && den == 2
                coeffs =coeffs_12
                coeffs_inverse =coeffs_m12
            elseif n == 15 && num == 1 && den == 2
                coeffs =coeffs_12_n15
                coeffs_inverse =coeffs_m12_n15
            elseif n == 15 && num == 1 && den == 4
                coeffs =coeffs_14
                coeffs_inverse =coeffs_m14
            elseif n == 10 && num == 1 && den == 4
                coeffs =coeffs_14_n10
                coeffs_inverse =coeffs_m14_n10
            elseif n == 15 && num == 1 && den == 8
                coeffs =coeffs_18
                coeffs_inverse =coeffs_m18
            elseif n == 10 && num == 1 && den == 8
                coeffs =coeffs_18_n10
                coeffs_inverse =coeffs_m18_n10
            elseif n == 10 && num == -1 && den == 2
                coeffs =coeffs_m12
                coeffs_inverse =coeffs_12
            elseif n == 15 && num == -1 && den == 2
                coeffs =coeffs_m12_n15
                coeffs_inverse =coeffs_12_n15
            elseif n == 15 && num == -1 && den == 4
                coeffs =coeffs_m14
                coeffs_inverse =coeffs_14
            elseif n == 10 && num == -1 && den == 4
                coeffs =coeffs_m14_n10
                coeffs_inverse =coeffs_14_n10
            elseif n == 15 && num == -1 && den == 8
                coeffs =coeffs_18
                coeffs_inverse =coeffs_18
            elseif n == 10 && num == -1 && den == 8
                coeffs =coeffs_18_n10
                coeffs_inverse =coeffs_18_n10
            else
                println("$y//$z with the order $n: coefficients for RHMC should be calculated")
                coeff_plus,coeff_minus = calc_coefficients(abs(num),den,n,lambda_low,lambda_high,precision=precision)
                if num > 0
                    coeffs =coeff_plus
                    coeffs_inverse =coeff_minus
                elseif num < 0
                    coeffs_inverse =coeff_plus
                    coeffs =coeff_minus
                end
            end

            println("the coefficients for x^{$num/$den}: ")
            display(coeffs)
            println("the coefficients for x^{-$num/$den}: ")
            display(coeffs_inverse)
            println("-------------------------------------------------------------")

            return new(num,den,coeffs,coeffs_inverse)
        end
    end

    function get_α(x::RHMC)
        return x.coeffs.α
    end

    function get_α0(x::RHMC)
        return x.coeffs.α0
    end

    function get_β(x::RHMC)
        return x.coeffs.β
    end

    function get_order(x::RHMC)
        return x.coeffs.n
    end

    function get_α_inverse(x::RHMC)
        return x.coeffs_inverse.α
    end

    function get_α0_inverse(x::RHMC)
        return x.coeffs_inverse.α0
    end

    function get_β_inverse(x::RHMC)
        return x.coeffs_inverse.β
    end



end