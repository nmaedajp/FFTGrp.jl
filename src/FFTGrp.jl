module FFTGrp

  using FFTW
    
  export ex_trend, dw, dws, w2, swin, swin_sq
  export kukan, FFTkukan, FourPwrAutoCo, CrsSpecCo, f_tau, orikaesi

  function ex_trend(y)
  # removing linea trend
  # x: data of horizontal axis, y: data of vertical axis
    n = length(y)                    # number of data
    x = 1.0:1.0:n
    A = hcat(x,ones(n))              # Jacobian
    ab = A\y
    yy = y - (ab[1].*x .+ab[2])      # removing trend
    avr = sum(yy)/n
    yy = yy .- avr                   # removing offset
  end
  function dws(L;w_type="ct", base=1)
    j=1:L
    dws = dw.(j, L; w_type=w_type, base=base)
    return dws
  end
  function w2(; w_type="ct")
    if w_type == "p" || w_type == "parzen"
        w2 = 151/560
    elseif  w_type == "t" || w_type == "triangle" 
        w2 = 1/3
    elseif  w_type == "r" || w_type == "rectangle" 
        w2 = 1.0
    elseif  w_type == "b" || w_type == "boxcar" 
        w2 = 1.0
    elseif w_type == "w" || w_type == "welch"
        w2 = 8/15
    elseif w_type == "ct" || w_type == "cosine taper"
        w2 = 7/8
    elseif w_type == "han" || w_type == "hanning" 
        w2 = 3/8
    else
        println("please set w_type at p, w, ct, or han.") 
    end
  end
  function j2t(j, L; base=1)
  # convert j to t
  # j: base to L-base+1,  t: -1 to 1.
  # when j=1, t=-1
    j2t =  (2/L)*(j-base)-1
  end
  function dw(j, L; w_type="ct", base=1)
    t = j2t(j,L,base=base)
    if w_type == "p" || w_type == "parzen"
        if abs(t) <= 0.5
            dw = 1.0-6.0*t^2 +6.0*abs(t)^3
        elseif abs(t) <= 1.0
            dw = 2.0*(1-abs(t))^3
        else 
            dw=0.0
        end
    elseif  w_type == "t" || w_type == "triangle" 
        if abs(t) < 1.0
            dw = 1.0-abs(t)
        else
            dw = 0.0
        end
    elseif  w_type == "r" || w_type == "rectangle" 
        if abs(t) < 1.0
            dw = 1.0
        else
            dw = 0.0
        end
    elseif  w_type == "b" || w_type == "boxcar" 
        if abs(t) < 1.0
            dw = 1.0
        else
            dw = 0.0
        end
    elseif w_type == "w" || w_type == "welch"
        if abs(t) < 1.0
            ww = 1.0-t^2
        else
            ww = 0.0
        end
    elseif w_type == "ct" || w_type == "cosine taper"
        if abs(t) < 0.8
            ww = 1.0
        else
            ww = 0.5+0.5*cos(5.0*pi*(abs(t)-0.8))
        end
    elseif w_type == "han" || w_type == "hanning" 
        if abs(t) < 1.0
            dw = (1.0 + cos(π*t))/2.0
        else
            dw = 0.0
        end
    else
        println("please set w_type at p, w, ct, or han.") 
    end
  end
  function swin_sq(Xf, band, N, dt)
    Xfs = sqrt.(swin(abs.(Xf).^2, band, N, dt))
    return Xfs
  end
  function swin(Pf, band, N, dt)
    T = N*dt; df = 1/T
    u = 151/280/band
    n = Int64(floor(band*4/df)+1)
    fw = 0.0:df:(n-1)*df
    pw  = pwin_s.(u,fw)    
#    println("band: ",band, " 2/u: ", 2/u, " df: ", df, " n: ", n)
    Pfs = zeros(N); df = 1/T
    if band == 0.0
        Pfs = Pf
    else 
        for k=1:N
            Pfs[k] = pw[1] * Pf[k] *df
            for g=2:n
                (kgm,kgp) = conidx(k, g, N)
                Pfs[k] = Pfs[k] + pw[g] *(Pf[kgm]+Pf[kgp])*df
#        println(kgm,"  ",kgp)
            end
        end
    end
    return Pfs
  end
  function pwin_s(u,f)
    Wf = (3/4)*u*dif(u*f)^4
    return Wf
  end
  function dif(x)
    t = pi*x/2
    if t == 0
        sinc = 1
    else
        sinc=sin(t)/t
    end
    return sinc
  end
  function conidx(k,g,N)
    k1 = k - 1; g1 = g-1
    kgm = k1 -g1
    kgp = k1 +g1
    n = (kgm+1, kgp+1)
    if kgm < 0
        n = (-kgm+1, kgp+1)
    end
    if kgp > N-1
        n = (kgm+1, 2N-1-kgp)
    end
    return n
  end

  include("FFTgrp22.jl")
  
end
