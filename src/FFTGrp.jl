module FFTGrp

  using FFTW
    
  export ex_trend, dw, dws, w2, swin, swin_sq
  export kukan, FFTkukan, FourPwr, FourPwrAutoCo, CrsSpecCo, f_tau, orikaesi
  export lwin, lwin_sq

  
  include("FFTGrp11.jl")
  include("FFTGrp22.jl")
  include("lwin.jl")
  
end
