# ここで計算したいもの
# 区間分け
# フーリエ変換，パワースペクトル，自己相関関数，自己相関係数
# クロススペクトル，相互相関関数，相互相関係数．

function kukan(nwave, L, M)
    kai = floor(Int64,(nwave - L)/M + 1) 
    k1 = Array{Int64}(undef,kai)
    k2 = Array{Int64}(undef,kai)
    for ikai = 1:kai
        k1[ikai] = (ikai-1)*M + 1
        k2[ikai] = (ikai-1)*M + L
    end
    return kai, k1, k2        
end

function FFTkukan(x, L, nch, kai, k1, k2; wtype="ct")
    Ck = Array{ComplexF64}(undef, L, nch, kai)
    # 区間，チャンネルごとに，リニアトレンドの除去，ウィンドウをかけて，FFT
    for ikai = 1:kai
        k11 = k1[ikai]; k22 = k2[ikai];
        for ich = 1:nch
           y = ex_trend(x[k11:k22, ich]) .* dws(L, w_type=wtype)
           Ck[1:L, ich, ikai] = fft(y) / sqrt(w2(w_type=wtype))
        end
    end
    return Ck
end

function FourPwrAutoCo(Ck, L, nch, hz)
    Xf = Ck / hz
    Pf = real.((conj.(Ck) .* Ck)/L/hz) 
    Cxx = Array{Float64}(undef, L, nch, kai)
    Rxx = Array{Float64}(undef, L, nch, kai)
    for ikai=1:kai
        for ich=1:nch
            Cxx[:,ich,ikai] = real.(ifft(Pf[:,ich,ikai]*hz))
        end
    end
    for ikai=1:kai
        for ich=1:nch
            Rxx[:,ich,ikai] = Cxx[:,ich,ikai]/Cxx[1,ich,ikai]
        end
    end
    return Xf, Pf, Cxx, Rxx
end

function CrsSpecCo(Ck, Cxx, L, nch, kai, hz; ich0=1)
    Pfxy = Array{ComplexF64}(undef, L, nch, kai)
    Cxy  = Array{Float64}(undef, L, nch, kai)
    Rxy  = Array{Float64}(undef, L, nch, kai)
    for ich =1:nch
        Pfxy[:, ich, :] = (conj.(Ck[:,ich0,:]) .* Ck[:,ich,:])/L/hz
    end
    for ikai=1:kai
        for ich=1:nch
            Cxy[:,ich,ikai] = real.(ifft(Pfxy[:,ich,ikai]*hz))
        end
    end
    for ikai=1:kai
        for ich=1:nch
            bunbo = sqrt(Cxx[1,ich0,ikai]*Cxx[1,ich,ikai])
            Rxy[:,ich,ikai] = Cxy[:,ich,ikai]/bunbo
        end
    end
return Pfxy, Cxy, Rxy
end

function f_tau(L, hz)
    dt = 1/hz; df = hz/L
    f   = [(i-1)*df for i=1:L]
    tau = [(i-1)*dt for i=1:L]
    tau2 = [(i-L/2)*dt for i=1:L]
    return f, tau, tau2
end

function orikaesi(x, L)
    x2 = copy(x)
    L2 = floor(Int64, L/2)
    x2[L2:L]   = x[1:L2+1] 
    x2[1:L2-1] = x[L2+2:L]   
    return x2
end