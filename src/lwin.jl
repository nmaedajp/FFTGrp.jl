# log-sinc^4 型の対数ウィンドウ
# f  : 周波数配列 (>0)
# f0 : 中心周波数（スムージングしたい周波数）
# b  : スムージングの強さ（20〜40 くらいがよく使われます）
# normalize=true なら、数値的に「∑ W ≈ 1」になるように正規化します。
function log_sinc4_window(f::AbstractVector,
                          f0::Real;
                          b::Real = 40.0,
                          normalize::Bool = true)
# チェック
    @assert all(>(0), f) "frequencies must be > 0"
    @assert f0 > 0 "f0 must be > 0"
#
# 対数軸での距離 t = b * log10(f / f0)
    r  = f ./ f0
    t  = b .* log10.(r)
    W = similar(f, Float64)
#
# t→0 のときは極限値 1 を使う
    for i in eachindex(f)
        ti = t[i]
        if abs(ti) < 1e-8
            W[i] = 1.0
        else
            s = sin(ti) / ti
            W[i] = s^4
        end
    end
#
    if normalize
        # 数値積分用に「その場で正規化」しておく
        s = sum(W)
        if s != 0
            W ./= s
        end
    end

    return W
end
#
#    log_smooth_spectrum(f, S; b=20.0)
#
# 対数型 log-sinc^4 ウィンドウを使って、スペクトル S(f) を
# 対数周波数軸上で平滑化する簡単な実装例。
#
function log_smooth_spectrum(f::AbstractVector,
              S::AbstractVector;
              b::Real = 40.0,
              logspec::Bool = true)
# f : 周波数配列（単調増加・>0 を想定）
# S : スペクトル（f と同じ長さ）
# b : スムージングの強さ
    N = length(f)
    @assert length(S) == N
    Ssm = similar(S, Float64)
    Ssm[1] = S[1]
    for k in 2:N
        f0 = f[k]
        # 実用上は、窓の主ローブ付近だけに絞ると速い：
        # |t| = |b*log10(f/f0)| ≤ π あたりまで使う
        logratio = b .* log10.(f ./ f0)
        mask = abs.(logratio) .<= π   # ここを変えると窓の有効幅が変わる
        fk  = f[mask]
        # if k<=10 
        #     @show k
        #     @show mask
        #     @show fk
        # end
        if logspec
            Sk  = log.(S[mask])
        else
            Sk  = S[mask]
        end
        # その近傍だけのウィンドウを計算し、その場で正規化
        Wk = log_sinc4_window(fk, f0; b=b, normalize=true)
        # 加重平均で平滑化
        if logspec
            Ssm[k] = exp.(sum(Wk .* Sk))
        else
            Ssm[k] = sum(Wk .* Sk)
        end
    end
    return Ssm
end
#
function lwin(f::AbstractVector,
              Pf::AbstractMatrix,
              nch::Int64;
              b::Real = 40.0,
              logspec::Bool = true)
#
# f  : 振動数配列．区間数と同じ．0を含んでいる．
# Pf : パワースペクトル密度．f と同じ長さ．
# b  : スムージングの強さ
# 返り値は、正の周波数成分に対応する平滑化されたパワースペクトル密度．
# 
Nyq = div(length(f), 2) + 1
Pfsm = similar(Pf, Float64, Nyq, nch)
for ich = 1:nch
    Pfsm[1:Nyq, ich] = log_smooth_spectrum(f[1:Nyq], Pf[1:Nyq, ich]; b=b, logspec = logspec)
end
return Pfsm
end
# フーリエスペクトル用の対数ウィンドウ平滑化
function lwin_sq(f::AbstractVector,
                 Xf::AbstractMatrix,
                 nch::Int64;
                 b::Real = 40.0,
                 logspec::Bool = true)
Nyq = div(length(f), 2) + 1
Xfsm = similar(Xf, Float64, Nyq, nch)
for ich = 1:nch
    Xfsm[1:Nyq, ich] = sqrt.(log_smooth_spectrum(f[1:Nyq], abs.(Xf[1:Nyq, ich].^2); b=b, logspec = logspec))
end
return Xfsm
end
