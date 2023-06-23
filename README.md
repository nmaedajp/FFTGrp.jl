# FFTGrp

## FFTを行う際の周辺関数

### ex_trend(y)
* データ y のリニアトレンドを取り除く

### dws(L; w_type="ct", base=1)
* L個のデータに対する，データウィンドウを返す．
* 初期設定は，cosine taper．
* データウィンドウを計算するとき，データの順序$j$を$-1 \leq t \leq 1$ の$t$に変換しているが，base = 1 は，配列の最初が１から始まるという意味で，$j=1$を$t=-1$に変換している．
* ウィンドウのタイプ
  * w_type = "p" || "parzen"    : parzen ウィンドウ
  * w_type = "t" || "triangle"  : 三角形ウィンドウ
  * w_type = "r" || "rectangle" : 長方形ウィンドウ
  * w_type = "b" || "boxcar"    : 長方形ウィンドウ
  * w_type = "w" || "welch"     : Welch ウィンドウ
  * w_type = "ct" || "cosine taper" : cosine taper
  * w_type = "han" || "hanning"  : ハニングウィンドウ

### swin_sq(Xf, band, N, dt)
* フーリエ変換に対するスペクトルウィンドウ
  * Xf：フーリエ変換
  * band: バンド幅
  * N：Xf のデータ数
  * dt：サンプリング周期
### swin(Pf, band, N, dt)
* パワースペクトルに対するスペクトルウィンドウ
  * Pf：パワースペクトル（two sided）
  * band: バンド幅
  * N：Xf のデータ数
  * dt：サンプリング周期

## FFTを実行するための関数

### kukan(nwave, L, M)
* (kai, k1, k2) = kukan(nwave, L, M)
* 波形データを区間分けし，区間の数，先頭と終わりの番号を返す．
  * nwave：波形データの数
  * L, M : 区間の長さとずらし幅．
  * kai：区間数
  * k1：先頭の番号，k2：終わりの番号 （整数型の配列．大きさkai）

### FFTkukan(x, L, nch, kai, k1, k2; wtype="ct")
* Ck = FFTkukan(x, L, nch, kai, k1, k2; wtype="ct")
* 戻り値は，区間ごとの FFTW.fftによる係数．
  * x : 波形データ．x(nwave,nch)
  * L : 区間長
  * kai：区間数
  * k1：先頭の番号，k2：終わりの番号 （整数型の配列．大きさkai）
  * wtype="ct"： ウィンドウのタイプ．初期設定は，cosine taper

### FourPwrAutoCo(Ck, L, nch, hz)
* (Xf, Pf, Cxx, Rxx) = FourPwrAutoCo(Ck, L, nch, hz)
  * Xf: フーリエ変換，Pf：パワースペクトル，Cxx：自己相関関数，Rxx．
  * Xf のみ複素数．他は実数．
  * 配列の大きさは，(L, nch, kai)
  
### CrsSpecCo(Ck, Cxx, L, nch, kai, hz; ich0=1)
* (Pfxy, Cxy, Rxy) = CrsSpecCo(Ck, Cxx, L, nch, kai, hz; ich0=1)
  * Pfxy：クロススペクトル，Cxy：相互相関関数，Rxy：相互相関係数
  * ich0=1：初期設定では，チャンネル１との相互相関関数，相互相関係数を計算する．

### f_tau(L, hz)
* (f, tau, tau2) = f_tau(L, hz)
* 振動数，タイムラグ（プラスのみ，正負）

### orikaesi(x, L)
* x2 = orikaesi(x, L)
  * 相互相関関数などを正負のタイムラグに対応させて，折り返す．

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nmaedajp.github.io/FFTGrp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nmaedajp.github.io/FFTGrp.jl/dev/)
[![Build Status](https://github.com/nmaedajp/FFTGrp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nmaedajp/FFTGrp.jl/actions/workflows/CI.yml?query=branch%3Amain)
